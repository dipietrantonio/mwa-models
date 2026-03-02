#!/usr/bin/env python3

"""
Query the ATNF PSRCAT 'proc_form.php' endpoint.

Examples:
  python atnf_psrcat_query.py --names b0525+21 --names J0437-4715 --out result.txt
  python atnf_psrcat_query.py --names b0525+21 --params Name,JName,DM,P0 --style "long with last digit error"
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path
from typing import Dict, List, Optional
from urllib.parse import urljoin
import requests

def query_psrcat(psr_name) -> dict:
    ENDPOINT = "https://www.atnf.csiro.au/research/pulsar/psrcat/proc_form.php"
    payload: Dict[str, str] = {
        "version": "2.7.0",
        "startUserDefined": "true",
        "sort_attr": "jname",
        "sort_order": "asc",
        "condition": "",
        "coords_unit": "raj/decj",
        "radius": "",
        "coords_1": "",
        "coords_2": "",
        "pulsar_names": f"{psr_name}\n",
        "ephemeris": "long",
        "ephemeris_submit": "",
        "style": "long with last digit error",
        "no_value": "*",
        "fsize": "3",
        "x_axis": "",
        "x_scale": "linear",
        "y_axis": "",
        "y_scale": "linear",
        "state": "query",
    }

    # Use POST to avoid URL-length issues and to mimic HTML form submission.
    r = requests.get(ENDPOINT, payload)
    r.raise_for_status()
    content = r.text
    istart = content.index("<pre>")
    iend = content.index("</pre>")
    content = content[istart + 6: iend]

    def parse_val(v):
        try:
            return int(v)
        except Exception:
            try:
                return float(v)
            except Exception:
                return v

    # Parse ASCII table in a dictionary
    data = {}
    for line in content.splitlines():
        cols = line.split()
        data[cols[0]] = parse_val(cols[1])
    return data


def estimate_peak_flux_jy(mean_flux, period_s, width_ms):
    """
    mean_flux = Pulse_fluence / Period = Width * Peak_flux / Period 

    Then, peak_flux = mean_flux * (Period / Width)
    """
    peak_mJy = mean_flux * period_s * 1000 / width_ms
    return peak_mJy / 1000



if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Query ATNF PSRCAT proc_form.php")
    parser.add_argument("--psr", type=str, required=True, help="Name of the reference pulsar.")
    parser.add_argument("--duration", type=float, default=None, help="Duration of the observation, in seconds.")
    parser.add_argument("--detected", type=int, default=None, help="Number of detected pulses.")
   
    args = vars(parser.parse_args())
   
    try:
        table_data = query_psrcat(args["psr"])
    except requests.RequestException as e:
        print(f"ERROR: request failed: {e}", file=sys.stderr)
        exit(1)
    
    peak_jy = estimate_peak_flux_jy(table_data["S150"], table_data["P0"], table_data["W10"])
    
    if args["duration"] is None or args["detected"] is None:
        print(f"Pulsar has a peak flux of {peak_jy:.4f} Jy")
    else:
        n_expected_pulses = args["duration"] / table_data["P0"]
        gain = args["detected"] / n_expected_pulses * 100
        print(f"Number of detected pulses is {gain}")

        # TODO: we need to add our detection threshold in Jy. We indeed are not supposed
        # to detect all the pulses, but just the ones above the SNR threshold. So the
        # question is: of all the pulses the pulsar emitted above the SNR threshold,
        # how many of those we detect? This will tell us how far our real threshold is
        # from the predicted/expected threshold based on the measured noise.
        # But is it the measured noise or the model/predicted noise.
        # Also, to get the expected number of pulses above threshold we need to know
        # the distribution!!! Can we use beamforming to get that distribution? Single
        # pulse searches might be hard even with beamforming.

        # TODO: do not use ATNF values!! USE SMART!!!!
