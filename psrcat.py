#!/usr/bin/env python3

from __future__ import annotations

import io
import re
from typing import Any, Dict, Iterable, Optional, Tuple, Union
import requests
import pandas as pd
from bs4 import BeautifulSoup
import argparse


ATNF_PROC_FORM_URL = "https://www.atnf.csiro.au/research/pulsar/psrcat/proc_form.php"


_REF_BRACKET_RE = re.compile(r"【\d+†[^】]+】")
_MULTI_WS_RE = re.compile(r"[ \t]+")


def query_atnf_psrcat(
    query_params: Dict[str, Any],
    *,
    base_url: str = ATNF_PROC_FORM_URL,
    timeout: Union[float, Tuple[float, float]] = (10.0, 60.0),
    user_agent: str = "python-requests (contact: you@example.com)",
) -> pd.DataFrame:
    """
    Query the ATNF pulsar catalogue proc_form.php endpoint and parse the returned table.

    Parameters
    ----------
    query_params:
        Dict of query parameters to send to proc_form.php (these become the URL's query string).
        Example keys (as in your URL): version, Name, JName, DM, startUserDefined, sort_attr,
        sort_order, condition, coords_unit, radius, coords_1, coords_2, pulsar_names,
        ephemeris, style, no_value, fsize, state, ...
    base_url:
        Endpoint URL (default is ATNF proc_form.php).
    timeout:
        Requests timeout, either a float or (connect_timeout, read_timeout).
    user_agent:
        Set a UA; some sites are friendlier with one.
    return_raw:
        If True, also return the raw response text used for parsing.

    Returns
    -------
    pandas.DataFrame (or (DataFrame, raw_text) if return_raw=True)

    Parsing strategy
    ---------------
    1) If the response *looks like* CSV (semicolon-delimited header), parse with pandas.read_csv.
       (ATNF "short csv ..." uses semicolons per the docs.)
    2) Otherwise, parse the plain-text table extracted from HTML.

    Raises
    ------
    requests.HTTPError on HTTP errors
    ValueError if no table could be parsed
    """
    headers = {"User-Agent": user_agent}

    resp = requests.get(base_url, params=query_params, headers=headers, timeout=timeout)
    resp.raise_for_status()

    # The endpoint typically returns HTML even for "csv" styles; table content is inside it.
    html = resp.text
    soup = BeautifulSoup(html, "html.parser").pre
    text = soup.text
    # Try CSV first (semicolon-separated is typical for ATNF "short csv ..." output).
    df = _parse_plaintext_table(text)
    return df



def _parse_plaintext_table(extracted_text: str) -> pd.DataFrame:
    """
    Parse the text table rendered in the HTML output.

    Looks for:
        dashed rule
        header line starting with '#'
        dashed rule
        numbered rows
        dashed rule

    Cleans ATNF's bracketed reference markers like: 
    """
    # Clean bracketed reference markers and normalize whitespace
    cleaned = _REF_BRACKET_RE.sub("", extracted_text).replace("\xa0", " ")
    lines = [ln.rstrip() for ln in cleaned.splitlines() if ln.strip()]

    # Locate a header line (starts with '#') and dashed separators nearby
    header_i = None
    for i, ln in enumerate(lines):
        if ln.lstrip().startswith("#"):
            header_i = i
            break
    if header_i is None:
        raise ValueError("Could not find a header line (starting with '#') in the response.")

    # Column names: split header after '#'
    header_line = lines[header_i].lstrip()[1:].strip()
    colnames = [c for c in _MULTI_WS_RE.split(header_line) if c]

    # Data usually begins after the next dashed rule following the header
    def is_dashed(s: str) -> bool:
        s2 = s.strip()
        return len(s2) >= 5 and set(s2) <= {"-", " "}

    start = None
    for j in range(header_i + 1, len(lines)):
        if is_dashed(lines[j]):
            start = j + 1
            break

    # Read rows until the next dashed rule or footer.
    rows = []
    for ln in lines[start:]:
        if is_dashed(ln):
            break
        if ln.strip().upper().startswith("PLEASE READ"):
            break

        # Normalize whitespace
        ln2 = _MULTI_WS_RE.sub(" ", ln.strip())

        # Expect leading row index like "1 ..."
        parts = ln2.split(" ")
        if not parts or not parts[0].isdigit():
            continue

        # Drop the leading row number
        parts = parts[1:]

        # Heuristic: Some styles include +/- error columns. If we have more fields than colnames,
        # keep the leftmost len(colnames) fields (common for link-duplication artifacts).
        # if len(parts) >= len(colnames):
        #     parts = parts[: len(colnames)]
        #     rows.append(parts)
        # Name PSRJ DM
        rows.append([parts[0], parts[2], parts[4]])

    return pd.DataFrame(rows, columns=colnames)


# ------------------- Example usage -------------------
if __name__ == "__main__":

    parser = argparse.ArgumentParser()

    parser.add_argument("--radius", type=float, default=1.0, help="Radius (in degrees) or the search area defined by the coordinates.")
    parser.add_argument("COORDINATES", nargs="+", help="Whitespace-delimited oordinates, in degrees, defining the centre of the search area.")

    args = vars(parser.parse_args())

    params = {
        "version": "2.7.0",
        "Name": "Name",
        "JName": "JName",
        "DM": "DM",
        "startUserDefined": "true",
        "sort_attr": "jname",
        "sort_order": "asc",
        "condition": "",
        "coords_unit": "rajd/decjd",
        "radius": str(args["radius"]),
        "coords_1": str(args["COORDINATES"][0]),
        "coords_2": str(args["COORDINATES"][1]),
        "pulsar_names": "",
        "ephemeris": "short",
        "style": "long with last digit error",
        "no_value": "*",
        "fsize": "3",
        "x_axis": "",
        "x_scale": "linear",
        "y_axis": "",
        "y_scale": "linear",
        "state": "query",
    }

    df = query_atnf_psrcat(params)

    print("==============================")
    print("Number of results:", len(df))
    print("==============================")
    if len(df) > 0:
        print(df)
