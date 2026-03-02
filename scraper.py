#!/usr/bin/env python3

import requests
import json
import argparse
import pandas as pd

def query_scraper(ra_deg, dec_deg, radius = 1):
    response = requests.get(
        f"https://pulsar.cgca-hub.org/api?type=search&ra={ra_deg}&dec={dec_deg}&radius={radius}"
    )
    data = response.json()
    results = pd.DataFrame(columns=["Name", "DM", "RA", "DEC", "Distance", "Survey"])
    for j in data:
        if j.startswith("search") or j == "nmatches": continue
        fields = data[j]
        row = (j, fields["dm"]["value"], fields["ra"]["value"], fields["dec"]["value"], fields["distance"]["value"], fields["survey"]["value"])
        results.loc[len(results)] = row
    return results.sort_values(by="Distance")


if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument("--radius", type=float, default=1.0, help="Radius (in degrees) or the search area defined by the coordinates.")
    parser.add_argument("--dm", type=float, default=-1, help="Filter results by the specified DM +/- 2.")
    parser.add_argument("COORDINATES", nargs="+", help="Whitespace-delimited oordinates, in degrees, defining the centre of the search area.")

    args = vars(parser.parse_args())
    print(query_scraper(args["COORDINATES"][0], args["COORDINATES"][1], args["radius"]))