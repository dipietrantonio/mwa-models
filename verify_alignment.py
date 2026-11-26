import numpy as np
import struct
from statistics import stdev, median, mean
import os
import shutil


def read_filterbank(filename):

    string_keys = {
        "rawdatafile", "source_name"
    }

    int_keys = {
        "machine_id", "telescope_id", "data_type",
        "barycentric", "pulsarcentric", "nbits",
        "nsamples", "nchans", "nifs", "nbeams",
        "ibeam"
    }

    double_keys = {
        "az_start", "za_start", "src_raj", "src_dej",
        "tstart", "tsamp", "fch1", "foff",
        "refdm", "period"
    }

    def read_string(f):
        length = struct.unpack("<i", f.read(4))[0]

        if not 0 < length < 65536:
            raise ValueError("Invalid header string length")

        return f.read(length).decode("ascii")

    header = {}

    with open(filename, "rb") as f:

        if read_string(f) != "HEADER_START":
            raise ValueError("Not a SIGPROC filterbank")

        while True:

            key = read_string(f)

            if key == "HEADER_END":
                break

            if key in string_keys:
                value = read_string(f)

            elif key in int_keys:
                value = struct.unpack("<i", f.read(4))[0]

            elif key in double_keys:
                value = struct.unpack("<d", f.read(8))[0]

            elif key == "signed":
                value = struct.unpack("<b", f.read(1))[0]

            else:
                raise ValueError(f"Unknown header key: {key}")

            header[key] = value

        offset = f.tell()

    nchans = header["nchans"]
    nbits = header["nbits"]
    nifs = header.get("nifs", 1)
    signed = header.get("signed", 0)

    dtypes = {
        8: "i1" if signed else "u1",
        16: "<i2" if signed else "<u2",
        32: "<f4"
    }

    if nbits not in dtypes:
        raise ValueError(f"Unsupported nbits: {nbits}")

    dtype = np.dtype(dtypes[nbits])

    bytes_per_sample = nchans * nifs * dtype.itemsize
    payload_size = os.path.getsize(filename) - offset

    nsamples, remainder = divmod(
        payload_size, bytes_per_sample
    )

    if remainder:
        raise ValueError("Incomplete spectrum at end of file")


    header["nsamples"] = nsamples
    header["header_size"] = offset
    return header


if __name__ == "__main__":

    import sys
    for filename in sys.argv[1:]:
        header = read_filterbank(filename)
        for k,v in header.items():
            print(k, v)

