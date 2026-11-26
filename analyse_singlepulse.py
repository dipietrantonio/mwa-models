#!/usr/bin/env python3
import matplotlib.pyplot as plt
import sys
import argparse
import os

def filekey(filename : str):
	parts = filename.split('_')
	return parts[1] + parts[2] + f"{int(parts[3][:-12]):02}"


def parse_singlepulse_file(filename):
	with open(filename, "r") as fd:
		lines = fd.readlines()
	lines = lines[1:] # skip header
	cands = []
	for line in lines:
		v = [float(x) for x in line.split()]
		cands.append((v[2], v[1]))
	return cands

if __name__ == "__main__":
	parser = argparse.ArgumentParser()
	parser.add_argument("singlepulse_files", nargs='+', type=str, help=".singlepulse file containing PRESTO candidates.")
	args = vars(parser.parse_args())

	time_s = []
	SNRs = []
	dm0_time_s = []
	dm0_SNRs = []

	offset_s = 0
	for file in sorted(args["singlepulse_files"], key=filekey):
		cands = parse_singlepulse_file(file)
		for c in cands:
			time_s.append(c[0] + offset_s)
			SNRs.append(c[1])
			if c[1] > 6.5:
				print(file, c)
		
		dm0_name = "dm0/" + file[:-12] + "_dm0.singlepulse"
		if os.path.exists(dm0_name) and False:
			cands_dm0 = parse_singlepulse_file(dm0_name)
			for c in cands_dm0:
				dm0_time_s.append(c[0] + offset_s)
				dm0_SNRs.append(c[1])
		if len(cands):
			offset_s += cands[-1][0]

	plt.scatter(time_s, SNRs)
	plt.scatter(dm0_time_s, dm0_SNRs, c='r')
	plt.title("Single pulse candidates")
	plt.xlabel("Time (s)")
	plt.ylabel("Sigma")
	plt.show()

