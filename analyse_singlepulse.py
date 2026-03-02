#!/usr/bin/env python3
import matplotlib.pyplot as plt
import sys

if __name__ == "__main__":
	if len(sys.argv) < 2:
		print(f"Usage: {sys.argv[0]} <singlepulse> file.")
		exit(0)
	with open(sys.argv[1], "r") as fd:
		lines = fd.readlines()
	
	lines = lines[1:] # skip header
	cands = []
	for line in lines:
		v = [float(x) for x in line.split()]
		cands.append((v[2], v[1]))
	x = [v[0] for v in cands]
	y = [v[1] for v in cands]
	plt.scatter(x, y)
	plt.title("Single pulse candidates")
	plt.xlabel("Time (s)")
	plt.ylabel("Sigma")
	plt.show()

