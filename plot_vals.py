#!/usr/bin/env python3
from sys import argv
from matplotlib import pyplot as plt
import os
import struct
import numpy as np
import array
import argparse
import re
from statistics import mean, stdev
from itertools import product
from matplotlib.animation import FuncAnimation
import pandas
import seaborn



def read_candidates_series(filename):
	with open(filename, "rb") as fp:
		candidate_id = 0
		batch_id = 0
		while True:
			header = fp.read(24)
			if len(header) == 0: break
			global_offset, side_size, n_dms, buffer_size, n_candidates, n_series = struct.unpack("iiiiii", header)
			header_info = (global_offset, side_size, n_dms, buffer_size, n_series)
			candidates = []
			for i in range(n_candidates):
				data = fp.read(28)
				candidate = [candidate_id, batch_id] + list(struct.unpack("iififff", data)) # ID, x, y, dm, time_step, value, avg, std
				candidate[5] += global_offset
				candidates.append(candidate)
				candidate_id += 1
			series = np.zeros((n_series, buffer_size))
			series_idxs = {}
			for i in range(n_series):
				series_coord = struct.unpack("iif", fp.read(12))
				series_idxs[series_coord] = i
				series_data = fp.read(4 * buffer_size)
				series[i, :] = np.array(array.array('f', series_data), dtype=float)
			yield (header_info, candidates, series_idxs, series)
			batch_id += 1


def plot_candidates(candidates, header_info):
	df = pandas.DataFrame(sorted(candidates, key = lambda x: x[6]),  columns=['ID', 'Batch ID', 'x', 'y', 'dm', 'Time step', 'Peak Flux (Jy)', 'Mean Flux (Jy)', 'std'])
	df.sort_values(by='Peak Flux (Jy)')
	pandas.set_option('display.max_rows', None)
	print(df)
	seaborn.scatterplot(df, x='x', y='y', hue='dm', size='Peak Flux (Jy)', sizes=(5, 200), palette='muted')
	ax = plt.gca()
	ax.set_xlim([0, header_info[1]])
	ax.set_ylim([0, header_info[1]])
	plt.title("Location and flux of all candidates present in the observation.") 
	plt.show()


def create_movie(filename, snr, cand_ids = None):
	res = list(read_candidates_series(filename))
	header_info = res[0][0]
	candidates = []
	if snr > -1:
		for x in res: candidates.extend([y for y in x[1] if ((y[6] - y[7])/y[8]) >= snr])
	else:
		for x in res: candidates.extend(x[1])
	df = pandas.DataFrame(sorted(candidates, key = lambda x: x[6]),  columns=['ID', 'Batch ID', 'x', 'y', 'dm', 'Time step', 'Peak Flux (Jy)', 'Mean Flux (Jy)', 'std'])
	sort = df.sort_values(by='Time step')
	figure = plt.figure()
	ax = plt.gca()
	ax.set_xlim([0, header_info[1]])
	ax.set_ylim([0, header_info[1]])
	plt.title("Location and flux of all candidates present in the observation.") 
	
	def update(frame):
		frame_at = sort[sort['Time step'] == frame]
		# ax.clear()
		ax.clear()
		ax.set_xlim([0, header_info[1]])
		ax.set_ylim([0, header_info[1]])
		scat = seaborn.scatterplot(frame_at, x='x', y='y', hue='dm', size='Peak Flux (Jy)', sizes=(5, 200), palette='muted')
		
	ani = FuncAnimation(figure, update, frames=200, interval=500)
	ani.save('my_animation.mp4', writer='ffmpeg', fps=5)
	plt.show()
	

def process_candidates(filename, snr, plot_candidates = True, cand_ids = None):
	res = list(read_candidates_series(filename))
	# filter by SNR
	candidate_set = []
	if snr > -1:
		for x in res: candidate_set.extend([y for y in x[1] if ((y[6] - y[7])/y[8]) >= snr])
	else:
		for x in res: candidate_set.extend(x[1])
	header_info = res[0][0]
	if plot_candidates:
		plot_candidates(candidate_set, header_info)
	time_res = 0.02
	if cand_ids is not None:
		legend = []
		for candidate in candidate_set:
			if candidate[0] in cand_ids:
				# plot time series
				batch_id = candidate[1]
				start_step = res[batch_id][0][0]
				series_length = res[batch_id][0][3]
				series_idxs = res[batch_id][2]
				(x, y, dm) = candidate[2:5]
				i = series_idxs[(x, y, dm)]
				time_axis = [time_res * i for i in range(start_step, start_step + series_length)]
				series = res[batch_id][3]
				plt.title("Time series")
				plt.ylabel("Flux (Jy)")
				plt.xlabel("Time (s)")
				legend.append(f"Pixel ({x}, {y}), DM = {dm}")
				plt.plot(time_axis, series[i, :])
		plt.legend(legend)
	plt.show()


def running_mean(time_series, alpha, SNR):
	running_mean = time_series[0]
	result = {'raw' : [running_mean], 'avg' : [running_mean], 'candidates' : []}
	for t, val in enumerate(time_series[1:], 1):
		val = abs(val)
		result['raw'].append(val)
		if val / running_mean >= SNR:
			result['candidates'].append((t, val))
		running_mean = (alpha) * running_mean + (1 -alpha) * val
		result['avg'].append(running_mean)
	return result

def simple_avg(time_series, SNR):
	avg = abs(mean(time_series))
	sd = stdev(time_series)
	result = {'raw' : [], 'avg' : [avg] * len(time_series), 'candidates' : []}
	for t, val in enumerate(time_series, 0):
		result['raw'].append(val)
		if (val - avg) / sd >= SNR:
			result['candidates'].append((t, val))
	return result



def parse_coordinates(coord_string):
    matches = re.findall(r'\(\s*([-+]?\d*\.?\d+)\s*,\s*([-+]?\d*\.?\d+)\s*\)', coord_string)
    coordinates = [(int(x), int(y)) for x, y in matches]
    
    return coordinates


if __name__ == "__main__":
	parser = argparse.ArgumentParser()
	parser.add_argument('--txt', action='store_true')
	parser.add_argument('--dm', type=str)
	parser.add_argument('--pixels', type=str)
	parser.add_argument('--alpha', type=float, default=0.9)
	parser.add_argument('--snr', type=float, default=-1)
	parser.add_argument('--plot', action='store_true')
	parser.add_argument('--peakfind', action='store_true')
	parser.add_argument('files', type=str, nargs='+')
	parser.add_argument('--candidates', action='store_true')
	parser.add_argument('--candidatesnew', action='store_true')
	parser.add_argument('--movie', action='store_true')
	parser.add_argument('--ids', type=str)
	
	args = vars(parser.parse_args())

	if args['candidates']:
		candidates_ids = None if args['ids'] is None else [int(x) for x in args['ids'].split(',')]
		process_candidates(args['files'][0], args['snr'], candidates_ids)
	elif args['movie']:
		create_movie(args['files'][0], args['snr'])
	else:
		print("No action specified.")

