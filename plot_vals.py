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

import pandas
import seaborn


def read_time_series(filename):
	with open(filename, "rb") as fp:
		while True:
			header = fp.read(12)
			if len(header) == 0: break
			side_size, n_dms, buffer_size = struct.unpack("iii", header)
			data_size = side_size * side_size * n_dms * buffer_size
			content = fp.read(data_size * 4)
			table = np.array(array.array('f', content), dtype=float) 
			table = table.reshape((side_size, side_size, n_dms, buffer_size))
			yield table


def read_candidates(filename):
	with open(filename, "rb") as fp:
		while True:
			header = fp.read(20)
			if len(header) == 0: break
			global_offset, side_size, n_dms, buffer_size, n_candidates = struct.unpack("iiiii", header)
			header_info = (global_offset, side_size, n_dms, buffer_size, n_candidates)
			candidates = []
			for i in range(n_candidates):
				data = fp.read(28)
				candidate = struct.unpack("iififff", data) # x, y, dm, time_step, value, avg, std
				candidates.append(candidate)
			yield (header_info, candidates)


def read_candidates_series(filename):
	with open(filename, "rb") as fp:
		timestep_offset = 0
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
				candidate[5] += timestep_offset
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
			timestep_offset += buffer_size
			batch_id += 1


def plot_candidates(candidates, header_info):
	df = pandas.DataFrame(sorted(candidates, key = lambda x: x[6]),  columns=['ID', 'Batch ID', 'x', 'y', 'dm', 'Time step', 'Peak Flux (Jy)', 'Mean Flux (Jy)', 'std'])
	df.sort_values(by='Peak Flux (Jy)')
	# pandas.set_option('display.max_rows', None)
	print(df)
	seaborn.scatterplot(df, x='x', y='y', hue='dm', size='Peak Flux (Jy)', sizes=(5, 200), palette='muted')
	ax = plt.gca()
	ax.set_xlim([0, header_info[1]])
	ax.set_ylim([0, header_info[1]])
	plt.title("Location and flux of all candidates present in the observation.") 
	plt.show()


def process_candidates_test(filename, snr, cand_ids = None):
	res = list(read_candidates_series(filename))
	# filter by SNR
	candidate_set = []
	if snr > -1:
		for x in res: candidate_set.extend([y for y in x[1] if ((y[6] - y[7])/y[8]) >= snr])
	else:
		for x in res: candidate_set.extend(x[1])
	header_info = res[0][0]
	plot_candidates(candidate_set, header_info)

	if cand_ids is not None:
		print(cand_ids)
		for candidate in candidate_set:
			if candidate[0] in cand_ids:
				# plot time series
				batch_id = candidate[1]
				series_idxs = res[batch_id][2]
				series = res[batch_id][3]
				tmp_coord = candidate[2:5]
				tmp_coord[2] = int(tmp_coord[2] - 50)
				i = series_idxs[tuple(tmp_coord)]
				plt.plot(series[i, :])
	plt.show()

def process_candidates(filename, snr):
	res = list(read_candidates(filename))
	#for time_step in range(header_info[3]):
	#	current_candidates = [c for c in candidate_set if c[3] == time_step]
	candidate_set = []
	if snr > -1:
		for x in res: candidate_set.extend([y for y in x[1] if ((y[4] - y[5])/y[6]) >= snr])
	else:
		for x in res: candidate_set.extend(x[1])
	plot_candidates(candidate_set, header_info = res[0][0])


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


def peak_finding(filename, alpha, SNR):
	results = []
	t_b = 0
	for table in read_time_series(filename):
		N_x, N_y, n_dms, n_steps = table.shape
		X = list(range(N_x))
		Y = list(range(N_y))
		for dm in range(n_dms):
			#for i, (x, y) in enumerate(product(X, Y)):
			x, y = 0, 0 # 58, 630
			# result = running_mean(table[x, y, dm, :], alpha, SNR)
			result = simple_avg(table[x, y, dm, :], SNR)
			result['dm'] = dm
			result['coord'] = (x, y)
			results.append(result)
		# t_b += n_steps
	return results


def plot_result(result):

	plt.plot(result['raw'], label='raw')
	plt.plot(result['avg'], label='avg')
	plt.scatter([result['candidates'][i][0] for i in range(len(result['candidates']))],
			[result['candidates'][i][1] for i in range(len(result['candidates']))])

	plt.legend()
	plt.show()




def compare_files(filename1, filename2):
	table1 = next(read_time_series(filename1)).flatten()
	table2 = next(read_time_series(filename2)).flatten()
	for i in range(table1.shape[0]):
		print(table1[i], table2[i])
		if i == 10: break


def get_time_series(filename):
	ts = []
	for table in read_time_series(filename):
		ts.extend(table[58, 630, 1, :])
	return ts



def plot_pixels(filename, dm_list, pixel_coord):
	arr = [list() for i in range(len(dm_list) * len(pixel_coord))]
	legend = [f"DM Idx = {dm}, coord = {c}" for dm in dm_list for c in pixel_coord]	
	for table in read_time_series(filename):
		for j, dm in enumerate(dm_list):
			for i, (x, y) in enumerate(pixel_coord):
				arr[j * len(pixel_coord) + i].extend(table[x, y, dm, :])
	
	for i in range(len(arr)): plt.plot([i * 0.02 for i in range(len(arr[i]))], arr[i])
	# plt.legend(legend)
	plt.title("Time series for selected candidate.")
	plt.xlabel("Time (s)")
	plt.ylabel("Flux density (Jy)")
	plt.show()


def plot_text_files(files):
	for txtfile in files:
		vals = [float(x) for x in open(txtfile).readlines()] 
		plt.plot(vals)
	plt.legend(files)
	plt.show()

def parse_coordinates(coord_string):
    matches = re.findall(r'\(\s*([-+]?\d*\.?\d+)\s*,\s*([-+]?\d*\.?\d+)\s*\)', coord_string)
    coordinates = [(int(x), int(y)) for x, y in matches]
    
    return coordinates


if __name__ == "__main__":
	parser = argparse.ArgumentParser()
	parser.add_argument('--txt', action='store_true')
	parser.add_argument('--dm', type=str)
	parser.add_argument('--pixels', type=str )
	parser.add_argument('--alpha', type=float, default=0.9)
	parser.add_argument('--snr', type=float, default=-1)
	parser.add_argument('--plot', action='store_true')
	parser.add_argument('--peakfind', action='store_true')
	parser.add_argument('files', type=str, nargs='+')
	parser.add_argument('--candidates', action='store_true')
	parser.add_argument('--candidatesnew', action='store_true')
	parser.add_argument('--ids', type=str)
	
	
 	
	args = vars(parser.parse_args())
	if args['peakfind']:
		for res in peak_finding(args['files'][0], args['alpha'], args['snr']):
			plot_result(res)
		
	elif args['plot']:	
		dm_idxs = [int(v) for v in args['dm'].split(',')]	
		coords = parse_coordinates(args['pixels'])
		plot_pixels(args['files'][0], dm_idxs, coords) 
	elif args['txt']:
		plot_text_files(args['files'])
	elif args['candidates']:
		process_candidates(args['files'][0], args['snr'])
	elif args['candidatesnew']:
		candidates_ids = None if args['ids'] is None else [int(x) for x in args['ids'].split(',')]
		process_candidates_test(args['files'][0], args['snr'], candidates_ids)
	else:
		compare_files(*args['files'])

