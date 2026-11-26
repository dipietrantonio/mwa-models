#!/usr/bin/env python3
import sys
import subprocess
import os
from collections import defaultdict
import csv

if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("Usage: process-followups.py <obsid>")
        exit(0)

    obsid = sys.argv[1]
    bash_script = f"""mkdir -p followup_trials; cd followup_trials;
scp -r "data-mover:/scratch/pawsey1154/cdipietrantonio/{obsid}_output_followup_dm*" . ;
find . -name "*.bin" -exec blink-candidates.py {{}} --csv {{}}.csv \;
find . -name "*.csv" > csv_list.txt ; 
"""
    try:
        process = subprocess.run(
            ["bash", "-c", bash_script], 
            capture_output=True, 
            text=True, 
            check=True
        )
    except subprocess.CalledProcessError as e:
        print(f"Command failed with error: {e.stderr}")
        exit(1)

    os.chdir("followup_trials")
    with open("csv_list.txt") as fp:
        dm_to_cands = defaultdict(list)
        for line in fp.readlines():
            with open(line.strip()) as csv_fp:
                dm = line.split('_')[3]
                trials = dm_to_cands[dm]
                current_trial = []
                reader = csv.reader(csv_fp, delimiter='\t')
                for row in reader:
                    if row[3] == 'x': continue # skip headers
                    current_trial.append(row[3:11])
                trials.append(current_trial)

        # now find common candidates among trials for each dm
        for dm, trials in dm_to_cands.items():
            cands_count = defaultdict(lambda : 0)
            for trial in trials:
                for x in trial:
                    cand_tuple = (int(float(x[2])), int(float(x[3])), int(x[7]))
                    cands_count[cand_tuple] += 1
            cands = [(k, v) for k, v in cands_count.items() if v > 1]
            if len(cands) > 0:
                print(f"\n====== CANDIDATE {dm} ========")
                for cand in cands:
                    ra = cand[0][0]
                    dec = cand[0][1]
                    ts = cand[0][2]
                    print(f"RA {ra} - DEC {dec} - time step {ts} - number of matches: {cand[1]}")
            
