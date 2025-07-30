from asvo import parse_asvo_results_xml
import matplotlib.pyplot as plt
from collections import defaultdict
from bisect import bisect_left

def freq_dist(results):
    centre_freqs = [x["center_frequency_mhz"] for x in results]
    plt.hist(centre_freqs)
    plt.title("Centre frequency distribution across MWAX VCS observations")
    plt.ylabel("Number of Observations")
    plt.xlabel("Frequency (MHz)")
    plt.show()



def proj_to_freq(results):
    proj_to_freq = defaultdict(lambda : defaultdict(lambda : 0))
    for x in results:
        pid = x["projectid"]
        frq = x["center_frequency_mhz"]
        cdict = proj_to_freq[pid]
        cdict[frq] = cdict[frq] + 1
    
    for x, v in proj_to_freq.items():
        freq_to_count_str = " - ".join(f"Frequency: {f},  count: {c}" for f, c in sorted(v.items(), key=lambda y: y[0]))
        print(f"Project code: {x:5}, {freq_to_count_str}")


def freq_to_project(results):
    freq_to_proj = defaultdict(lambda : defaultdict(lambda : 0))
    for x in results:
        pid = x["projectid"]
        frq = x["center_frequency_mhz"]
        cdict = freq_to_proj[frq]
        cdict[pid] = cdict[pid] + 1
    
    for x, v in sorted(freq_to_proj.items(), key= lambda f: f[0]):
        freq_to_count_str = " - ".join(f"Project: {p},  count: {c}" for p, c in sorted(v.items(), key=lambda y: y[0]))
        print(f"Central frequency: {x:5}, {freq_to_count_str}")

def filter_project_obs(proj_id, min_dur):
    total = 0

    list_of_obs_to_save = set()
    for x in results:
        pid = x["projectid"]
        duration = x["duration"]
        frq = x["center_frequency_mhz"]
        if  pid != proj_id or frq < 180 or duration < 1200: continue
        # if x["calibration"] == 0: continue
        obsid = x["obs_id"]
        conf = x["mwa_array_configuration"]
        # if not "Extended" in conf: continue
        bytes = x["total_archived_data_bytes"] / 1000**4
        total += bytes
        cal_id = find_cal(obsid, cal_obs, frq)
        if cal_id == None: continue
        if cal_id['center_frequency_mhz'] != frq: raise Exception
        list_of_obs_to_save.add(obsid)
        list_of_obs_to_save.add(cal_id['obs_id'])
        print(f"{obsid:10} {x['starttime_utc']} {conf:20} {frq:3.2f} {duration:5} {bytes:5.2f} {cal_id['obs_id']}")

        
    print("Total space: ", total)
    print("List of observations to save: ", sorted(list_of_obs_to_save))

def find_cal(obsid, cals, frq):
    idx = bisect_left(cals, obsid, key=lambda x: x['obs_id'])
    while idx < len (cals) and cals[idx]["center_frequency_mhz"] != frq:
        idx += 1
    if idx == len(cals): return None
    return cals[idx]



if __name__ == "__main__":
    # results = parse_asvo_results_xml("data/mwaxvcs-2013-2025.xml") 
    results = parse_asvo_results_xml("data/voltagestart-2013-2020.xml") + parse_asvo_results_xml("data/voltagestart-2020-2025.xml")
    cal_obs = parse_asvo_results_xml("data/cal_obs.xml")

    filter_project_obs('G0024', 1000)
    # freq_to_project(results)


