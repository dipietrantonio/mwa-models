import xml.etree.ElementTree as ET
from matplotlib import pyplot as plt
from collections import defaultdict
from statistics import mode

plt.rcParams.update({'font.size': 20})

def __parse_value(val : str):
    if val is None: return None
    if val.isdecimal():
        return int(val)
    else:
        try:
            return float(val)
        except Exception as e:
            return val


def parse_asvo_results_xml(filename):
    tree = ET.parse(filename)
    root = tree.getroot()
    table = root[0][0]
    header = [table[i].attrib["name"] for i in range(len(table) -1)]
    data = table[-1][0]

    results = list()
    for row in data:
        results.append(dict(zip(header, 
                (__parse_value(row[c].text) for c in range(len(row))))))
    return results


def print_general_stats(obs):
    #no_cal_obs 
    print(len(obs))
    print("Total duration (hours) of SMART VCS observations:", 
          sum(x['duration'] for x in obs) / 3600 )

    print("Total space (PiB) occupied by SMART VCS:",
          sum(x['total_archived_bytes'] for x in obs) / 1024**5)


def histogram_per_year(obs):
    years = []
    year_to_projects = defaultdict(list)
    year_to_label = dict()

    def get_project_desc(pid):
        for o in obs:
            if o['projectid'] == pid:
                return o['projectshortname'], o['projectdescription']

    for o in obs:
        year = int(o["starttime_utc"].split("-")[0])
        years.append(year)
        year_to_projects[year].append(o["projectid"])

    for year in sorted(year_to_projects.keys()):
        projects = year_to_projects[year]
        pmode = mode(projects)
        name, desc = get_project_desc(pmode)
        print(f"Year {year} - Top project: {pmode} ({projects.count(pmode)}) - {name}")
    
    plt.hist(years)
    plt.title("Number of VCS observations per year")
    plt.ylabel("Number of observations")
    plt.xlabel("Year")
    plt.show()

def pie_chart_hour_project(obs):
    project_to_hours = defaultdict(lambda : 0)
    total = 0
    for o in obs:
        project_to_hours[o['projectid']] += o['duration']
        total += o['duration']
    
    for p in project_to_hours:
        project_to_hours[p] /= total
    
    labels = []
    values = []
    for p, v in project_to_hours.items():
        labels.append(p)
        values.append(v)
    
    plt.title("Hours of observation (%) per project")
    plt.pie(values, labels=labels)
    plt.show()

def histogram_category(obs, cat):

    cat_to_count = defaultdict(lambda : 0)
    for o in obs:
        cat_to_count[o[cat]] += 1
    
    x, y = [], []
    for a, b in cat_to_count.items():
        x.append(a)
        y.append(b)
    
    plt.bar(x, y)
    plt.ylabel("Number of observations")
    plt.xlabel(cat)
    plt.show()

def print_all_smart(obs):
    for x in obs:
        if x['projectid'] == 'G0057': print(x['obs_id'])


if __name__ == "__main__":
    # do not get cal observations
    obs = [x for x in parse_asvo_results_xml("data/all_vcs.xml") if not x["calibration"] and not x["deleted_flag"]]
    print_all_smart(obs)
    exit(0)
    for x in obs:
        if x['center_frequency_mhz'] < 120: print(x['obs_id'])

    # print_general_stats(obs)
    # pie_chart_hour_project(obs)
    histogram_category(obs, "mode")
    histogram_category(obs, "center_frequency_mhz")