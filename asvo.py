import xml.etree.ElementTree as ET

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

if __name__ == "__main__":
    smart_obs = parse_asvo_results_xml("data/asvo-results-smart.xml")
    
    print("Total duration (hours) of SMART VCS observations:", 
          sum(x['duration'] for x in smart_obs) / 3600 )
