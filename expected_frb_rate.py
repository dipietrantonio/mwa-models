from basics import expected_frb_daily_rate
from matplotlib import pyplot as plt
plt.rcParams.update({'font.size': 20})
labels = ["CHIME", "ASKAP", "UTMOST"]
references = [ {'ref_fluence' : 5, 'ref_rate': 820, 'ref_freq' : 600e6},
     {'ref_fluence' : 26, 'ref_rate': 37, 'ref_freq' : 1400e6},
     {'ref_fluence' : 8, 'ref_rate': 98, 'ref_freq' : 843e6}     
]
for ref in references:
    F = [100 + i for i in range(100)]
    R = [expected_frb_daily_rate(ref['ref_rate'], ref['ref_freq'], ref['ref_fluence'], 150e6, f, 0) for f in F]
    SMART_R = [r/2 * (1.5/24) for r in R]
    print(sum(SMART_R))
    plt.plot(F, SMART_R)


plt.title("Expected number of FRBs in SMART derived from FRB rates measured \nby other telecopes. (⍺ = 0)")
plt.xlabel("Fluence (Jy ms)")
plt.ylabel("# FRBs")
plt.legend(labels)
plt.show()