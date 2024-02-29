import matplotlib.pyplot as plt
from models import MWA_PHASE_1, Correlator, Imager
from plotting import make_barplot

CORRELATOR = Correlator(MWA_PHASE_1, 0.05, 4)
IMAGER = Imager(CORRELATOR)


#########################################################################
#                      DATA RATE STUDY
#########################################################################
vcs_rate =  MWA_PHASE_1.data_rate / 1024**3 / 8
corr = CORRELATOR.data_rate  / 1024**3 / 8
imaging = IMAGER.data_rate  / 1024**3 / 8

def data_rates_plot():
    plt.rcParams.update({'font.size': 23})
    title = "Significant data rates of the FRB search pipeline."

    make_barplot(["VCS", "Correlation", "Imaging"], [vcs_rate, corr, imaging], 
        ylabel="Output data rate (GiB/s)", title=title)
    plt.show()



if __name__ == "__main__":
    data_rates_plot()