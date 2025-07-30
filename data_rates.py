import matplotlib.pyplot as plt
from argparse import ArgumentParser
from models import MWA_PHASE_1, Correlator, Imager, Dedispersion, dispersive_delay
from plotting import make_barplot



def display_data_requirements(integration_time_s, channel_avg_factor, image_side, bits_per_pixel, n_dm_trials, repr_dm, bwcentre_mhz):
    # dispersive delay for 30MHz of bandwidth centered at 150MHz with a DM of 600
    lb_ghz = (bwcentre_mhz - 30.72/2) / 1000
    hb_ghz = (bwcentre_mhz + 30.72/2) / 1000
    delay_s = dispersive_delay(repr_dm, lb_ghz, hb_ghz) / 1000
    CORRELATOR = Correlator(MWA_PHASE_1, integration_time_s, channel_avg_factor)
    IMAGER = Imager(image_side, bits_per_pixel, CORRELATOR)
    DEDISP = Dedispersion(n_dm_trials, IMAGER)

    vcs_gibps =  MWA_PHASE_1.data_rate / 1024**3 / 8
    corr_gibps = CORRELATOR.data_rate  / 1024**3 / 8
    imaging_gibps = IMAGER.data_rate  / 1024**3 / 8
    dedispersion_gibps = DEDISP.data_date / 1024**3 / 8
    print(f"Input parameters\n================\n"
          f"Integration time (s):        {integration_time_s}\n"
          f"N. fine channels:            {3072 // 4}\n"
          f"Fine channel bandwith (kHz): {channel_avg_factor * 10}\n"
          f"Image side size (n. pixels): {image_side}\n"
          f"N. DM trials:                {n_dm_trials}\n"  
          )
    print("Data produced per second of observation (GiB) \n=============================================\n")
    print(f"VCS:           {vcs_gibps:0.2f} GiB")
    print(f"Visibilities:  {corr_gibps:0.2f} GiB")
    print(f"Images:        {imaging_gibps:0.2f} GiB")
    print(f"Dedisp images: {dedispersion_gibps:.2f} GiB\n")

    print("Data required to perform de-dispersion \n=====================================\n")

    print(f"Dispersive delay for DM = {repr_dm} pc cm-3 for 30.72 MHz of bandwidth centered at {bwcentre_mhz} MHz is {delay_s:.2f} seconds.")
    print(f"Dynamic spectra (images) total volume: {imaging_gibps * delay_s:.2f} GiB")
    print(f"DM-Arrival time (time series) total volume: {dedispersion_gibps * delay_s:.2f} GiB")
    
    print()

    # def data_rates_plot():
    #     plt.rcParams.update({'font.size': 23})
    #     title = "Significant data rates of the FRB search pipeline."

    #     make_barplot(["VCS", "Correlation", "Imaging", "Dedispersion"], [vcs_rate, corr, imaging, dedispersion], 
    #         ylabel="Output data rate (GiB/s)", title=title)
    #     plt.show()



if __name__ == "__main__":

    parser = ArgumentParser()

    parser.add_argument("--inttime", "-t", type=float, required=True, help="Integration time (in seconds).")
    parser.add_argument("--avg", "-c", required=True, type=int, help="Fine channel averaging factor.")
    parser.add_argument("--imageside", "-s", type=int, required=True, help="Image side size.")
    parser.add_argument("--bpp", "-b", type=int, default=32, help="Bits used to represent an image pixel.")
    parser.add_argument("--dmtrials", "-d", required=True, type=int, help="Number of DM trials.")
    parser.add_argument("--bwcentre", type=float, default=150, help="Centre of the 30.72MHz MWA bandwidth (in MHz).")
    parser.add_argument("--dm", type=float, default=600, help="Representative DM used to compute the dispersive delay (in pc cm-3).")
    args = vars(parser.parse_args())

    display_data_requirements(
        integration_time_s=args["inttime"],
        channel_avg_factor=args["avg"],
        image_side=args["imageside"],
        bits_per_pixel=args["bpp"],
        n_dm_trials=args["dmtrials"],
        repr_dm=args["dm"],
        bwcentre_mhz=args["bwcentre"]
    )
