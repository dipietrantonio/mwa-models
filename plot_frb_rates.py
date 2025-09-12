# see also : /home/msok/Desktop/GRANTS/2021/LIEF/SECTIONS/logbook//20220416_FRB_limits_summary.odt
# 
from optparse import OptionParser

# plotting :
# from pylab import *
import numpy as np
import numpy
import math
import matplotlib.pyplot as plt
import matplotlib.dates as md  # to convert unix time to date on X-axis 
import matplotlib.patches as patches
# from matplotlib.ticker import LogLocator
from matplotlib.ticker import FormatStrFormatter
from matplotlib.ticker import (MultipleLocator, AutoMinorLocator)
import sys
import os
import copy
import time
import math
from scipy.interpolate import interp1d # for interpolation 
from scipy import stats


import matplotlib
matplotlib.rc('xtick', labelsize=25) 
matplotlib.rc('ytick', labelsize=25) 

# import radec2azim

# time operations :
from astropy.coordinates import SkyCoord, EarthLocation
from astropy.time import Time
from datetime import datetime
MWA_POS=EarthLocation.from_geodetic(lon="116:40:14.93",lat="-26:42:11.95",height=377.8)

 
def parse_options(idx):
   global debug_level

   usage="Usage: %prog [options]\n"
   usage += "Plotting FRB rates using results of different papers scaled down to 150 or 200 MHz :\n"
   
   parser = OptionParser(usage=usage,version=1.00)

   # General parameters :
#   parser.add_option('-s','--station_name','--station',dest="station_name",default="AAVS2", help="Station name [default %default]")

   # TEST : azim=0, za=0, lst=15.4 
   parser.add_option('-g','--grid','--show_grid',action="store_true",dest="show_grid",default=False, help="Show grid lines [default %default]")
   parser.add_option('-T','--tingay','--show_tingay',action="store_true",dest="show_tingay",default=False, help="Show Tingay et al. (2015) upper limit [default %default]")
#   parser.add_option('-a','--azim','--azim_deg',dest="azim_deg",default=None, help="Pointing direction azimuth in degrees [default %default]",metavar="float",type="float")
   parser.add_option('-f','--freq','--freq_mhz',dest="freq_mhz",default=200.00, help="Survey frequency in MHz [default %default MHz]",metavar="float",type="float")
   parser.add_option('--include_alpha_minus2','--alpha_minus2',action="store_true",dest="include_alpha_minus2",default=False, help="Include curves for spectral index alpha -2 [default %default]")
   parser.add_option('--legend_with_curves','--include_curves',action="store_true",dest="legend_with_curves",default=False, help="Include curves in the legend [default %default]")
   parser.add_option('-F','--fluence_threshold','--fluence_threshold','--fluence_cutoff',dest="fluence_threshold",default=200.00, help="Fluence threshold in Jy ms [default %default Jy ms]",type="float")
   parser.add_option('--description','--add_plot_text','--text',action="store_true",dest="add_plot_text",default=False, help="Add text to plot [default %default]")
   (options, args) = parser.parse_args(sys.argv[idx:])

   print("#####################################################################################################")
   print("PARAMETERS:")
   print("#####################################################################################################")
   print("Frequency = %.2f [MHz]" % (options.freq_mhz))
   print("Include alpha = -2 : %s" % (options.include_alpha_minus2))
   print("Legend with curves = %s" % (options.legend_with_curves))
   print("Fluence threshold  = %.4f [Jy ms]" % (options.fluence_threshold))
   print("Add text to plot   = %s" % (options.add_plot_text))
   print("#####################################################################################################")
   
   return (options, args)



#####################################################################################################

#####################################################################################################
def frb_rate_euclidean( fluence, fluence_ref, frb_rate_ref , debug=True ) :
   if debug : 
      print("DEBUG : frb_rate_euclidean(%.4f,%.4f,%.4f)" % (fluence, fluence_ref, frb_rate_ref))
#   frb_rate = frb_rate_ref*math.pow( fluence/fluence_ref , -1.5 )
#   frb_rate = 1
   frb_rate = float(frb_rate_ref)*math.pow( (float(fluence)/float(fluence_ref)) , -1.5 )
   
   if debug : 
      print("DEBUG :    %.8f " % (frb_rate))
   
   return frb_rate

#####################################################################################################
# Given certain FRB rate at, for example, 700 FRBs / day / sky above fluence 700 Jy ms , it returns 
# FRB rates for various fluences in the range min_fluence to max_fluence 
#####################################################################################################
def calc_frb_rates_vs_fluence( freq_target,                    # target frequency
                               fluence_ref, frb_rate_ref ,     # reference fluence and FRB rate 
                               freq_ref=154, spectral_index=0, # reference frequency and spectral index for scaling, 0 - flat spectrum 
                               min_fluence=1, max_fluence=10000, 
                               fluence_step=1, # range to be outputted 
                               options=None, tag=""
                             ): 
   out_fluence = []
   out_rates   = []
   out_rate    = -1

   fluence_freq = fluence_ref*math.pow( (float(freq_target)/float(freq_ref)) , spectral_index )
   
   fluence = min_fluence
   while fluence < max_fluence :
      frb_rate = frb_rate_euclidean( fluence, fluence_freq, frb_rate_ref )
      
      out_fluence.append( fluence )
      out_rates.append( frb_rate )
      
      if math.fabs(fluence - options.fluence_threshold) < 0.00001 :
         out_rate = frb_rate
   
      fluence += fluence_step   

#      if fluence > 10 :
#         fluence_step = 10

#      if fluence > 100 :
#         fluence_step = 100

#      if fluence > 1000 :
#         fluence_step = 1000
         
   if out_rate > 0 :
      print("FRB rate for these data and fluence threshold = %.4f Jy ms is %.4f FRBs/day/sky (%s)" % (options.fluence_threshold,out_rate,tag))
   else :
      print("WARNING : could not determine FRB rate for the fluence threshold = %.4f Jy ms (%s)" % (options.fluence_threshold,tag))

   return (out_fluence,out_rates,out_rate)   
  
# https://matplotlib.org/3.5.0/api/_as_gen/matplotlib.pyplot.legend.html   
# legend()
# legend(handles, labels)
# legend(handles=handles)
# legend(labels)
def do_plots(options) :
   fig_size_x=20
   fig_size_y=10
   legend_list = []
   plot_list   = []
   legend_location = "upper right"

   plt.figure( figsize=( fig_size_x , fig_size_y ) )
   fig = plt.gcf()
   plt.xlim(( 0.8, 10000 ))   
   # plt.yticks( numpy.array([0.01,0.05,0.1,0.5,1.00,5,10.,100,1000,10000,100000]) )
#   plt.yticks( [0.01,0.05,0.1,0.5,1.00,5,10.,100,1000,10000,100000] )   
   axes = plt.gca()
   axes.set_yscale('log')
   plt.yticks( [0.01,0.1,10] )
#   plt.tick_params(axis='y', which='minor')
#   ax_x.yaxis.set_minor_formatter(FormatStrFormatter("%.1f"))
   # Make a plot with major ticks that are multiples of 20 and minor ticks that
   # are multiples of 5.  Label major ticks with '.0f' formatting but don't label
   # minor ticks.  The string is used directly, the `StrMethodFormatter` is
   # created automatically.
#   axes.yaxis.set_major_locator(MultipleLocator(0.1))
   # ax_x.yaxis.set_major_formatter('{x:.0f}')
   
   
   survey_threshold = []
   survey_frb_rate = []
   
   ##################################################################################################################################################################       
   #                                   Tingay et al (2015) : /home/msok//Desktop/GRANTS/2020/LIEF/SECTIONS/references/LOW_FREQ_NON-DETECTIONS/Tingay_2015_AJ_150_199_MS.pdf   
   # Antonina R. : RFRB < 82 /sky/day at 182 MHz, above a fluence of F > 7980 Jy ms (not a strong limit -> not plotted)
   fluence_limit_tingay = 700
   frb_rate_tingay = 700
   point_tingay = 'v'
   if options.show_tingay :
      (fluence_tingay,frb_rates_tingay,our_rate) = calc_frb_rates_vs_fluence( options.freq_mhz, fluence_limit_tingay, frb_rate_tingay, spectral_index=0, options=options, tag="Tingay et al" )   
      ax1 =  plt.plot( fluence_tingay,frb_rates_tingay, linestyle='dashed', color='green', linewidth=1 )   
#   ax2 =  plt.plot( fluence_tingay,frb_rates_tingay, point_tingay, color='green', linewidth=2, markersize=5 )
#   legend_list.append("Tingay et al.,2015, 154 MHz")
      ax3 = plt.plot( [fluence_limit_tingay] , [frb_rate_tingay] , point_tingay , linestyle='dashed', color='green', linewidth=2, markersize=10 )
      plot_list.append(ax3[0])   
      legend_list.append("MWA, 154 MHz, Tingay et al.,2015")
#      plt.xlim(( 1, 10000 ))
      # ax1.yaxis.set_major_locator(MultipleLocator(0.1))
   else :
      print("WARNING : upper limit from Tingay et al. 2015 not shown")
#   ax1.set_xlim(1,10000)
#   legend_list.append("Tingay et al.,2015, 154 MHz")  
   ##################################################################################################################################################################
   
   ##################################################################################################################################################################   
   #                             Shannon et al. 2018 page 5 in /home/msok/Desktop/GRANTS/2021/LIEF/SECTIONS/references/FRB/FRB_rates/ASKAP/Shannon_et_al_2018.pdf
   #
   # 37 +/- 8 /day/sky 
   # above threshold = 26 Jy ms (w/1.26 ms)^{-1/2}
   # Frequency spectral index scaling_index=-2.1
   fluence_limit_rshannon = 26
   frb_rate_rshannon = 37
   spectral_index_rshannon = -2.1
   freq_shannon = 1400.00
   point_rshannon = 's'
   
   survey_threshold.append( fluence_limit_rshannon )
   survey_frb_rate.append( frb_rate_rshannon )
   
   # add Ryan's point iteself (without any frequency scaling, i.e. spectral index = 0 )
   ax_rshannon_point = plt.plot( [fluence_limit_rshannon] , [frb_rate_rshannon] , point_rshannon , color='orange', markersize=15 )
   plot_list.append(ax_rshannon_point[0])
   legend_list.append(r'ASKAP, 1.4 GHz, Shannon et al.,2018, (measured)')

   if options.include_alpha_minus2 :
      # using Ryan Shannon et al, (2018)  spectral index from ASKAP :
      (fluence_rshannon,frb_rates_rshannon,our_rate) = calc_frb_rates_vs_fluence( options.freq_mhz, fluence_limit_rshannon, frb_rate_rshannon, spectral_index=spectral_index_rshannon, freq_ref=freq_shannon, options=options, tag="Shannon et al (2018), \\alpha=-2" )
      ax_rshannon =  plt.plot( fluence_rshannon,frb_rates_rshannon, linestyle='--', color='orange', linewidth=2, markersize=12 )
      if options.legend_with_curves :      
         plot_list.append(ax_rshannon[0])
         legend_list.append(r'ASKAP, 1.4 GHz, Shannon et al.,2018, $\alpha$=-2.1')
   
   # using spectral index = -1
   (fluence_rshannonM1,frb_rates_rshannonM1,our_rate) = calc_frb_rates_vs_fluence( options.freq_mhz, fluence_limit_rshannon, frb_rate_rshannon, spectral_index=-1 , freq_ref=freq_shannon, options=options, tag="Shannon et al (2018), \\alpha=-1" )
   ax_rshannonM1 =  plt.plot( fluence_rshannonM1,frb_rates_rshannonM1, linestyle='-.', color='orange', linewidth=2, markersize=12 )
   if options.legend_with_curves :
      plot_list.append(ax_rshannonM1[0])
      legend_list.append(r'ASKAP, 1.4 GHz, Shannon et al.,2018, $\alpha$=-1')
#   legend_list.append("Shannon et al.,2018, #alpha=-1")

   # using flat (0) spectral index :
   (fluence_rshannon0,frb_rates_rshannon0,our_rate) = calc_frb_rates_vs_fluence( options.freq_mhz, fluence_limit_rshannon, frb_rate_rshannon, spectral_index=0 , freq_ref=freq_shannon, options=options, tag="Shannon et al (2018), \\alpha=0 (use)" )
   ax_rshannon0 =  plt.plot( fluence_rshannon0,frb_rates_rshannon0, linestyle='-', color='orange', linewidth=2, markersize=12 )
   if options.legend_with_curves :
      plot_list.append(ax_rshannon0[0])
      legend_list.append(r'ASKAP, 1.4 GHz, Shannon et al.,2018, $\alpha$=0')
      
   # save to text file :
   f = open("shannon_et_al_2018.txt","w")
   for i in range(0,len(fluence_rshannon0)) :
      line = "%.6f %.6f\n" % (fluence_rshannon0[i],frb_rates_rshannon0[i])
      f.write(line)
   f.close()   

#   legend_list.append("Shannon et al.,2018, #alpha=0")
   ##################################################################################################################################################################

   ##################################################################################################################################################################
   #                                     Lovell telescop rate is similar to GBT < 5500 FRBs / day / sky at 332 MHz (not plotted) 
   # /home/msok/Desktop/GRANTS/2021/LIEF/SECTIONS/references/FRB/FRB_Newsletter_202203/2203.04890.pdf 
   ##################################################################################################################################################################

   ##################################################################################################################################################################
   #                                     Parent et al (2020) /home/msok/Desktop/GRANTS/2021/LIEF/SECTIONS/references/FRB/FRB_rates/GBT/acroread/Parent_et_al_2022.pdf
   # Rate = 3400+15400-3300 above fluence 1.44 Jy ms 
   #
   # How to convert flux limit to fluence limit ? Is it assuming Gaussian of peak = Peak flux ?
   #   
   ##################################################################################################################################################################
   fluence_limit_parent_gbt = 1.44
   frb_rate_parent_gbt = 3400
   spectral_index_parent_gbt  = 0. # 350 MHz is very close to we can assume pessimistic (specrtal index = 0) scenario
   freq_parent_gbt = 350.00
   point_parent_gbt = 'o'
   
   survey_threshold.append( fluence_limit_parent_gbt )
   survey_frb_rate.append( frb_rate_parent_gbt )


   # add GBT point itself -  (without any frequency scaling, i.e. spectral index = 0 )
   ax_parent_point = plt.plot( [fluence_limit_parent_gbt] , [frb_rate_parent_gbt] , point_parent_gbt , color='green', markersize=15 )
   plot_list.append(ax_parent_point[0])
   legend_list.append(r'GBT, 350 MHz, Parent et al., 2020 (measured)')

   if options.include_alpha_minus2 :
      (fluence_parent_gbt,frb_rates_parent_gbt,our_rate) = calc_frb_rates_vs_fluence( options.freq_mhz, fluence_limit_parent_gbt, frb_rate_parent_gbt, spectral_index=-2, freq_ref=freq_parent_gbt, options=options, tag="Parent et al (2020), \\alpha=-2" )
      ax_parent_gbtM2 =  plt.plot( fluence_parent_gbt, frb_rates_parent_gbt, linestyle='--', color='green', linewidth=2, markersize=12 )
      if options.legend_with_curves :
         plot_list.append(ax_parent_gbtM2[0])
         legend_list.append(r'GBT, 350 MHz, Parent et al., 2020, $\alpha$=-2')

   (fluence_parent_gbt,frb_rates_parent_gbt,our_rate) = calc_frb_rates_vs_fluence( options.freq_mhz, fluence_limit_parent_gbt, frb_rate_parent_gbt, spectral_index=-1, freq_ref=freq_parent_gbt, options=options, tag="Parent et al (2020), \\alpha=-1" )
   ax_parent_gbtM1 =  plt.plot( fluence_parent_gbt,frb_rates_parent_gbt, linestyle='-.', color='green', linewidth=2, markersize=12 )
   if options.legend_with_curves :
      plot_list.append(ax_parent_gbtM1[0])
      legend_list.append(r'GBT, 350 MHz, Parent et al., 2020, $\alpha$=-1')

   (fluence_parent_gbt,frb_rates_parent_gbt,our_rate) = calc_frb_rates_vs_fluence( options.freq_mhz, fluence_limit_parent_gbt, frb_rate_parent_gbt, spectral_index=0.00, freq_ref=freq_parent_gbt, options=options, tag="Parent et al (2020), \\alpha=0 (use)" )
   ax_parent_gbt =  plt.plot( fluence_parent_gbt,frb_rates_parent_gbt, linestyle='-', color='green', linewidth=2, markersize=12 )
   if options.legend_with_curves :
      plot_list.append(ax_parent_gbt[0])
      legend_list.append(r'GBT, 350 MHz, Parent et al., 2020, $\alpha$=0')   
 
   (fluence_parent_gbt_lower,frb_rates_parent_gbt_lower,our_rate) = calc_frb_rates_vs_fluence( options.freq_mhz, fluence_limit_parent_gbt, frb_rate_parent_gbt-3300, spectral_index=0.00, freq_ref=freq_parent_gbt, options=options, tag="Parent et al (2020) - lower limit, \\alpha=0" )     
   (fluence_parent_gbt_higher,frb_rates_parent_gbt_higher,our_rate) = calc_frb_rates_vs_fluence( options.freq_mhz, fluence_limit_parent_gbt, frb_rate_parent_gbt+15400, spectral_index=0.00, freq_ref=freq_parent_gbt, options=options, tag="Parent et al (2020) - upper limit, \\alpha=0" )
   plt.fill_between( fluence_parent_gbt, frb_rates_parent_gbt_lower, frb_rates_parent_gbt_higher, color='green' , alpha=0.1 )
   
   # save to text file :
   f = open("GBT_Parent_et_al_2020.txt","w")
   for i in range(0,len(fluence_parent_gbt)) :
      line = "%.6f %.6f\n" % (fluence_parent_gbt[i], frb_rates_parent_gbt[i])
      f.write(line)
   f.close()   

   f = open("GBT_Parent_et_al_2020_LOWER.txt","w")
   for i in range(0,len(fluence_parent_gbt)) :
      line = "%.6f %.6f\n" % (fluence_parent_gbt[i], frb_rates_parent_gbt_lower[i])
      f.write(line)
   f.close()   

   f = open("GBT_Parent_et_al_2020_HIGHER.txt","w")
   for i in range(0,len(fluence_parent_gbt)) :
      line = "%.6f %.6f\n" % (fluence_parent_gbt[i], frb_rates_parent_gbt_higher[i])
      f.write(line)
   f.close()   

   
#   lower_error = numpy.ones(len(fluence_parent_gbt))*3300
#   upper_error = numpy.ones(len(fluence_parent_gbt))*15400
#   asymmetric_error = [lower_error, upper_error]
#   plt.errorbar( fluence_parent_gbt, frb_rates_parent_gbt, yerr=asymmetric_error, fmt='.', color='green')
#   plt.fill_between( fluence_parent_gbt, frb_rates_parent_gbt-lower_error, frb_rates_parent_gbt+upper_error, color='green' , alpha=0.1 )
   ##################################################################################################################################################################
   
   ##################################################################################################################################################################
   # CHIME :   /home/msok/Desktop/GRANTS/2021/LIEF/SECTIONS/references/FRB/FRB_rates/CHIME 
   #             https://ui.adsabs.harvard.edu/abs/2021ApJ...923....1P/abstract
   # 
   #          In this assumption, the derived FRB rate was RFRB = 820 /sky/day above a fluence F > 5 Jy ms. in /home/msok/Desktop/GRANTS/2021/LIEF/SECTIONS/references/FRB/FRB_Newsletter_202203/2203.04890.pdf 
   #              https://ui.adsabs.harvard.edu/abs/2021ApJS..257...59C/abstract
   # " We infer a sky rate of [820 +/- 60(stat.) +220 -200(sys.)]/sky/day above a fluence of 5 Jy ms at 600 MHz, with a scattering time at 600 MHz under 10 ms and DM above 100 pc cm-3. "
   #
   ##################################################################################################################################################################
   ##################################################################################################################################################################
   fluence_limit_chime = 5
   frb_rate_chime = 820
   freq_chime = 600.00
   point_chime = 'D'
   chime_color='red'
   
   survey_threshold.append( fluence_limit_chime ) 
   survey_frb_rate.append( frb_rate_chime )

   # 
   # add Ryan's point iteself (without any frequency scaling, i.e. spectral index = 0 )
   ax_chime_point = plt.plot( [fluence_limit_chime] , [frb_rate_chime] , point_chime , color=chime_color, markersize=10 )
   plot_list.append(ax_chime_point[0])
   legend_list.append(r'CHIME, 600 MHz, Amiri et al., 2021 (measured)')

   (fluence_chime,frb_rates_chime,our_rate) = calc_frb_rates_vs_fluence( options.freq_mhz, fluence_limit_chime, frb_rate_chime, spectral_index=0, freq_ref=freq_chime, options=options, tag="CHIME, 600 MHz, Amiri et al., 2021, \\alpha=0 (use)" )
   ax_chime =  plt.plot( fluence_chime, frb_rates_chime, linestyle='-', color=chime_color, linewidth=2, markersize=12 )
   if options.legend_with_curves :
      plot_list.append(ax_chime[0])
      legend_list.append(r'CHIME, 600 MHz, Amiri et al., 2021, $\alpha$=0')
      
   # save to text file :
   f = open("chime.txt","w")
   for i in range(0,len(fluence_chime)) :
      line = "%.6f %.6f\n" % (fluence_chime[i],frb_rates_chime[i])
      f.write(line)
   f.close()   
   
   # using spectral index = -1
   (fluence_chimeM1,frb_rates_chimeM1,our_rate) = calc_frb_rates_vs_fluence( options.freq_mhz, fluence_limit_chime, frb_rate_chime, spectral_index=-1 , freq_ref=freq_chime, options=options, tag="CHIME, 600 MHz, Amiri et al., 2021, \\alpha=-1" )
   ax_chimeM1 =  plt.plot( fluence_chimeM1,frb_rates_chimeM1, linestyle='-.', color=chime_color, linewidth=2, markersize=12 )
   if options.legend_with_curves :
      plot_list.append(ax_chimeM1[0])
      legend_list.append(r'CHIME, 600 MHz, Amiri et al., 2021, $\alpha$=-1')
#   legend_list.append("Shannon et al.,2018, #alpha=-1")

   # using flat (0) spectral index :
   if options.include_alpha_minus2 :
      (fluence_chimeM2,frb_rates_chimeM2,our_rate) = calc_frb_rates_vs_fluence( options.freq_mhz, fluence_limit_chime, frb_rate_chime, spectral_index=-2, freq_ref=freq_chime, options=options, tag="CHIME, 600 MHz, Amiri et al., 2021, \\alpha=-2" )
      ax_chimeM2 =  plt.plot( fluence_chimeM2, frb_rates_chimeM2, linestyle='--', color=chime_color, linewidth=2, markersize=12 )
      if options.legend_with_curves :
         plot_list.append(ax_chimeM2[0])
         legend_list.append(r'CHIME, 600 MHz, Amiri et al., 2021, $\alpha$=-2')
   
#   lower_error = numpy.ones(len(fluence_chime))*200
#   upper_error = numpy.ones(len(fluence_chime))*220   
#   asymmetric_error = [lower_error, upper_error]
#   plt.fill_between( fluence_chime, frb_rates_chime-lower_error, frb_rates_chime+upper_error, color=chime_color , alpha=0.1 )
   (fluence_chime_lower,frb_rates_chime_lower, our_rate) = calc_frb_rates_vs_fluence( options.freq_mhz, fluence_limit_chime, frb_rate_chime-200, spectral_index=0, freq_ref=freq_chime, options=options, tag="CHIME, 600 MHz, Amiri et al., 2021, lower limit \\alpha=0" )
   (fluence_chime_higher,frb_rates_chime_higher, our_rate) = calc_frb_rates_vs_fluence( options.freq_mhz, fluence_limit_chime, frb_rate_chime+220, spectral_index=0, freq_ref=freq_chime, options=options, tag="CHIME, 600 MHz, Amiri et al., 2021, upper limit \\alpha=0" )
   plt.fill_between( fluence_chime_lower, frb_rates_chime_lower, frb_rates_chime_higher, color=chime_color , alpha=0.1 )


   ##################################################################################################################################################################
   # LOFAR :   /home/msok/Desktop/GRANTS/2021/LIEF/SECTIONS/references/FRB/FRB_rates/LOFAR/2012.08348_MS.pdf
   #           pm = Pastor-Marazuela, Ines : https://ui.adsabs.harvard.edu/abs/2021Natur.596..505P/abstract
   #           FRB rate based on FRB20180916B burst rate at 150 megahertz, we find there are 3-450 FRBs in the sky per day above 50 Jy ms.
   #           By extending this limit with Euclidean fluence scaling, the rate becomes RFRB = 90  1400 sky1 day1.
   ##################################################################################################################################################################
   fluence_limit_lofar_pm_low = 50
   frb_rate_lofar_pm_low = 3   
   fluence_limit_lofar_pm_high = 50
   frb_rate_lofar_pm_high = 450
   frb_rate_lofar_pm = (frb_rate_lofar_pm_low + frb_rate_lofar_pm_high) / 2.00
   fluence_limit_lofar_pm = (fluence_limit_lofar_pm_low+fluence_limit_lofar_pm_high)/2.00
         
   spectral_index_lofar_pm  = 0. # 150 MHz is very close to we can assume pessimistic (specrtal index = 0) scenario
   freq_lofar_pm = 150.00
   point_lofar_pm = 's'

   ax_lofar_point = plt.plot( [(fluence_limit_lofar_pm_low+fluence_limit_lofar_pm_high)/2.00] , [frb_rate_lofar_pm] , point_lofar_pm , color='blue', markersize=15 )
   plot_list.append(ax_lofar_point[0])
   legend_list.append(r'LOFAR, 150 MHz, Pastor-Marazuela et al., 2021 (measured)')

   (fluence_lofar_pm,frb_rates_lofar_pm,our_rate) = calc_frb_rates_vs_fluence( options.freq_mhz, fluence_limit_lofar_pm, frb_rate_lofar_pm, spectral_index=0, freq_ref=freq_lofar_pm, options=options, tag="LOFAR, 150 MHz, Pastor-Marazuela et al., 2021  \\alpha=0" )
   ax_parent_lofar_pm_sp0 =  plt.plot( fluence_lofar_pm,frb_rates_lofar_pm, linestyle='-', color='blue', linewidth=2, markersize=12 )
   if options.legend_with_curves :
      plot_list.append(ax_parent_lofar_pm_sp0[0])
      legend_list.append(r'LOFAR, 150 MHz, Pastor-Marazuela et al., 2021, $\alpha$=0')

   # add shaded error indicating error range :
   (fluence_lofar_low_pm,frb_rates_lofar_low_pm,our_rate) = calc_frb_rates_vs_fluence( options.freq_mhz, fluence_limit_lofar_pm_low, frb_rate_lofar_pm_low, spectral_index=0, freq_ref=freq_lofar_pm, options=options, tag="LOFAR, 150 MHz, Pastor-Marazuela et al., 2021 (lower limit)  \\alpha=0" )   
   (fluence_lofar_high_pm,frb_rates_lofar_high_pm,our_rate) = calc_frb_rates_vs_fluence( options.freq_mhz, fluence_limit_lofar_pm_high, frb_rate_lofar_pm_high, spectral_index=0, freq_ref=freq_lofar_pm, options=options, tag="LOFAR, 150 MHz, Pastor-Marazuela et al., 2021 (upper limit) \\alpha=0" )   
#   lower_error = numpy.ones(len(fluence_lofar_pm))*frb_rate_lofar_pm_low
#   upper_error = numpy.ones(len(fluence_lofar_pm))*frb_rate_lofar_pm_high
#   asymmetric_error = [lower_error, upper_error]
#   plt.fill_between( fluence_lofar_pm, lower_error, upper_error,  color='blue', alpha=0.1 )
   plt.fill_between( fluence_lofar_low_pm, frb_rates_lofar_low_pm, frb_rates_lofar_high_pm, color='blue', alpha=0.1 )

   ##################################################################################################################################################################
   # LOFAR OLD : /home/msok/Desktop/GRANTS/2021/LIEF/SECTIONS/references/FRB/FRB_Newsletter_202203/2203.04890.pdf
   # 
   #   LPPS : < 150/day/sky at F > 107 Jy ms at frequency ??? MHz ???
   #   Rawlings Array : 29/day/sky at F > 310 Jy ms at 145 MHz 
   #   FRATS : <1400/day/sky F > 6000 Jy ms  between 119 and 151 MHz    
   ##################################################################################################################################################################

   ##################################################################################################################################################################
   # UTMOST : /home/msok/Desktop/GRANTS/2021/LIEF/SECTIONS/references/FRB/FRB_Newsletter_202203/2203.04890.pdf
   #          /home/msok/Desktop/GRANTS/2021/LIEF/SECTIONS/references/FRB_rates/UTMOST
   #  1/ RFRB < 1000/sky/day at 843 MHz above a fluence F > 23 Jy ms
   #  2/ RFRB = 78 /sky/day above a fluence F > 11 Jy ms at 800 MHz ( Caleb et al, 2017 https://ui.adsabs.harvard.edu/abs/2017MNRAS.468.3746C/abstract )
   #  3/ RFRB = 98 /sky/day above a fluence F > 8 Jy ms at 800 MHz  ( Farah et al, 2019 https://ui.adsabs.harvard.edu/abs/2019MNRAS.488.2989F/abstract )
   ##################################################################################################################################################################
   fluence_limit_utmost = 8
   frb_rate_utmost = 98
   freq_utmost = 843.00
   point_utmost = 'P'
   utmost_color='pink'
   
   survey_threshold.append( fluence_limit_utmost )
   survey_frb_rate.append( frb_rate_utmost )

   # 
   # add Ryan's point iteself (without any frequency scaling, i.e. spectral index = 0 )
   ax_utmost_point = plt.plot( [fluence_limit_utmost] , [frb_rate_utmost] , point_utmost , color=utmost_color, markersize=20 )
   plot_list.append(ax_utmost_point[0])
   legend_list.append(r'UTMOST, 843 MHz. Farah et al. (2019), (measured)')
   
   (fluence_utmost,frb_rates_utmost,our_rate) = calc_frb_rates_vs_fluence( options.freq_mhz, fluence_limit_utmost, frb_rate_utmost, spectral_index=0, freq_ref=freq_utmost, options=options, tag="UTMOST, 843 MHz. Farah et al. (2019) \\alpha=0 (use)" )
   ax_utmost =  plt.plot( fluence_utmost, frb_rates_utmost, linestyle='-', color=utmost_color, linewidth=2, markersize=12 )
   if options.legend_with_curves :
      plot_list.append(ax_utmost[0])
      legend_list.append(r'UTMOST, 843 MHz, Farah et al. (2019), $\alpha$=0')
   
   # using spectral index = -1
   (fluence_utmostM1,frb_rates_utmostM1,our_rate) = calc_frb_rates_vs_fluence( options.freq_mhz, fluence_limit_utmost, frb_rate_utmost, spectral_index=-1 , freq_ref=freq_utmost, options=options, tag="UTMOST, 843 MHz. Farah et al. (2019) \\alpha=-1" )
   ax_utmostM1 =  plt.plot( fluence_utmostM1,frb_rates_utmostM1, linestyle='-.', color=utmost_color, linewidth=2, markersize=12 )
   if options.legend_with_curves :
      plot_list.append(ax_utmostM1[0])
      legend_list.append(r'UTMOST, 843 MHz, Farah et al. (2019), $\alpha$=-1')
#   legend_list.append("Shannon et al.,2018, #alpha=-1")

   # using flat (0) spectral index :
   if options.include_alpha_minus2 :
      (fluence_utmostM2,frb_rates_utmostM2,our_rate) = calc_frb_rates_vs_fluence( options.freq_mhz, fluence_limit_utmost, frb_rate_utmost, spectral_index=-2, freq_ref=freq_utmost, options=options, tag="UTMOST, 843 MHz. Farah et al. (2019) \\alpha=-2" )
      ax_utmostM2 =  plt.plot( fluence_utmostM2, frb_rates_utmostM2, linestyle='--', color=utmost_color, linewidth=2, markersize=12 )
      if options.legend_with_curves :
         plot_list.append(ax_utmostM2[0])
         legend_list.append(r'UTMOST, 843 MHz, Farah et al. (2019), $\alpha$=-2')

   # save to text file :
   f = open("UTMOST_Farash_et_al_2019.txt","w")
   for i in range(0,len(fluence_utmost)) :
      line = "%.6f %.6f\n" % (fluence_utmost[i], frb_rates_utmost[i])
      f.write(line)
   f.close()   




   ##################################################################################################################################################################
   # PARKES : see /home/msok/Desktop/GRANTS/2021/LIEF/SECTIONS/logbook//20220416_FRB_limits_summary.odt
   #   /home/msok/Desktop/GRANTS/2021/LIEF/SECTIONS/references/FRB/Parkes/
   #   Shivani Bhendari : https://ui.adsabs.harvard.edu/abs/2018MNRAS.475.1427B/abstract
   #   Not good : https://ui.adsabs.harvard.edu/abs/2019MNRAS.488.2989F/abstract
   ##################################################################################################################################################################
   fluence_limit_parkes = 2
   frb_rate_parkes = 1.7*1000
   freq_parkes = 1350.00 # 1200 - 1500 MHz 
   point_parkes = 'd'
   parkes_color='magenta'
   
   survey_threshold.append( fluence_limit_parkes )
   survey_frb_rate.append( frb_rate_parkes )

   
   # https://ui.adsabs.harvard.edu/abs/2019MNRAS.488.2989F/abstract 81 FRBs X. Yang 
   # fluence_limit_parkes = 0.33
   # frb_rate_parkes 

   # 
   # add Ryan's point iteself (without any frequency scaling, i.e. spectral index = 0 )
   ax_parkes_point = plt.plot( [fluence_limit_parkes] , [frb_rate_parkes] , point_parkes , color=parkes_color, markersize=15 )
   plot_list.append(ax_parkes_point[0])
   legend_list.append(r'Parkes, 1350 MHz. Bhandari et al. (2018), (measured)')
   
   (fluence_parkes,frb_rates_parkes,our_rate) = calc_frb_rates_vs_fluence( options.freq_mhz, fluence_limit_parkes, frb_rate_parkes, spectral_index=0, freq_ref=freq_parkes, options=options, tag="Parkes, 1350 MHz. Bhandari et al. (2018), \\alpha=0 (use)" )
   ax_parkes =  plt.plot( fluence_parkes, frb_rates_parkes, linestyle='-', color=parkes_color, linewidth=2, markersize=12 )
   if options.legend_with_curves :
      plot_list.append(ax_parkes[0])
      legend_list.append(r'Parkes, 1350 MHz. Bhandari et al. (2018), $\alpha$=0')
#   ax_parkes.yaxis.set_major_locator(LogLocator(base=1))      
   
   # using spectral index = -1
   (fluence_parkesM1,frb_rates_parkesM1,our_rate) = calc_frb_rates_vs_fluence( options.freq_mhz, fluence_limit_parkes, frb_rate_parkes, spectral_index=-1 , freq_ref=freq_parkes, options=options, tag="Parkes, 1350 MHz. Bhandari et al. (2018), \\alpha=-1" )
   ax_parkesM1 =  plt.plot( fluence_parkesM1,frb_rates_parkesM1, linestyle='-.', color=parkes_color, linewidth=2, markersize=12 )
   if options.legend_with_curves :
      plot_list.append(ax_parkesM1[0])
      legend_list.append(r'Parkes, 1350 MHz. Bhandari et al. (2018), $\alpha$=-1')
#   legend_list.append("Shannon et al.,2018, #alpha=-1")

   # using flat (0) spectral index :
   if options.include_alpha_minus2 :
      (fluence_parkesM2,frb_rates_parkesM2,our_rate) = calc_frb_rates_vs_fluence( options.freq_mhz, fluence_limit_parkes, frb_rate_parkes, spectral_index=-2, freq_ref=freq_parkes, options=options, tag="Parkes, 1350 MHz. Bhandari et al. (2018), \\alpha=-2" )
      ax_parkesM2 =  plt.plot( fluence_parkesM2, frb_rates_parkesM2, linestyle='--', color=parkes_color, linewidth=2, markersize=12 )
      if options.legend_with_curves :
         plot_list.append(ax_parkesM2[0])
         legend_list.append(r'Parkes, 1350 MHz. Bhandari et al. (2018), $\alpha$=-2')

   # save to text file :
   f = open("Parkes_Bhandari_et_al_2018.txt","w")
   for i in range(0,len(fluence_parkes)) :
      line = "%.6f %.6f\n" % (fluence_parkes[i], frb_rates_parkes[i])
      f.write(line)
   f.close()   

   # Add Ian's plot here :
   # read data
   # plt.plot
   # etc 


#   plt.legend( legend_list, loc=legend_location, fontsize=20)
   fontsize=10
   legend_ncol=2 
   if not options.legend_with_curves :
      legend_ncol=1
   plt.legend( plot_list, legend_list, loc=legend_location, fontsize=fontsize, ncol=legend_ncol) # 14

   # mark region : https://matplotlib.org/3.1.1/api/_as_gen/matplotlib.pyplot.html
   # plt.axhspan(ymin, ymax[, xmin, xmax])   
#   plt.axhspan( 0.0001, 200, 90, 1000, alpha=0.3, facecolor='0.5' )
#   plt.axvspan( 200, 10000, ymin=0.235, ymax=0.535, alpha=0.5, facecolor='red' )
   # legend of the span :
# OLD AT TOP  plt.axvspan( 1, 30, ymin=0.94, ymax=1.1, alpha=0.5, facecolor='red' )
   if options.add_plot_text :
      fontsize=20
      start_region=0.9
      plt.axvspan( start_region, 1400, ymin=0.01, ymax=0.22, alpha=0.5, facecolor='red' )
#   esc = (r' Red colour marks the region with fluences (F) $\ge$100 Jy ms.') # Red marks the region above a threshold of 100 Jy ms,
      desc = (r' The existing measurements (data points) and their errors (shaded regions)') # Red marks the region above a threshold of 100 Jy ms,
      plt.text(start_region,0.005,desc,fontsize=fontsize)
      desc = (r' scaled to %d MHz using spectral index $\alpha=$0 (solid lines) and $\alpha=$-1 (dashed ' % (options.freq_mhz))
      plt.text(start_region,0.002,desc,fontsize=fontsize)
      desc = (r' dotted lines). The red colour marks the region with fluences (F) $>$100 Jy ms')
      plt.text(start_region,0.0008,desc,fontsize=fontsize)
      desc = (r' where the FRB rate is between 0.2 and 180 per day ($\sim$360 - 65000 per year)')
      plt.text(start_region,0.0003,desc,fontsize=fontsize)
      desc = (r' and decreases according to $\propto F^{-3/2}$ scaling for the Euclidean Universe.')
      plt.text(start_region,0.0001,desc,fontsize=fontsize)
# OLD AT TOP  plt.text(1,3500000,desc,fontsize=15)
#   desc = (r'measurements, $\alpha$ <= -1 and Euclidean Universe')
# OLD AT TOP   plt.text(1,1700000,desc,fontsize=15)
 
#   points = [[100, 2], [100,600], [10000, 0.6], [10000, 0.002] ]
   points = [[100, 1], [100,170], [9000, 0.18], [9000, 0.001] ]
   polygon= plt.Polygon(points,  fill=True, color='red' , alpha=0.5, edgecolor='r', figure=fig)
   # plt.add_patch(polygon)
   plt.gca().add_patch(polygon)


   # fit line to data points :
   # need to log10 data points 
   # stats.linregress(x_arr,y_arr)
   x_arr = numpy.log10( survey_threshold )
   y_arr = numpy.log10( survey_frb_rate )
#   slope,intercept, rval, pval, stderr = stats.linregress(x_arr,y_arr)
   coeffs = np.polyfit(x_arr,y_arr,deg=1)
   slope = coeffs[0]
   intercept = coeffs[1]
   print("DEBUG FIT : slope = %.4f , intercept = %.4f" % (slope,intercept))
   max_fluence = 5000
   x_arr_plot = numpy.arange( 1 , max_fluence+1 )
   y_arr_plot = numpy.ones( max_fluence , dtype=numpy.float64 ) # dtype=
   y_arr_plot_log = numpy.ones( max_fluence , dtype=numpy.float64 )

   
   i=0
   c = math.pow(10.00,intercept)
   print("DEBUG_FIT c = %.8f" % (c))
   for x in x_arr_plot :
      y_log = slope*math.log10(x) + intercept
      y_arr_plot_log[i] = y_log
      # y_arr_plot[i] = math.pow( 10.00, y_log )
      p = math.pow(x,slope)
      y_arr_plot[i] = c*p
      # y_log = math.log10(y_arr_plot[i])
      if int(x) % 10 == 0 :
         frbs_per_day = y_arr_plot[i]
         frbs_per_year = frbs_per_day*365.00
         print("DEBUG : %.4f %.4f -> %.8f %.8f (from %.8f*%.8f) -> %.8f FRBs / day -> %.4f FRBs / year" % (x,math.log10(x),y_arr_plot[i],y_log,c,p,frbs_per_day,frbs_per_year))
      i += 1
   ax_fit =  plt.plot( x_arr_plot, y_arr_plot, linestyle='-', color='black', linewidth=2 )
#   matplotlib.pyplot.loglog( numpy.log10( x_arr_plot ), y_arr_plot_log )
   
   plt.xlabel('Fluence [Jy ms]' , fontsize=30 )
   plt.ylabel('#FRBs / day / sky' , fontsize=30 )
   plt.xscale('log')
   plt.yscale('log')
      
   if options.show_grid :
      plt.grid()
   plt.ylim(( 0.00006,  1e6 ))      
   y_axis_ticks  = numpy.array([0.01,0.05,0.1,0.3,0.5,1.00,5,10.,100,1000,10000,100000])
   y_axis_labels = [ '{}'.format(i) for i in y_axis_ticks ]
   plt.yticks( y_axis_ticks , y_axis_labels , fontsize=15 )
   # ax_fit.yaxis.set_major_locator(LogLocator(base=1))
   
   plt.savefig( "frb_rates.png" )
   plt.show()   



if __name__ == "__main__":
    (options, args) = parse_options(1)
    
    do_plots(options)
    
        
