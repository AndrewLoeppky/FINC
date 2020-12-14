# ============================================================
#
# FINC Calibration
# code for comparing LAUDA bath temperature sensor and thermocouple
# array placed in 9 wells in a pcr tray
#
# Author: Andrew Loeppky
# NBD Group, ETH Zurich, Dec 11/2020
#
# Make sure to change default temp controller format .XLS to .csv
#
# ============================================================

# %%
import numpy as np
from matplotlib import pyplot as plt
from pathlib import Path
import pandas as pd
import matplotlib.pyplot as plt
from scipy.stats import linregress
import datetime

# %%
def count_lines(txt_file):
    """
    counts the number of lines in a given txt file
    input the folder path, returns int
    """
    count = -1 # include 0 as index for creating temp and time arrays
    for line in open(txt_file):
        count += 1
    return count 


def get_finc_ramp(file_path):
    '''
    Input: ramp.txt from FINCrunner.m
    specify absolute file path

    Out: dictionary containing datetime and lauda temp
    '''
    finc_time = np.empty(count_lines(file_path),dtype='datetime64[us]')
    finc_temp = np.empty(count_lines(file_path))
    i = 0

    for line in open(file_path):
        try:
            the_str = line.split(',')
            # get datetime and reformat to np.datetime
            # https://stackabuse.com/converting-strings-to-datetime-in-python/
            dt = datetime.datetime.strptime(the_str[0], '%d-%b-%Y %H:%M:%S')
            finc_time[i] = np.datetime64(dt)
            #npdt = np.datetime64(dt)
            #finc_time[i] = npdt
            
            #get the temp
            finc_temp[i] = the_str[1]
            i += 1
        except ValueError:
            # omit metadata lines
            pass
        
    
    plt.plot(finc_time,finc_temp)
    


def get_tc_ramp(file_path):
    '''
    Input: Save temp controller file as .csv and change
    nothing else. Specify absolute path to where file is saved.
    
    Output: 12 thermocouple time series and a datetime list
    as a dictionary.
    '''
    tc_panda = pd.read_csv(file_path)
    tc_dict = {}
    tc_dict['datetime'] = tc_panda['Date'] + ' ' + tc_panda['Time']
    print(tc_dict)



# finc data
finc_ramp_path = "C:/Users/Owner/UBC_F2020/FINC/outputs/20201209_calibrate/ramp.txt"
# thermocouple data
tc_ramp_path = "C:/Users/Owner/UBC_F2020/FINC/outputs/20201209_calibrate/TMB01001.csv"

#get_tc_ramp(tc_ramp_path)

get_finc_ramp(finc_ramp_path)



# %%
