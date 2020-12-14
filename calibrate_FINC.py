# ============================================================
#
# FINC Calibration
# code for comparing LAUDA bath temperature sensor and thermocouple
# array placed in 9 wells in a pcr tray
#
# Author: Andrew Loeppky
# NBD Group, ETH Zurich, Dec 11/2020
#
# Make sure to change default temp controller format .XLS to .txt
# My favorite format for most data is dictionaries full of np arrays
#
# ============================================================

# %%
import numpy as np
from matplotlib import pyplot as plt
from pathlib import Path
import matplotlib.pyplot as plt
from scipy.stats import linregress
import datetime

# %%
def count_lines(txt_file):
    """
    counts the number of lines in a given txt file
    input the folder path, returns int
    """
    count = -1  # include 0 as index for creating temp and time arrays
    for line in open(txt_file):
        count += 1
    return count


def get_finc_ramp(file_path):
    """
    Input: ramp.txt from FINCrunner.m
    specify absolute file path

    Out: dictionary containing datetime and lauda temp
    """
    finc_time = np.empty(count_lines(file_path), dtype="datetime64[us]")
    finc_temp = np.empty(count_lines(file_path))
    i = 0

    for line in open(file_path):
        try:
            the_str = line.split(",")
            # get datetime and reformat to np.datetime
            # https://stackabuse.com/converting-strings-to-datetime-in-python/
            dt = datetime.datetime.strptime(the_str[0], "%d-%b-%Y %H:%M:%S")
            finc_time[i] = np.datetime64(dt)

            # get the temp
            finc_temp[i] = the_str[1]
            i += 1
        except ValueError:
            # omit metadata lines
            pass

    # create and return the dictionary full of np arrays
    out_dict = {}
    out_dict["finc_time"] = finc_time
    out_dict["finc_temp"] = finc_temp
    return out_dict


def get_tc_ramp(file_path):
    """
    Input: Save temp controller file as .txt (tab delimited) and change
    nothing else. Specify absolute path to where file is saved.

    Output: 12 thermocouple time series and a datetime list
    as a dictionary.
    """

    tc_time = np.empty(count_lines(file_path), dtype="datetime64[us]")
    tc1 = np.empty(count_lines(file_path))
    tc2 = np.empty_like(tc1)
    i = 0
    omitted_lines = 0

    for line in open(file_path):
        # get the time
        try:
            the_str = line.split(",")
            raw_date = the_str[1] + " " + the_str[2]
            dt = datetime.datetime.strptime(raw_date, "%Y-%m-%d %H:%M:%S")
            tc_time[i] = np.datetime64(dt)
            print(tc_time[i])
            # get 12 thermocouple temperatures
            tc1[i] = the_str[3]
            tc2[i] = the_str[4]
            i += 1
        
        except ValueError:
            # omit lines containing metadata
            omitted_lines += 1

    print(f'Omitted {omitted_lines} lines')

    
        


# finc data
finc_ramp_path = "C:/Users/Owner/UBC_F2020/FINC/outputs/20201209_calibrate/ramp.txt"
# thermocouple data
tc_ramp_path = "C:/Users/Owner/UBC_F2020/FINC/outputs/20201209_calibrate/TMB01001.csv"

get_tc_ramp(tc_ramp_path)

# a = get_finc_ramp(finc_ramp_path)


# %%
