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

    Out: time, temp as nparrays
    """
    # this is horrible, please change to less memory abusive method ASAP
    finc_time = np.array([],dtype="datetime64[us]")
    finc_temp = np.array([])

    for line in open(file_path):
        try:
            the_str = line.split(",")
            # get datetime and reformat to np.datetime
            # https://stackabuse.com/converting-strings-to-datetime-in-python/
            dt = datetime.datetime.strptime(the_str[0], "%d-%b-%Y %H:%M:%S")
            finc_time = np.append(finc_time, np.datetime64(dt))

            # get the temp
            finc_temp = np.append(finc_temp, the_str[1])
        except ValueError:
            # omit metadata lines
            pass

    return finc_time, finc_temp


def get_tc_ramp(file_path):
    """
    Input: Save temp controller file as .csv and change
    nothing else. Specify absolute path to where file is saved.

    Output: datetime nparray, 12 thermocouple time series and a as a dictionary.
    """

    # initialize time and thermocouple temps
    tc_time = np.array([], dtype="datetime64[us]")
    tc1 = np.array([])
    tc2 = np.array([])
    tc3 = np.array([])
    tc4 = np.array([])
    tc5 = np.array([])
    tc6 = np.array([])
    tc7 = np.array([])
    tc8 = np.array([])
    tc9 = np.array([])
    tc10 = np.array([])
    tc11 = np.array([])
    tc12 = np.array([])


    

    i = 0
    for line in open(file_path):
        # get the time
        try:
            the_str = line.split(",")
            raw_date = the_str[1] + " " + the_str[2]
            dt = datetime.datetime.strptime(raw_date, "%Y-%m-%d %H:%M:%S")
            tc_time = np.append(tc_time, np.datetime64(dt))

            # get 12 thermocouple temperatures
            tc1 = np.append(tc1,float(the_str[3]))
            tc2 = np.append(tc2,float(the_str[5]))
            tc3 = np.append(tc3,float(the_str[7]))
            tc4 = np.append(tc4,float(the_str[9]))
            tc5 = np.append(tc5,float(the_str[11]))
            tc6 = np.append(tc6,float(the_str[13]))
            tc7 = np.append(tc7,float(the_str[15]))
            tc8 = np.append(tc8,float(the_str[17]))
            tc9 = np.append(tc9,float(the_str[19]))
            tc10 = np.append(tc10,float(the_str[21]))
            tc11 = np.append(tc11, float(the_str[23]))
            tc12 = np.append(tc12, float(the_str[25]))

            i += 1

        except ValueError:
            # omit lines containing metadata
            pass

    # package results in dictionary
    tc_dict = {
        "tc1": tc1,
        "tc2": tc2,
        "tc3": tc3,
        "tc4": tc4,
        "tc5": tc5,
        "tc6": tc6,
        "tc7": tc7,
        "tc8": tc8,
        "tc9": tc9,
        "tc10": tc10,
        "tc11": tc11,
        "tc12": tc12,
    }
    return tc_time, tc_dict


# this will become main()
# ====================================================================================

# finc data
finc_ramp_path = "C:/Users/Owner/UBC_F2020/FINC/outputs/20201209_calibrate/ramp.txt"
# thermocouple data
tc_ramp_path = (
    "C:/Users/Owner/UBC_F2020/FINC/outputs/20201209_calibrate/TMB01001_fixed.csv"
)



finc_time, finc_temp = get_finc_ramp(finc_ramp_path)
cal_time, cal_temp = get_tc_ramp(tc_ramp_path)

fig, ax = plt.subplots(figsize=(10,10))
ax.plot(finc_time, finc_temp, label='lauda')
#plt.gca().invert_yaxis()
for key in cal_temp:
    ax.plot(cal_time, cal_temp[key], label=key)
plt.legend()


# %%
