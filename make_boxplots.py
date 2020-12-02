"""
===============================================================
DATA PROCESSING LIBRARY - FINC OUTPUTS

Author: Andrew Loeppky
UBC/ETH Atmospheric Science, NBD Group
Last Update Dec 1/2020

***Currently under development, don't trust this code yet***

functions for processing raw output from Freezing Ice Nuclei 
Counter (FINC) and generating plots.

TODO:
modify functions to return plot objects instead of just 
printing the plots

modify the remaining functions to accept the new dictionary 
format

create a prototype differential concentration plot fcn
===============================================================
"""
# %%
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import scipy.io as sio
import csv

# %%
def get_mat(the_file):
    """
    retrieves freezing temperature vector from matlab workspace
    and returns it as a 288x1 numpy array

    'the_file' (str) is the name of the experiment folder
    """
    the_path = "C:/Users/Owner/UBC_F2020/FINC/outputs/"
    load_this = the_path + the_file + "/workspace.mat"
    mat = sio.loadmat(load_this)
    frzTall = mat["FrzTall"]

    return frzTall



# %%
def get_csv(the_file):
    """
    NOT WORKING ATM

    output a freezing temperature dictionary exptname:nparray()
    from a csv. CSVs are assumed to be in the format:

    metadata1, ..., metadata2, ...
    name1,     name2, ...,   nameN
    frztemp,   frztemp, ..., frztemp
    .
    .
    .
    As per ETH repository entries found at:
    https://www.research-collection.ethz.ch/handle/20.500.11850/438875
    """
    reader = csv.reader(open(the_file))
    print(reader)
    """
    result = {}
    for row in reader:
        key = row[1]
        if key in result:
            # implement your duplicate row handling here
            pass
        result[key] = row[1:]
    return result
    """


# %%
def make_hist(data):
    """
    Returns a histogram of freezing events at temp T given a dictionary
    of samplename:([freezetemp1,freezetemp2, etc])

    WARNING: NOT NORMALIZED TO N(T)!!! Strongly dependent on
    droplet size
    """
    fig, ax = plt.subplots()
    for key in data.keys():
        ax.hist(data[key], bins=40, label=key, alpha=0.6)

    ax.legend()
    ax.set_ylabel("Freezing Temp ($^oC$)")
    ax.set_title("Freezing Temp Distrubutions")
    ax.set_ylabel("$N_{wells}$")


# %%
def make_boxplot(data):
    """
    given one or more input freezing vectors, make boxplots of
    freezing temperature
    """
    fig, ax = plt.subplots()
    ax.boxplot(data.values(), showmeans=True, whis=[10, 90])
    ax.set_xticklabels(data.keys(), rotation=-60)
    ax.set_ylabel("Freezing Temp ($^oC$)")
    ax.set_title("Freezing Temp Distrubutions")


# %%
def make_ff_curve(data):
    """
    creates a frozen fraction (FF) curve as a function of
    temperature
    """
    fig, ax = plt.subplots()

    for key in data.keys():
        length = len(data[key])
        ind = np.linspace(length, 1, length) / length
        dat = np.sort(data[key], axis=0)
        ax.plot(dat, ind, label=key)

    ax.legend()
    ax.set_xlabel("Temp ($^oC$)")
    ax.set_ylabel("FF")
    ax.set_title("Frozen Fraction")


# %%
def make_big_K(data, norm=1, Vwell=10):
    """
    creates a cumulative concentration curve from Vali eq.

    default normalization constant (eg surface area) = 1
    default well volume = 10 microlitre
    """
    # make frozen fraction curve and calculate K from it
    fig, ax = plt.subplots()

    for key in data.keys():
        length = len(data[key])
        ind = np.linspace(length, 1, length) / length
        dat = np.sort(data[key], axis=0)

        K = -np.log(1 - ind) / (norm * Vwell * 10 ** -6)  # Vali 2018 eq. 4
        ax.plot(dat, K, label=key)

    ax.legend()
    ax.set_yscale('log')
    ax.set_xlabel("Temp ($^oC$)")
    ax.set_ylabel("K(T) ($L^{-1}$)")
    ax.set_title(f"Cumulative Freezing Spectra ({Vwell}$\mu L$ droplet vol)")


# %%
def make_heatmap(data, ntrays=3):
    """
    Generate heatmap of freezing temperatures to troubleshoot
    FINC -- redundant code to MATLAB function in autowell.m
    (produces identical output)
    """
    for key in data.keys():
        fig, ax = plt.subplots()
        the_map = data[key].reshape(6 * ntrays, 16).T
        ax.set_title("FINC Run: " + key)
        ax.set_xticks([])
        ax.set_yticks([])
        fig.colorbar(ax.imshow(the_map, cmap="seismic"))


# %%
def make_small_k(data, Tint=1.0, norm=1, Vwell=10):
    """
    Calculates differential freezing spectra, default delta T = 1.0oC
    See Vali 2018 for calcualtion

    defaults: normalization constant = 1
              well volume = 10 microliter
              temperature interval = 1oC
    """
    fig, ax = plt.subplots()

    for key in data.keys():
        length = len(data[key])
        ind = np.linspace(length, 1, length) / length
        dat = np.sort(data[key], axis=0)

        mintemp = round(min(data[key]))  # create a temp array with Tint spacing
        maxtemp = round(max(data[key]))
        temps = np.arange(mintemp, maxtemp, Tint)

        def big_k(ind):
            big_k = -np.log(1 - ind) / (norm * Vwell * 10 ** -6)  # Vali 2018 eq. 4
            return big_k

        lil_k = np.empty_like(temps)
        for T in temps:
            print(key + " " + str(T))
            # lil_k[i] = (big_k(i)-big_k(i-1)) / Tint # Vali 2018 eq. 14, forward difference

        # ax.plot(dat,lil_k,label=key)

    ax.legend()
    ax.set_yscale("log")
    ax.set_xlabel("Temp ($^oC$)")
    ax.set_ylabel("K(T) ($L^{-1}$)")
    ax.set_title(f"Differential Freezing Spectra ({Vwell}$\mu L$ droplet vol)")


# %%
def main():
    the_filenames = [
        "20201028_lignin2",
        "20201106_lignin3",
        #"20201106_lignin4",
    ]

    # create a dictionary samplename:data, reshape vals to 1D nparray
    the_data = {}
    for file in the_filenames:
        the_data[file[9:]] = np.reshape(get_mat(file), len(get_mat(file)))

    

    # freestyle code to get lignin data ===========================================
    annas_data = "C:/Users/Owner/UBC_F2020/FINC/outputs/20201202_lignin_datadump/anna.mat"
    mat = sio.loadmat(annas_data)
    for key in mat.keys():
        if str(key)[0] == 'F':
            the_data[key] = mat[key]



    jons_data = "C:/Users/Owner/UBC_F2020/FINC/outputs/20201202_lignin_datadump/jon.mat"
    mat = sio.loadmat(jons_data)
    for key in mat.keys():
        if str(key)[0] == 'A':
            the_data[key] = mat[key]
        elif str(key)[0] == 'S':
            the_data[key] = mat[key]
        
    
    sophies_data = "C:/Users/Owner/UBC_F2020/FINC/outputs/20201202_lignin_datadump/sophie.mat"
    mat = sio.loadmat(sophies_data)
    for key in mat.keys():
        for i in range(8):
            if str(key) ==  str('lignin_sophie' + str(i)):
                the_data[key] = mat[key]


    # ==============================================================================
    make_hist(the_data)
    # make_big_K(the_data)
    # make_small_k(the_data)
    # make_boxplot(the_data)
    # make_ff_curve(the_data)
    # make_heatmap(the_data)


if __name__ == "__main__":
    main()


# %%
