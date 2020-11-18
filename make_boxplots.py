"""
===============================================================
DATA PROCESSING LIBRARY - FINC OUTPUTS

Author: Andrew Loeppky
UBC/ETH Atmospheric Science, NBD Group
Last Update Nov 18/2020

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
    '''
    output a freezing temperature dictionary filename:nparray()
    from a csv. CSVs are assumed to be in the format:
    
    metadata, metadata, ...
    name1,   name2, ...,   nameN
    frztemp, frztemp, ..., frztemp
    .
    .
    .
    As per ETH repository entries found at: 
    https://www.research-collection.ethz.ch/handle/20.500.11850/438875
    '''
    reader = csv.DictReader(open(the_file))

    result = {}
    for row in reader:
        for column, value in row.items(): 
            result.setdefault(column, []).append(value)
    return result

# %%
def make_hist(*frzdata):
    """
    Returns a histogram of freezing events at temp T

    WARNING: NOT NORMALIZED TO N(T)!!! Strongly dependent on
    droplet size
    """
    for data in frzdata:
        plt.hist(data, bins=100)

    plt.xlabel("Freezing Temperature ($^oC$)")
    plt.ylabel("$N_{wells}$")
    plt.title("Frozen Well Count as a Function of Temperature")


# %%
def make_boxplot(data):
    """
    given one or more input freezing vectors, make boxplots of
    freezing temperature
    """
    """
    for key in data.keys():
        data[key].reshape(len(data[key]))
    """
    fig, ax = plt.subplots()
    ax.boxplot(data.values(), showmeans=True)
    ax.set_xticklabels(data.keys(), rotation=-60)
    ax.set_ylabel("Freezing Temp ($^oC$)")
    ax.set_title("Freezing Temp Distrubutions")


# %%
def make_ff_curve(data):
    """
    creates a frozen fraction (FF) curve as a function of
    temperature
    """
    for dat in data.values():
        dat = np.sort(dat, axis=0)
    fig, ax = plt.subplots()

    length = len(data)
    ind = np.linspace(length, 1, length) / length

    for dat in data:
        dat = np.sort(dat, axis=0)
        plt.plot(dat, ind)

    plt.xlabel("Temp ($^oC$)")
    plt.ylabel("FF")
    plt.title("Frozen Fraction")


# %%
def make_big_K(*data, norm=1, Vwell=10):
    """
    creates a cumulative concentration curve from Vali eq.

    default normalization constant (eg surface area) = 1
    default well volume = 10 microlitre
    """
    # make frozen fraction curve and calculate K from it
    length = len(data[1])
    ind = np.linspace(length, 1, length) / length

    K = -np.log(1 - ind) / (norm * Vwell * 10 ** -6)  # Vali eq.

    for dat in data:
        dat = np.sort(dat, axis=0)
        plt.plot(dat, K)

    plt.yscale("log")
    plt.xlabel("Temp ($^oC$)")
    plt.ylabel("K(T) ($L^{-1}$)")
    plt.title(f"Cumulative Freezing Spectra ({Vwell}$\mu L$ droplet vol)")


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
        ax.set_title('FINC Run: ' + key)
        ax.set_xticks([])
        ax.set_yticks([])
        fig.colorbar(ax.imshow(the_map, cmap="seismic"))


# %%
def make_small_k(*data, Tint=1.0):
    """
    Calculates differential freezing spectra, default delta T = 1.0oC
    See Vali 2018 for calcualtion
    """

    return None


# %%
def main():
    the_filenames = [
        "20201007_SA1",
        #"20201007_SA2", this one is bunk!
        "20201007_SA3",
        "20201020_SA2",
        "20201020_SA3",
        "20201111_SArepeat-1",
        "20201117_SA-1",
        "20201117_SA-2",
        "20201117_SA-3",
    ]
    
    # create a dictionary samplename:data, reshape vals to 1D nparray
    the_data = {}
    for file in the_filenames:
        the_data[file[4:]] = np.reshape(get_mat(file), len(get_mat(file)))

    # make_hist(*the_data)
    # make_big_K(*the_data)
    make_boxplot(the_data)
    # make_ff_curve(the_data)
    # make_heatmap(the_data)


if __name__ == "__main__":
    main()


# %%
