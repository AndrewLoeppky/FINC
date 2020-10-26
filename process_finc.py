# %%
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import scipy.io

# %%
def get_mat(the_path):
    mat = scipy.io.loadmat(the_path)
    frz_t = pd.DataFrame(np.hstack((mat['FrzTall'])))
    plt.plot(frz_t)
    
# %%
get_mat('C:/Users/Owner/UBC_F2020/FINC/plots/finc_20201007_SA1/20201007_SA_andrew1.mat')


# %%
