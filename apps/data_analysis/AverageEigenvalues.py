# %% [markdown] # Average Eigenvalues from Simulation output The output of the simulation are currently the raw Dirac
# operators. The spectral properties of the Dirac operators (specific properties of the eigenvalues) are closely
# related to the geometry they describe. So we would like to measure the eigenvalues from the simulation data. This
# is done alone with the auto-correlation function of the eigenvalues.

# %%
import numpy as np
import matplotlib.pyplot as plt
import h5py
from autocorrelation import *
# from pyarma import *
# %%
f=h5py.File(
    "/Users/pauldruce/Dev/RanNCGinC++/output/Simulation_1_1_N_6/dirac_matrices_1_1_N_6/Dirac_Matrices_1_1_N_6_g2_2.400.h5", 'r')
# %%
diracs = []
for i in range(1,2000):
    dirac_matrix = f['dirac_{}'.format(i)]
    diracs.append( dirac_matrix['real'] + 1.0j*dirac_matrix['imag'])
# %%
all_eigenvalues = []
for d in diracs:
    eigenvalues =  np.linalg.eigvals(d)
    eigenvalues.sort()
    all_eigenvalues.append(eigenvalues)
# %%
import pandas as pd

# %%
df = pd.DataFrame(all_eigenvalues)

# %%
np_all_eigenvalues = np.array(all_eigenvalues)

# %%
av_eigenvalues = []
for i in range(len(np_all_eigenvalues[0,:])):
    av_eigenvalues.append(np.average(np_all_eigenvalues[:,i]))
# %%
av_eigenvalues
# %% [markdown]
#  ## Autocorrelation Function for the eigenvalues.
#  In order to calculate the statical error associated with the measurements of the eigenvalues, the autocorrelation function needs to be calculated and the integrated autocorrelation time needs to be calculated. The
#  The autocorrelation function below makes use of the Fourier transform to speed up calculation. It is taken from here: https://emcee.readthedocs.io/en/stable/tutorials/autocorr/
#  
#  TODO: Note sure this is the correct way to calculate the autocorrelation function, as there seems to be a lot of noise compared to the various other plots I've seen.

# %%
acfs = []
for i in range(len(np_all_eigenvalues[0,:])):
    eigenvalues = np_all_eigenvalues[:,i]
    acf  = autocorr_func_1d(eigenvalues)
    acfs.append(acf)
# %%
x_values = []
for i in range(1,2000):
    x_values.append(i*20)
# %%
fig = plt.figure(figsize=(5,4), dpi=400)
axes = fig.add_axes([0,len(acf),0.9,0.9])
axes.set_ylim(-1,1)
axes.set_xlim(0,2000*20)
axes.axhline(0, color='black', linestyle='--', linewidth=1)
axes.set_title("Autocorrelation") 
for i in range(len(np_all_eigenvalues[0,:])):
    axes.plot(x_values, acfs[i])
# %%

# %%
