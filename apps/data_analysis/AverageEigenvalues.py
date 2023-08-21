# %% [markdown] 
# # Average Eigenvalues from Simulation output The output of the simulation are currently the raw Dirac
# operators. The spectral properties of the Dirac operators (specific properties of the eigenvalues) are closely
# related to the geometry they describe. So we would like to measure the eigenvalues from the simulation data. This
# is done alone with the auto-correlation function of the eigenvalues.

# %%
import matplotlib.pyplot as plt
import h5py
from autocorrelation import *
import os


# %%
g2start = 0
g2end = -4.000
g2step = -0.04
g2values = []
g2 = g2start
while g2>g2end:
    g2values.append(g2)
    g2+= g2step


# %%
def create_png_of_av_eigenvalues(g2):
    if not os.path.exists(output_dir+simulation_subdir+"images/"):
        os.makedirs(output_dir+simulation_subdir+"images/")

    eigenvalues_file = "Eigenvalues_{}_{}_N_{}_g2_{:.3f}.h5".format(p,q,N,-1*g2) 
    filename = output_dir+simulation_subdir+eigenvalues_dir + eigenvalues_file
    f=h5py.File(filename, 'r')

    all_eigenvalues = []
    for i in range(1,2000):
        eigenvalues = f['eigenvalues_{}'.format(i)]
        # print(np.shape(np.array(eigenvalues[0])))
        all_eigenvalues.append( eigenvalues[0])

    np_all_eigenvalues = np.array(all_eigenvalues)

    av_eigenvalues = []
    for i in range(len(np_all_eigenvalues[0,:])):
        av_eigenvalues.append(np.average(np_all_eigenvalues[:,i]))
    
    #  ## Autocorrelation Function for the eigenvalues.
    #  In order to calculate the statical error associated with the measurements of the eigenvalues, the autocorrelation function needs to be calculated and the integrated autocorrelation time needs to be calculated. The
    #  The autocorrelation function below makes use of the Fourier transform to speed up calculation. It is taken from here: https://emcee.readthedocs.io/en/stable/tutorials/autocorr/
    #  
    #  TODO: Note sure this is the correct way to calculate the autocorrelation function, as there seems to be a lot of noise compared to the various other plots I've seen.


    acfs = []
    integrated_auto_corrs  = []
    variances = []
    errors = []
    for i in range(len(np_all_eigenvalues[0,:])):
        eigenvalues = np_all_eigenvalues[:,i]

        # Variances
        variance = np.var(eigenvalues)
        variances.append(variance)

        # Autocorrelations
        acf  = autocorr_func_1d(eigenvalues)
        acfs.append(acf)

        # Integrated autocorrelation
        int_auto_cor = 1 + 2*np.sum(acf)
        integrated_auto_corrs.append(int_auto_cor)

        num_samples = len(np_all_eigenvalues[:,i])

        error = np.sqrt(2*variance * int_auto_cor/num_samples)
        errors.append(error)

    x_values = []
    for i in range(1,2000):
        x_values.append(i*20)
    
    # fig = plt.figure(figsize=(5,4), dpi=400)
    # axes = fig.add_axes([0,len(acf),0.9,0.9])
    # axes.set_ylim(-1,1)
    # axes.set_xlim(0,2000*20)
    # axes.axhline(0, color='black', linestyle='--', linewidth=1)
    # axes.set_title("Autocorrelation") 
    # for i in range(len(np_all_eigenvalues[0,:])):
    #     axes.plot(x_values, acfs[i])
    
    # some directives for prettier plots
    # standardsize=18
    # ticksize=16
    # legendsize=12
    # axeslabelsize=18
    # plt.rc('text', usetex=True)  ##this is the magic
    # font = {'family' : 'serif',
    #         'weight' : 'normal',
    #         'size'   : standardsize}
    # plt.rc('font', **font)
    # plt.rcParams['xtick.labelsize'] = ticksize
    # plt.rcParams['ytick.labelsize'] = ticksize
    # plt.rcParams.update({'figure.autolayout': True})

    x_val = []
    for i in range(len(av_eigenvalues)):
        x_val.append(i)

    fig = plt.figure(figsize=(5,4),dpi=400)
    axes = fig.add_subplot()
    axes.set_ylim([-1.5, 1.5])

    axes.errorbar(x_val, av_eigenvalues, yerr=errors, fmt=".",ecolor="black", capsize=5)
    fig.savefig(output_dir+simulation_subdir+"images/av_eigenvalues_{}_{}_N_{}_g2_{:.3f}.png".format(p,q,N, -1*g2))
    plt.close(fig)


# %%
# %%
p=1
q=1
output_dir = "/Users/pauldruce/Dev/RanNCGinC++/output/"
for N in range(5,11,1):
    simulation_subdir = "Simulation_{}_{}_N_{}/".format(p,q,N)
    eigenvalues_dir = "eigenvalues_{}_{}_N_{}/".format(p,q,N)
    for g2 in g2values:
        create_png_of_av_eigenvalues(g2)
# %%
