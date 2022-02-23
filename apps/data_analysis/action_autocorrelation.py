# %%
import numpy as np
import matplotlib.pyplot as plt
import csv
import autocorrelation
# %%
filename = "/Users/pauldruce/Dev/RanNCGinC++/output/Simulation_1_1_N_6/action_data1_1_N_6/action_data_1_1_N_6_g2_0.000_2022-02-20.csv" 
# %%
action_data = []
with open(filename, 'r') as f:
    csv_reader = csv.reader(f, delimiter=' ')
    for row in csv_reader:
        action_data.append([float(row[0]),float(row[1])])
action_data = np.array(action_data)
# %%
print(action_data)
# %%
action_values = action_data[:][:,0]
acceptance_rates = action_data[:][:,1]
# %%
acf = autocorrelation.autocorr_func_1d(action_values[:5000],True)

# %%
fig = plt.figure(figsize=(5,4), dpi=400)
axes = fig.add_axes([0,1,0.9,0.9])
axes.axhline(0, color='black', linestyle='--', linewidth=1)
axes.set_title("Auto-correlation")
axes.plot(acf)
# %%
