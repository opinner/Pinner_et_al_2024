import matplotlib.pyplot as plt
import numpy as np
import matplotlib.colors as mcolors
import pandas as pd
import scipy.stats as ss
from matplotlib.markers import MarkerStyle

# read Thorpe and strain results data
#-------------------------------------------------------------------
data = np.load("/home/ole/Desktop/Mooring_Analysis/energy_levels/data/Thorpe_result.npz", allow_pickle=True)
mab = data["mab"]

data = np.load("/home/ole/Desktop/CTD/mixsea/Thorpe_vertical_eps.npz")
thorpe_z = data["z"]
thorpe_eps = data["eps"]

data = np.load("/home/ole/Desktop/CTD/mixsea/Strain_vertical_eps.npz")
strain_z = data["z"]
strain_eps = data["eps"]


#TODO!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!TODO#
strain_eps = strain_eps*2.694 #Correction from Rw =3 to Rw = 7
#TODO!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!TODO#


#Horizontally average of the mooring results
#-------------------------------------------------------------------
energy_levels = pd.read_csv("/home/ole/Desktop/Mooring_Analysis/energy_levels/wave_energy_result.csv")
mab_ranges = [(20,60),(120,160),(320,400)]
mabs = []
IGWs = []
for mab_range in mab_ranges:
    inbetween_mabs = []
    inbetween_IGWs = []
    for mab, IGW in zip(energy_levels["mab"], energy_levels["eps_IW"]):
        if not (mab > mab_range[0] and mab < mab_range[1]):
            continue
        inbetween_mabs.append(mab)        
        inbetween_IGWs.append(IGW) 
    mabs.append(np.mean(inbetween_mabs))
    IGWs.append(np.mean(inbetween_IGWs))
    
    
cm = 1/2.54  # centimeters in inches
ONE_COLUMN_WIDTH = 8.3
GOLDEN_RATIO = 1.61
fig, ax = plt.subplots(1, figsize=(ONE_COLUMN_WIDTH*cm, GOLDEN_RATIO*ONE_COLUMN_WIDTH*cm))

ax.semilogx(thorpe_eps, thorpe_z, c = "tab:red", label = "$\\langle\\varepsilon_{\\mathrm{total, Thorpe}}\\rangle$")
ax.semilogx(IGWs, mabs, "o", c = "k", alpha = 0.6, label = "$\\langle\\varepsilon_{\\mathrm{IGW, IDEMIX}}\\rangle$", zorder = 10)    
ax.semilogx(strain_eps, strain_z, linestyle='-', marker='.', c = "tab:blue", label = "$\\langle\\varepsilon_{\\mathrm{IGW, strain}}\\rangle$")


ax.set_ylabel("Distance from sea floor (m)")

ax.set_xlabel("$\\langle\\varepsilon\\rangle\\,$(W kg$^{-1})$")
ax.set_ylim(-10,450)
ax.legend(loc = "upper left", framealpha = 0.6, labelspacing = 0.4, ncols = 1, fontsize="9")
#lims = ax.get_xlim()
#ax.set_xlim((0.95*lims[0],lims[1]))
fig.tight_layout()
fig.savefig("./avrg_vertical_eps.png", dpi = 400, bbox_inches = "tight")
fig.savefig("./avrg_vertical_eps.svg", bbox_inches = "tight")
plt.show()




   





