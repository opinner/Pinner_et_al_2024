import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

# read Thorpe and strain results data
#-------------------------------------------------------------------
#data = np.load("/home/ole/Desktop/Mooring_Analysis/energy_levels/data/Thorpe_result.npz", allow_pickle=True)
#mab = data["mab"]

data = np.load("../scripts/thorpe_scales/method_data/horizontally_averaged_Thorpe_eps.npz")
thorpe_z = data["z"]
thorpe_eps = data["eps"]

data = np.load("../scripts/shear_strain_parametrization/method_data/Strain_vertical_eps.npz")
strain_z = data["z"]
strain_eps = data["eps"]


#TODO!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!TODO#
strain_eps = strain_eps*2.694 #Correction from Rw =3 to Rw = 7
#TODO!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!TODO#


#Horizontally average of the mooring results
#-------------------------------------------------------------------
energy_levels = pd.read_csv("../scripts/IDEMIX_parametrization/method_data/eps_IGW_IDEMIX_results.csv")
mab_ranges = [(20, 60),(120, 160),(320, 400)]
mabs = []
IGWs = []
IGW_errors = []
for mab_range in mab_ranges:
    inbetween_mabs = []
    inbetween_IGWs = []
    inbetween_IGW_errors = []
    for mab, IGW, IGW_error in zip(energy_levels["rounded_mab"], energy_levels["eps_IGW"], energy_levels["eps_IGW_error"]):
        if not (mab > mab_range[0] and mab < mab_range[1]):
            continue
        inbetween_mabs.append(mab)        
        inbetween_IGWs.append(IGW)
        inbetween_IGW_errors.append(IGW_error)
    mabs.append(np.mean(inbetween_mabs))
    IGWs.append(np.mean(inbetween_IGWs))
    IGW_errors.append(np.mean(inbetween_IGW_errors))

cm = 1/2.54  # centimeters in inches
ONE_COLUMN_WIDTH = 8.3
GOLDEN_RATIO = 1.61
fig, ax = plt.subplots(1, figsize=(ONE_COLUMN_WIDTH*cm, GOLDEN_RATIO*ONE_COLUMN_WIDTH*cm))

ax.semilogx(thorpe_eps, thorpe_z, c = "tab:red", label = "$\\langle\\varepsilon_{\\mathrm{total, Thorpe}}\\rangle$")
ax.semilogx(IGWs, mabs, "o", c = "k", label = "$\\langle\\varepsilon_{\\mathrm{IGW, IDEMIX}}\\rangle$", zorder = 10)

label = "$\\langle\\varepsilon_{\\mathrm{IGW, strain}}\\rangle$"
for x, z in zip(strain_eps, strain_z):
    spacing = np.abs(np.mean(np.diff(strain_z)))
    ax.vlines(x, z-spacing/2, z+spacing/2, lw = 3,  colors = "tab:blue", label = label)
    label = None # use only the label from the first plot instance
    ax.fill_betweenx([z-spacing/2, z+spacing/2], x/5, x*5, color = "tab:blue", alpha = 0.3)
#ax.semilogx(strain_eps, strain_z, 'o', ms = 6, mec = "xkcd:darkblue")

#

# Plot Errors
for IGW, IGW_error, mab in zip(IGWs, IGW_errors, mabs):
    ax.plot([IGW-IGW_error, IGW+IGW_error], [mab,mab], lw = 3, c = "k", alpha = 0.5)
ax.fill_betweenx(thorpe_z, thorpe_eps/5, thorpe_eps*5, color = "tab:red", alpha = 0.5)



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




   





