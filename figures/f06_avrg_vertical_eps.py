import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

plt.style.use('./paper.mplstyle')


# read Thorpe and strain results data
#-------------------------------------------------------------------

data = np.load("../scripts/thorpe_scales/method_results/horizontally_averaged_Thorpe_eps.npz")
thorpe_z = data["z"]
thorpe_eps = data["eps"]

data = np.load("../scripts/finestructure/method_results/Strain_vertical_eps2.npz")
strain_mab = data["mab"]
strain_eps = data["eps"]

#Horizontally average of the mooring results
#-------------------------------------------------------------------
energy_levels = pd.read_csv("../scripts/IDEMIX_parametrization/method_results/eps_IGW_IDEMIX_results.csv")

error = energy_levels["eps_IGW_mult_error"]
print("E Error", np.mean(error), np.std(error), np.max(error), np.min(error))

mab_ranges = [(20, 60), (120, 160), (320, 400)]
mabs = []
IGWs = []
IGW_errors = []
for mab_range in mab_ranges:
    inbetween_mabs = []
    inbetween_IGWs = []
    inbetween_IGW_errors = []
    for mab, IGW, IGW_error in zip(energy_levels["rounded mab"], energy_levels["eps_IGW"],
                                   energy_levels["eps_IGW_mult_error"]):
        if not (mab > mab_range[0] and mab < mab_range[1]):
            continue
        inbetween_mabs.append(mab)
        inbetween_IGWs.append(IGW)
        inbetween_IGW_errors.append(IGW_error)
    mabs.append(np.mean(inbetween_mabs))
    IGWs.append(np.mean(inbetween_IGWs))
    IGW_errors.append(np.mean(inbetween_IGW_errors))

cm = 1 / 2.54  # centimeters in inches
ONE_COLUMN_WIDTH = 8.3
GOLDEN_RATIO = 1.61
fig, ax = plt.subplots(1, figsize=(ONE_COLUMN_WIDTH * cm, GOLDEN_RATIO * ONE_COLUMN_WIDTH * cm))

ax.semilogx(thorpe_eps, thorpe_z, c="tab:red", label="$\\langle\\varepsilon_{\\mathrm{total, Thorpe}}\\rangle$")
ax.semilogx(IGWs, mabs, "o", c="k", label="$\\langle\\varepsilon_{\\mathrm{IGW, IDEMIX}}\\rangle$", zorder=10)

label = "$\\langle\\varepsilon_{\\mathrm{IGW, fine}}\\rangle$"
for x, z in zip(strain_eps, strain_mab):
    spacing = np.abs(np.mean(np.diff(strain_mab)))
    assert spacing == 125
    upper = z + spacing / 2
    lower = z - spacing / 2
    #if lower < 0: lower = 0
    ax.vlines(x, lower, upper, lw=3, colors="tab:blue", label=label)
    label = None  # use only the label from the first plot instance
    ax.fill_betweenx([lower, upper], x / 5, x * 5, color="tab:blue", alpha=0.3)
#ax.semilogx(strain_eps, strain_z, 'o', ms = 6, mec = "xkcd:darkblue")

# Plot multiplicative errors
for IGW, IGW_error, mab in zip(IGWs, IGW_errors, mabs):
    ax.plot([IGW / IGW_error, IGW * IGW_error], [mab, mab], lw=3, c="k", alpha=0.5)
ax.fill_betweenx(thorpe_z, thorpe_eps / 5, thorpe_eps * 5, color="tab:red", alpha=0.5)

ax.set_ylabel("Height above bottom (m)")

ax.set_xlabel("$\\langle\\varepsilon\\rangle\\,$(W kg$^{-1})$")
ax.set_ylim(-10, 450)
ax.legend(loc="upper left", framealpha=0.6, labelspacing=0.4, ncols=1, fontsize="9")
ax.set_xlim(3e-11,5e-8)
fig.tight_layout()
#fig.savefig("./avrg_vertical_eps.png", dpi=400, bbox_inches="tight")
fig.savefig("./avrg_vertical_eps.pdf")
plt.show()
