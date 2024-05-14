import matplotlib.pyplot as plt
import numpy as np
import matplotlib.colors as mcolors
import pandas as pd
import cmocean



def main():

    ONE_COLUMN_WIDTH = 8.3
    TWO_COLUMN_WIDTH = 12
    GOLDEN_RATIO = 1.61
    cm = 1/2.54  # centimeters in inches



    data = np.load("../scripts/IDEMIX_parametrization/method_data/results_available_energy.npz", allow_pickle = True)     
    energy_levels = pd.DataFrame(data = {
        "lat":data["lat"], 
        "lon":data["lon"], 
        "tidal":data["tidal_energies"],
        "rounded_mab": data["mab"],
        "cats":data["cats"],
        "barotropic":data["barotropic"],
    })

    fig,ax = plt.subplots(1, figsize=(TWO_COLUMN_WIDTH*cm, 0.8*TWO_COLUMN_WIDTH*cm))


    ax.plot(energy_levels["lon"],  energy_levels["cats"], c = "tab:blue", lw = 3, label = "CATS model")
    ax.plot(energy_levels["lon"], energy_levels["barotropic"], lw = 3, color = "k", label = "assumed barotropic energy", )
    
    # plot x-shifted tidal energy level
    groups = energy_levels.groupby("lon")
    for i,lon in enumerate(set(energy_levels["lon"])):
    
        group = groups.get_group(lon)
        group = group.sort_values("rounded_mab")
        #print(group["rounded_mab"])
        nr_of_points = len(group.index)
    
        #to declutter the figure, points are shifted centered around the original x value. 
        # magnitude of the shift in units of the x coordinate
        xstep = 0.1
        xshift = np.arange(0,nr_of_points*xstep, xstep) 
        #distinction between 1 point (no shifts), odd and even points
        if nr_of_points%2 != 0 and len(xshift)!= 1:
            xshift = xshift[:-1]
        #centering    
        xshift = xshift - np.median(xshift)
    
        #for a, b in zip(xshift, group["rounded_mab"]):
        #    print(a,b)
    
        label = None
        if i == 0:
            label = "measured"
        ax.scatter(group["lon"] + xshift, group["tidal"], marker = ".", c = "tab:gray", ls = "None", s = 200, label = label, edgecolor = "k", zorder = 5)    


    #ax.plot(longitudes,  available_energies, "D")

    #ax.fill_between(energy_levels["lon"],  [0]*len(energy_levels["barotropic"]), energy_levels["barotropic"], label = "assumed barotropic energy", hatch= "//", facecolor = "none")
    ax.legend()
    ax.ticklabel_format(axis = "y", scilimits = (0,0), style = "scientific", useMathText=True)  
    ax.yaxis.get_offset_text().set_x(-0.1)

    ax.set_title("Energy at semidurnal tidal frequencies")
    ax.set_xlabel("Longitude (°)")
    ax.set_ylabel(r"$E$ / (m$^2$ s$^{-2}$)")  
    fig.tight_layout()
    fig.savefig("./barotropic_tide.png", dpi = 400, bbox_inches = "tight")
    
    
if __name__ == '__main__':
    main()
    plt.show()
