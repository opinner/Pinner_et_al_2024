import matplotlib.pyplot as plt
#import matplotlib as mpl
#mpl.rcParams.update(mpl.rcParamsDefault)
#plt.rcParams['text.usetex'] = True
import numpy as np
import cmocean as cm
import helper_functions as help

"""
def path_as_footnote():
    # Add a footnote below and to the right side of the chart
    import os
    file_path = os.path.realpath(__file__)
    #plt.figtext(0.99, 0.5, file_path, horizontalalignment='right', va = "center",\
    #    figure = plt.gcf(), rotation = "vertical", fontsize = 6)
    plt.figtext(0.99, 0, file_path, horizontalalignment='right', va = "bottom",\
        figure = plt.gcf(), fontsize = 6)
"""   
                        
def main(): 

    data = np.load("./integrated_levels.npz", allow_pickle = True) 
    barotropic = data["barotropic"]
    baroclinic = data["baroclinic"]
    IW = data["IW"]
    lons = data["lons"] 
    rounded_depths = data["rounded_depths"]


    #Load temperature information
    data = np.load("./data/Temperature_at_velocity_locations.npz")#, allow_pickle = True) 
    uv_lon = data["lon"] 
    uv_mab = data["mab"] 
    uv_depth = data["depth"] 
    uv_temperature = data["temperature"]

    #synchronize the temperature information with the information about the energy levels
    temperatures = []
    for lon, depth in zip(lons, rounded_depths):
        flag = False
        for uvl, uvd, uvt in zip(uv_lon, uv_depth, uv_temperature):
            #print(depth,f"{uvd:.0f}")
            if depth == f"{uvd:.0f}":
                temperatures.append(uvt)
                flag = True
                break    
    temperatures = np.array(temperatures) 
                    
    #print(np.shape(lons),np.shape(barotropic),np.shape(baroclinic),np.shape(IW))
    #print(lons,"\n",barotropic,"\n",baroclinic,"\n",IW)
    
    f,ax = plt.subplots(2,sharex = True)
    temperature_boundary = -0.7
    below_boundary = temperatures<=-0.7
    ax[0].plot(lons, barotropic,"D", c = "tab:blue", label = "Barotropic", alpha = 0.8)      
    ax[0].plot(lons[~below_boundary], baroclinic[~below_boundary],"o", markeredgecolor = "k", markeredgewidth=1.5, markerfacecolor = "None", label = r"Baroclinic ($T>0.7\,$°C)")  
    ax[0].plot(lons[below_boundary], baroclinic[below_boundary],"x", c = "k", label = r"Baroclinic ($T<0.7\,$°C)")  
        
    ax[1].plot(lons[~below_boundary], IW[~below_boundary],"o", markeredgecolor = "k", markeredgewidth=1.5, markerfacecolor = "None", label = r"IW, $T>0.7\,$°C")   
    ax[1].plot(lons[below_boundary], IW[below_boundary],"x", c = "k", label = r"IW, $T<0.7\,$°C")
    
    ax[0].set_title("Tides")
    ax[1].set_title(r"Internal wave continuum ($f$ - $N$)")
    
    ax[0].legend()
    ax[1].legend()

    ax[0].set_ylabel("Energy / (m²/s²)")
    ax[1].set_ylabel("Energy / (m²/s²)")
            
    ax[1].set_xlabel("Longitude")
    f.tight_layout()
    
    import os
    path = os.path.realpath(__file__)           
    help.Plot.path_as_footnote(fig = f, path = path)
    #f.savefig(f"./integrated_energy_levels",dpi = 300)
    
    """
    f,ax = plt.subplots(1)#,sharex = True)
    ax.semilogy(lons, barotropic,"D", label = "barotropic")      
    ax.semilogy(lons, baroclinic,"o", label = "baroclinic")  
    ax.semilogy(lons, IW,"X", c = "tab:red", label = "IW")
    
    ax.legend()

    ax.set_ylabel("Energy / (m²/s²)")
            
    ax.set_xlabel("Longitude")

    help.Plot.path_as_footnote(fig = f)
    f.savefig(f"./integrated_energy_levels_log",dpi = 300)
    """
    
    
if __name__ == '__main__':
    main()
    plt.show()
