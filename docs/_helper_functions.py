#TODO add comments
import scipy.io as sio
import scipy.signal as signal
import scipy as sc
import matplotlib.pyplot 
import numpy as np
import datetime as dt
import matplotlib.dates as mdates #import DateFormatter
import matplotlib.pyplot as plt

#-------------------------------------------
#Functions for loading data 
#-------------------------------------------

class IO:
    import pickle
    
    @staticmethod 
    def save_as_pickle(figure, name= "figure"):
        # save the mpl figure as pickle format
        with open(f"{name}.pkl", 'wb') as fs:
            pickle.dump(figure, fs)
        return

    @staticmethod
    def load_pickle(name = None):
        with open(f"{name}.pkl", 'rb') as fs:
            return pickle.load(fs)

    @staticmethod                 
    def get_filepaths_from_directory(directory = None, inclusive = "", exclusive = ()):

        #endswith method requires string or tuple of strings
        if type(inclusive) is not str and type(inclusive) is not tuple:
            inclusive = tuple(inclusive)
        if not type(inclusive) is not str and type(inclusive) is not tuple:
            exclusive = tuple(exclusive)

        def all_files(folder):
            import os
            for path, dirs, files in os.walk(folder):
                for f in files:
                    yield os.path.join(path, f)

        list_of_files = [f for f in all_files(directory)]

        return [f for f in list_of_files if f.endswith(inclusive) and not f.endswith(exclusive)]

    
    
    def meta_data_lookuptable():

        ID_keys = [1133701,1135001,1270602] #Serial number of the sensor

        meta_data_keys = ["ID", "lat", "lon", "max_depth", "mab", "depth"]  
        meta_data = [[1133701, -76,  -31,    470,  18, 470-18],
                    [ 1135001, -76,  -31,    470, 118, 470-118],
                    [ 1270602, -76,  -31.5,  600, 200, 600-200]]

        meta_data = [dict(zip(meta_data_keys, e)) for e in meta_data]
        meta_data_dict = dict(zip(ID_keys, meta_data)) # zip together the lists and make a dict from it
        
        return meta_data_dict
        
        
        
    def get_meta_data_for_file(file_path):
        """
        Assumes the ID as filename
        """
        ID = int(file_path.split("/")[-1][:-4])
        meta_data_dict = meta_data_lookuptable()
        meta_data = meta_data_dict[ID]
        lat = meta_data["lat"]
        lon = meta_data["lon"]
        mab = meta_data["mab"]
        depth = meta_data["depth"]
        
        return lat,lon, mab, depth


    def read_Aquadopp(file_path):
        #define date columns
        #subtract 1 to account for starting counts at 0
        Month_column  =  1 - 1
        Day_column    =  2 - 1
        Year_column   =  3 - 1
        Hour_column   =  4 - 1
        #Minute_column =  5 - 1
        #Second_column =  6 - 1
        u_column      =  9 - 1
        v_column      = 10 - 1

        data = np.genfromtxt(file_path, dtype=("|S2","|S2","|S4","|S2", float,float), filling_values = np.nan, usecols = (Month_column,Day_column,Year_column,Hour_column,u_column,v_column))

        #convert read in data to numpy array
        data = np.asarray([[x for x in row] for row in data])

        #define east/west and north/south velocities and remove bad values at the end
        data = np.asarray(data) #convert to numpy array to slice
        u = data[:,-2].astype(float) 
        v = data[:,-1].astype(float)
        
        #join date columns for the parsing
        #test = " ".join(np.asarray(data[0][0:4]).astype(str))
        #print(test)

        #converts the first 4 values from bytes into strings, joins them and parses them to a datetime object
        time = [dt.datetime.strptime(" ".join(np.asarray(row[0:4]).astype(str)), '%m %d %Y %H') for row in data]

        lat, lon, mab, depth = get_meta_data_for_file(file_path)
        
        return {"time":time,"u":u,"v": v,"lat": lat, "lon": lon, "mab": mab, "depth":depth}
        
        
    def read_mat_file(file_path):  
        
        data_struct = sio.loadmat(file_path) 
        data = data_struct["S"][0][0]

        matlab_time = np.squeeze(data["DATETIME"])
        # convert Matlab variable "t" into list of python datetime objects
        time = [matlab2datetime(tval) for tval in matlab_time]

        u = np.squeeze(data["UC"])
        v = np.squeeze(data["VC"])
        
        #fill in NaNs
        u = interpolate_nans(u)
        v = interpolate_nans(v)

        #if velocities are too high, convert to m/s instead of cm/s
        if max(abs(u))>20 or max(abs(v))>20:
            u = u/100
            v = v/100
            
        assert np.shape(u) == np.shape(time)

        lat = data["LAT"][0][0]
        lon = data["LON"][0][0]
        mab = np.nan
        depth = np.nan
        
        return {"time":time,"u":u,"v": v,"lat": lat, "lon": lon, "mab": mab, "depth":depth}
    


    
def get_uv_from_tidal_model(lat,lon, time_span, directory = "/media/sf_VM_Folder/data/tidal_model_data"):

    import locale
    locale.setlocale(locale.LC_ALL,('en_US', 'UTF-8'))
    # HELP! Why does dt.datetime.strptime need this change and why exactly here? Outside the function doesn"t work
    # If I dont set locale here strptime and matplotlib use german month names all of a sudden.
    # Even if this file is imported as an external library
    
    
    paths = get_filepaths_from_directory(directory = directory, inclusive = ".model.mat")
        
    #check if filename contains the coordinates with one decimal place accurateness (rounded, not truncated)
    for p in paths:
        if f"u_{lat:.1f}_{lon:.1f}" in p:
            u_path = p
        if f"v_{lat:.1f}_{lon:.1f}" in p:
            v_path = p
          
    u_data = sio.loadmat(u_path) 
    u = u_data["TimeSeries"][0] #in units of cm/s
    u_time_as_strings = u_data["Time"]
    u_read_in_lat = u_data["lat"][0][0]
    u_read_in_lon = u_data["lon"][0][0]
    
    v_data = sio.loadmat(v_path) 
    v = v_data["TimeSeries"][0] #in units of cm/s
    v_time_as_strings = v_data["Time"]
    v_read_in_lat = v_data["lat"][0][0]
    v_read_in_lon = v_data["lon"][0][0]

    assert np.all(u_time_as_strings == v_time_as_strings) #Check if u and v share the same timestamps  
    assert np.allclose([u_read_in_lat,u_read_in_lon],[v_read_in_lat, v_read_in_lon], atol = 0.1) #Check if lat, lon from u and v files are equal
    assert np.allclose([lat,lon],[v_read_in_lat, v_read_in_lon], atol = 0.1)  #Check if desired and returned lat, lon are equal
    

    #converts the first 4 values from bytes into strings, joins them and parses them to a datetime object 
    time = np.asarray([dt.datetime.strptime(string, '%d-%b-%Y %H:%M:%S') for string in u_time_as_strings])
    
    #convert from cm/s to m/s  
    u = u/100
    v = v/100
    
    #trim the model time series to the length of the measured time series
    start_index = np.argmin(np.abs(time-time_span[0]))
    stop_index = np.argmin(np.abs(time-time_span[-1]))

    return time[start_index:stop_index], u[start_index:stop_index], v[start_index:stop_index]

#-------------------------------------------
#Functions for Conversions
#-------------------------------------------
class Conversions:  
    @staticmethod  
    def matlab2datetime(matlab_datenum):
        day = dt.datetime.fromordinal(int(matlab_datenum))
        dayfrac = dt.timedelta(days=matlab_datenum%1) - dt.timedelta(days = 366)
        return day + dayfrac    
    
    @staticmethod    
    def timedelta_to_list_of_ints(td):
        #returns days, hours, minutes, seconds as ints
        return td.days, td.seconds//3600, (td.seconds//60)%60, td.seconds%60
    
    @staticmethod    
    def cart2pol(x, y):
        rho = np.sqrt(x**2 + y**2)
        phi = np.arctan2(y, x)
        return(rho, phi)

    @staticmethod
    def pol2cart(rho, phi):
        x = rho * np.cos(phi)
        y = rho * np.sin(phi)
        return(x, y)
    
    
#-------------------------------------------
#working with data arrays
#-------------------------------------------    
class Data:
    @staticmethod
    def cut_trailing_nans(array):
        #if array has no nans, just return the input
        if not np.any(np.isnan(array)):
            return array  

        first_nan_index = np.argmax(np.isnan(array)) 
        assert first_nan_index != 0
        #Sanity check  
        assert np.sum(np.isnan(array)) == len(array)-first_nan_index,"There seems to be more Nans outside the padding"
        return array[:first_nan_index]
        
    def nan_helper(y):
        """Helper to handle indices and logical indices of NaNs.

        Input:
            - y, 1d numpy array with possible NaNs
        Output:
            - nans, logical indices of NaNs
            - index, a function, with signature indices= index(logical_indices),
              to convert logical indices of NaNs to 'equivalent' indices
        Example:
            >>> # linear interpolation of NaNs
            >>> nans, x= nan_helper(y)
            >>> y[nans]= np.interp(x(nans), x(~nans), y[~nans])
        """

        return np.isnan(y), lambda z: z.nonzero()[0]
        
    def interpolate_nans(array):
        nans, helper_index = nan_helper(array)
        array[nans]= np.interp(helper_index(nans), helper_index(~nans), array[~nans])
        return array    


    
    #code from https://scipy-cookbook.readthedocs.io/items/ButterworthBandpass.html



    def butter_bandpass(lowcut, highcut, fs, order=5):
        from scipy.signal import butter
        nyq = 0.5 * fs
        low = lowcut / nyq
        high = highcut / nyq
        b, a = butter(order, [low, high], btype='band')
        return b, a


    def butter_bandpass_filter(data, lowcut, highcut, fs, order=5):
        """
        # Sample rate and desired cutoff frequencies (in Hz).
        fs = 1/7200 #sampling period of 2 hours
        lowcut = 1/(14*60*60) #period of 18 hours
        highcut = 1/(10*60*60) #period of 10 hours
        filtered_u = butter_bandpass_filter(patched_u, lowcut, highcut, fs, order=6)
        filtered_v = butter_bandpass_filter(patched_v, lowcut, highcut, fs, order=6)
        """
        from scipy.signal import lfilter
        b, a = butter_bandpass(lowcut, highcut, fs, order=order)
        y = lfilter(b, a, data)
        return y



#-------------------------------------------
#Data Analysis
#-------------------------------------------  
class Constants:
    @staticmethod
    def get_tidal_frequencies_in_hours(tide_type = "semidiurnal"):

        #tidal_constit = {"M2": 12.42, "S2": 12.00, "N2":12.66}
        #the diurnal constituents, K1, O1, P1, Q1, and S1, with periods of 23.93, 25.82, 24.07, 26.87, and 24.00 h, 
        #semidiurnal constituents M2, S2, N2, and S2, with periods of 12.42, 12.00, 12.66, and 11.97 h, 
        if tide_type == "Padman":
            return {"M2": 12.42, "S2": 12.00, "N2":12.66, "K2": 11.96, "K1": 23.93, "O1": 25.82, "P1":24.07, "Q1":26.87, "Mf": 327.86, "Mm": 661.31}

        elif tide_type == "semidiurnal":
             return {"M2": 12.42, "S2": 12.00, "N2":12.66, "K2": 11.96}            
        
        else:
            tidal_constit = {"M2": 12.42, "K1": 23.93, "S2": 12.00, "O1": 25.82, "P1":24.07, "N2":12.66} #, "K2"}       
            return tidal_constit
        
        raise NotImplementedError
        
        
    @staticmethod
    def get_coriolis_frequency(lat, unit = "rad/s", absolute = False):
        
        if absolute == True: lat = abs(lat)
        
        Omega = 7.2921*1e-5 #unit: rad/s
        
        f = 2 * Omega * np.sin(lat * np.pi / 180.)
        
        if unit == "rad/s":
            return f
        
        if unit == "Hz":
            return f / (2*np.pi)

        if unit == "cpd":
            return f / (2*np.pi) *3600 * 24
                  
        raise NotImplementedError("Specified unit is not implemented")
    
    @staticmethod
    def cut_to_f_N_interval(freq,psd,f,N):
        f_index = np.argmin(np.abs(freq-f))
        N_index = np.argmin(np.abs(freq-N))
        
        fN_freq = freq[f_index:N_index]
        fN_psd = psd[f_index:N_index]
        
        if np.shape(fN_freq) == np.shape([]):
            raise AssertionError
        
        return fN_freq, fN_psd
    
        
def tidal_analysis(time,u,v,lat, out_style = "classic"):

    import ttide as tt
    
    complex_uv = u + 1j * v
    days,hours,minutes,seconds = timedelta_to_ints(np.mean(np.diff(time)))
    sample_period = hours + minutes/60 + seconds/3600 #in unit of hours
    
    errcalc = "cboot"
    #errcalc has to be changed to allow a dt > 1
    if sample_period > 1:
        errcalc = "wboot"

    out  = tt.t_tide(xin = complex_uv, dt = sample_period, out_style = out_style, lat = lat, errcalc = errcalc)

    xout = out["xout"]

    tide_u = np.real(xout.squeeze())
    #u_res = u-np.real(xout.squeeze())

    tide_v = np.imag(xout.squeeze())
    #v_res = v-np.imag(xout.squeeze())

    return tide_u, tide_v


def get_bathymetry_data(path = "/media/sf_VM_Folder/data/bathymetry.mat"):
    data = sio.loadmat(path) 
    #print(data.keys())
    #print("\n")
    lat = data["lat"]
    lon = data["lon"]
    z = data["z"]

    return lat,lon,z

class Spectrum:


    @staticmethod
    def total_multitaper(complex_velocity,dt = 1/12,P=10):
        """
        
        In:
          complex_velocity [m/s] technically the unit is arbitrary, the psd will just reflect the original unit
          dt = 1/12 [days]       time duration between measurements, unit is chosen to produce cpd frequency
          P = 10                 time_bandwidth_product, determines the smoothing
        
        Out:
          freq [cpd]
          total_psd [m$^2$/s$^2$ days]
        """
        import scipy.signal as sg  #Package for signal analysis
        import spectrum
    

        f, _S = sg.periodogram(complex_velocity-np.mean(complex_velocity), fs=1/dt, return_onesided=False) #fs = sampling frequency (cyclic)
        psi, eigs = spectrum.mtm.dpss(np.size(complex_velocity), NW=P, k=2*P-1)
        Zk, _weights, _eigenvalues = spectrum.mtm.pmtm(complex_velocity-np.mean(complex_velocity), k=2*P-1,  NFFT=np.size(complex_velocity), v=psi, e=eigs, method="unity");
        S=np.mean(np.transpose(np.abs(Zk)**2), axis=1) * dt
        freq = f[np.where(f>0)]
        #add positive and negative frequencies to a total PSD
        #negative side has to be reversed to align with the positive side
        try:
            total_psd = S[np.where(f<0)][::-1] + S[np.where(f>0)]
        #if len(S) is odd, the last data point is removed
        except ValueError:
            total_psd = S[np.where(f<0)][::-1][:-1] + S[np.where(f>0)]
        return freq, total_psd       
    
    """
    @static_method
    def rotary_multitaper(cvs, dt = 1/12, labels = None, ax = None, **Kwargs):    
        #TODO TODO 
        if ax == None:
            fig,ax = plt.subplots(1)
            ax.set_xlabel('Frequency (cycles/day)')
            ax.set_ylabel('PSD (cm$^2$/s$^2$ days)')
        
        if not isinstance(cvs, list):
            cvs = [cvs]
        
        for i,cv in enumerate(cvs):
        
            cv = cv[~np.isnan(cv)] #remove 

            f, _S = sg.periodogram(cv-np.mean(cv), fs=1/dt) #fs = sampling frequency (cyclic)

            P = 10
            psi, eigs = spectrum.mtm.dpss(np.size(cv), NW=P, k=2*P-1)

            Zk, weights, eigenvalues = spectrum.mtm.pmtm(cv-np.mean(cv), k=2*P-1,  NFFT=np.size(cv), v=psi, e=eigs, method="unity");
            S=np.mean(np.transpose(np.abs(Zk)**2), axis=1) * dt


            if labels is not None:
                if type(labels[i]) not in [str,int]: label = f"{labels[i]:.0f}"
                else: label = str(labels[i])
                ax.loglog(f[np.where(f>0)],S[np.where(f>0)], label = label, **Kwargs)
            else:
                ax.loglog(f[np.where(f>0)],S[np.where(f>0)], **Kwargs)

            
        return f[np.where(f>0)],S[np.where(f>0)]
    """
     
    @staticmethod 
    def integrate_psd_interval(freq,psd,a = None, b = None):
        """
        Integration between a und b using the trapezoidal integration method
        """
        if a == None: a = np.min(freq)
        if b == None: b = np.max(freq)
        
        lower = np.argmin(np.abs(freq-a)).astype(int)
        upper = np.argmin(np.abs(freq-b)).astype(int)
        return np.trapz(y = psd[lower:upper], x = freq[lower:upper]) 



    
                        
#-------------------------------------------
#Plotting
#-------------------------------------------     

class Plot:
    @staticmethod
    def path_as_footnote(fig = None,path = "", rot = "vertical"):
        """
        Add the script origin as a footnote below or to the right side of a figure
        """        
    
        if fig == None:
            fig = plt.gcf()

        if path == None:
            import os
            try:
                path = os.path.realpath(__file__)
            except NameError: 
                # if run in Jupyter Notebook __file__ is not defined 
                # and the path is only the passed string
                path = os.path.abspath('')
        
        if rot == "vertical":
            plt.figtext(0.99, 0.5, f"{path}", horizontalalignment='right', va = "center",\
                figure = fig, rotation = "vertical", fontsize = 6)
        else:
            plt.figtext(0.99, 0, f"{path}", horizontalalignment='right', va = "bottom",\
                figure = fig, fontsize = 6)

    @staticmethod
    def add_tidal_and_f_ticks(ax,f_cpd, lower_bound = 1e-6, upper_bound = 1e-3):
        freqs = np.array(
            [24 / (14 * 24), 24 / 12.4, 2 * 24 / 12.4, 4 * 24 / 12.4, f_cpd, 2 * f_cpd, 1]
        )
        freq_labels = ["14 days", "M2", "2M2", "4M2", " \nf", " \n2f", "K1"]
        
        axis = ax.get_ylim()
        for freq in freqs:
            if freq == 24 / (14 * 24): continue
            ax.vlines(freq, lower_bound, upper_bound, color="k", alpha=0.5, linestyle="-", linewidth=0.75)
        ax.set_ylim(axis)
        
        ax2 = ax.secondary_xaxis(location="bottom")
        ax2 = Plot._axstyle(ax2, ticks="in", grid=False, spine_offset=40) 
        ax2.xaxis.set_ticks([])
        ax2.xaxis.set_ticklabels([])
        ax2.minorticks_off()
        ax2.xaxis.set_ticks(freqs)
        ax2.xaxis.set_ticklabels(freq_labels)

        return ax2

        
    @staticmethod    
    def _axstyle(
        ax=None,
        fontsize=9,
        nospine=False,
        grid=True,
        ticks="off",
        ticklength=2,
        spine_offset=5,
    ):
        """
        Apply own style to axis.

        Parameters
        ----------
        ax : AxesSubplot (optional)
            Current axis will be chosen if no axis provided

        Returns
        -------
        ax : AxesSubplot
            Axis handle
        """
        # find out background color - if this is set to ayu dark, adjust some axis
        # colors
        figcolor = plt.rcParams["figure.facecolor"]
        dark = True if figcolor == "#0d1318" else False

        if ax is None:
            ax = plt.gca()

        # Remove top and right axes lines ("spines")
        spines_to_remove = ["top", "right"]
        for spine in spines_to_remove:
            ax.spines[spine].set_visible(False)

        # Remove bottom and left spines as well if desired
        if nospine:
            more_spines_to_remove = ["bottom", "left"]
            for spine in more_spines_to_remove:
                ax.spines[spine].set_visible(False)

        # For remaining spines, thin out their line and change
        # the black to a slightly off-black dark grey
        almost_black = "#262626"
        # if figure background is dark, set this close to white
        if dark:
            almost_black = "#ebe6d7"

        if ticks == "off":
            # Change the labels to the off-black
            ax.tick_params(
                axis="both",
                which="major",
                labelsize=fontsize,
                colors=almost_black,
            )
            # Get rid of ticks.
            ax.xaxis.set_ticks_position("none")
            ax.yaxis.set_ticks_position("none")
        elif ticks == "in":
            # Change the labels to the off-black
            ax.tick_params(
                axis="both",
                which="major",
                labelsize=fontsize,
                colors=almost_black,
                direction="in",
                length=ticklength,
            )

        spines_to_keep = ["bottom", "left"]
        for spine in spines_to_keep:
            ax.spines[spine].set_linewidth(0.5)
            ax.spines[spine].set_color(almost_black)
            ax.spines[spine].set_position(("outward", spine_offset))

        # Change the labels to the off-black
        ax.yaxis.label.set_color(almost_black)
        ax.yaxis.label.set_size(fontsize)
        ax.yaxis.offsetText.set_fontsize(fontsize)
        ax.xaxis.label.set_color(almost_black)
        ax.xaxis.label.set_size(fontsize)
        ax.xaxis.offsetText.set_fontsize(fontsize)

        # Change the axis title to off-black
        ax.title.set_color(almost_black)
        ax.title.set_size(fontsize + 1)

        # turn grid on
        if grid:
            ax.grid(
                which="major",
                axis="both",
                color="0.5",
                linewidth=0.25,
                linestyle="-",
                alpha=0.8,
            )

        # change legend fontsize (if there is one)
        try:
            plt.setp(ax.get_legend().get_texts(), fontsize=fontsize)
        except AttributeError:
            noleg = 1

        return ax        

