import scipy.signal as sg  #Package for signal analysis
import numpy as np

def _total_multitaper(complex_velocity,dt = 1/12,P=10):
    from multitaper import MTSpec  #using German Prieto's multitaper package, https://github.com/gaprieto/multitaper

    fo, So = sg.periodogram(complex_velocity-np.mean(complex_velocity), fs=1/dt) #fs = sampling frequency (cyclic)

    spec = MTSpec(complex_velocity-np.mean(complex_velocity), nw=P, dt=dt, iadapt=0, nfft=len(complex_velocity))
    S = np.ravel(spec.spec)
    f = np.ravel(spec.freq)

    #add positive and negative frequencies to a total PSD
    #negative side has to be reversed to align with the positive side
    try:
        total_psd = S[np.where(f<0)][::-1] + S[np.where(f>0)]
    #if len(S) is odd, the last data point is removed
    except ValueError:
        total_psd = S[np.where(f<0)][::-1][:-1] + S[np.where(f>0)]
    return freq, total_psd   

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
 
def integrate_psd_interval(freq,psd,a = None, b = None):
    """
    Integration between a und b using the trapezoidal integration method
    """
    if a == None: a = np.min(freq)
    if b == None: b = np.max(freq)
    
    lower = np.argmin(np.abs(freq-a)).astype(int)
    upper = np.argmin(np.abs(freq-b)).astype(int)
    return np.trapz(y = psd[lower:upper], x = freq[lower:upper]) 


