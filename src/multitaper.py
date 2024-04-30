import multitaper as mt #Prieto 2022 https://github.com/gaprieto/multitaper
import scipy.signal as sg  #signal analysis

fs = 1/dt

def multitaper(complex_velocity,product_bandwidth, sampling_frequency):
    """


    """
    fo, So = sg.periodogram(cv-np.mean(cv), fs=1/dt) #fs = sampling frequency (cyclic)

    #using somewhat different settings than the default for MTSpec
    #setting nfft, number of points in the FFT, to number of data points, so no zero-padding
    #setting adapt to 1 for standard average multitaper spectrum for illustrative purposes — 
    #    though in general I do prefer the adaptive version as standard
    P = 16
    spec = mt.MTSpec(cv-np.mean(cv), nw=P, dt=dt, iadapt=1, nfft=len(cv))
    #vars(spec)  #use this to inspect the spec object output by MTSpec
    S = np.ravel(spec.spec)
    f = np.ravel(spec.freq)


    #verify variance is approximately recovered
    #use numpy is_approximately equal
    #TODO assert f[1]-f[0])*np.sum(S) is close to np.var(cv)

    #print((f[1]-f[0])*np.sum(S)) #verify variance is approximately recovered
    raise NotImplementedError
    
def get_total_psd(S,f):
    """
    add positive and negative frequencies to a total PSD
    """
    

    if len(S) % 2 == 0:
        #negative side has to be reversed to align with the positive side
        total_psd = S[np.where(f<0)][::-1] + S[np.where(f>0)]
   
    else: #if len(S) is odd,
        #the last data point of the negative side is removed
        total_psd = S[np.where(f<0)][::-1][:-1] + S[np.where(f>0)]
        
    freq = f[np.where(f>0)]
            
    return freq, total_psd    




#define a function for computing confidence intervals.

def mconf(K,gamma,str):
    """
    Compute symmetric confidence intervals for multitaper spectral estimation.
    
    Args: 
        K: Number of tapers used in the spectral estimate, normally 2*P-1
        gamma: confidence level, e.g., 0.95
        str: 'lin' to return confidence intervals for linear axis 
             'log' to return confidence intervals for log10 axis 
    
    Returns:
        ra,rb:  Ratio factors for confidence interval
        
    If S0 is the true value of the spectrum and S is the spectral estimate,
    then the confindence interval is defined such that 
    
        Probability that ra < S/S0 < ra = gamma     (linear case)
        Probability that ra < log10(S)/log10(S0) < ra = gamma (log10 case)

    The confidence interval will be S*ra to S*rb for the spectral values
    in linear space, or S*10^ra to S*10^rb in log10 space.  If log10(S) 
    is plotted rather than S with a logarithmic axis, the latter would 
    become log10(S)+ra to log10(S)*rb.    
    
    taken from http://jmlilly.net/course/labs/html/SpectralAnalysis-Python.html
    """

    dx=0.0001
    #Compute pdf symmetrically about unity

    if str=='lin':
    
        x1=np.arange(1-dx/2,-1,step=-dx)*2*K  #from one to minus one
        x2=np.arange(1+dx/2,3,step=dx)*2*K #from 1 to three 
        fx=(chi2.pdf(x1,2*K)+chi2.pdf(x2,2*K))*2*K
        sumfx=np.cumsum(fx)*dx

        ii=np.where(sumfx>=gamma)[0][0]
        ra=x1[ii]/2/K
        rb=x2[ii]/2/K
    
    elif str=='log':
        
        xo=np.log(2)+np.log(1/2/K)+digamma(K) #see pav15-arxiv
        c=np.log(10)
        xo=xo/c  #change of base rule
    
        x1=np.power(10,np.arange(xo-dx/2,xo-2,step=-dx))   
        x2=np.power(10,np.arange(xo+dx/2,xo+2,step=dx))  
    
        fx1=c*2*K*np.multiply(x1,chi2.pdf(2*K*x1,2*K))
        fx2=c*2*K*np.multiply(x2,chi2.pdf(2*K*x2,2*K))

        sumfx=np.cumsum(fx1+fx2)*dx
    
        ii=np.where(sumfx>=gamma)[0][0]
        ra=np.log10(x1[ii])
        rb=np.log10(x2[ii])
    
    return ra,rb    



###################################
# Previous Code
###################################
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
