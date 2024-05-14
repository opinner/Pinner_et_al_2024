#from dataclasses import dataclass
#import numpy as np

"""
@dataclass
class CTDCast():
    lat: float 
    lon: float 
    #pressure: np.ndarray
    #in_situ_temperature: np.ndarray
    #Practical_Salinity: np.ndarray
"""
    
import pandas as pd

#@pd.api.extensions.register_dataframe_accessor("oze")
class CTDCast(pd.DataFrame):

    # normal properties
    _metadata = ["location","name", "date"]

    @property
    def _constructor(self):
        return CTDCast
       
     

