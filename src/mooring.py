import pandas as pd

@pd.api.extensions.register_dataframe_accessor("oze")
class Mooring(pd.DataFrame):

    # normal properties
    _metadata = ["location","time_delta", "max_length"]

    @property
    def _constructor(self):
        return Mooring
       
    def get_shallowest_time_series(self):
        """
        Assumes the names of the columns can be converted to numbers
        """
        shallowest_depth = min([int(col) for col in self.columns[1:]])
        return self[str(shallowest_depth)]
        
    def get_deepest_time_series(self):
        """
        Assumes the names of the columns can be converted to numbers
        """
        deepest_depth = max([int(col) for col in self.columns[1:]])
        return self[str(deepest_depth)]    
