import pytest
import pytest-regressions
import pandas as pd

def test_thorpe_output():
  data_frame = pd.DataFrame.from_dict({
      'U_gas': U[0][positions],
      'U_liquid': U[1][positions],
      'gas_vol_frac [-]': vol_frac[0][positions],
      'liquid_vol_frac [-]': vol_frac[1][positions],
      'P': Pa_to_bar(P)[positions],
  })
  
  dataframe_regression.check(data_frame)

def test_finestructure_output():
    pass
  
def test_idemix_output():
    pass


  
