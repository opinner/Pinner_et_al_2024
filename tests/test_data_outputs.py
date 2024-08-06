import pytest
# uses pytest-regressions
import pandas as pd

def test_thorpe_output():
    thorpe_eps_df = pd.read_pickle("../scripts/thorpe_scales/method_results/Thorpe_eps_df_with_mab.pkl")
    pytest.dataframe_regression.check(thorpe_eps_df)

def test_finestructure_output():
    eps_IGW_strain_df = pd.read_csv("../scripts/shear_strain_parametrization/method_results/binned_strain_eps.csv")
    pytest.dataframe_regression.check(eps_IGW_strain_df)

def test_idemix_output():
    eps_IGW_IDEMIX_df = pd.read_csv("../scripts/IDEMIX_parametrization/method_results/eps_IGW_IDEMIX_results.csv")
    pytest.dataframe_regression.check(eps_IGW_IDEMIX_df)

def test_neutral_density_output():
    binned_neutral_density_df = pd.read_csv("../data/binned_gamma_n.csv")
    pytest.dataframe_regression.check(binned_neutral_density_df)
  
