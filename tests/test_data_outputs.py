import pytest_regressions
import pandas as pd
# tests assumes the default tolerance level of numpy.isclose


def test_thorpe_dissipation_rate_output(dataframe_regression):
    thorpe_eps_df = pd.read_pickle("./scripts/thorpe_scales/method_results/Thorpe_eps_df_with_mab.pkl")
    dataframe_regression.check(thorpe_eps_df, default_tolerance=dict(rtol=1e-05, atol=1e-10))


def test_thorpe_dissipation_rate_derived_data(dataframe_regression):
    thorpe_eps_df = pd.read_csv("./derived_data/binned_thorpe_dissipation.csv")
    dataframe_regression.check(thorpe_eps_df, default_tolerance=dict(rtol=1e-05, atol=1e-10))


def test_finestructure_dissipation_rate_output(dataframe_regression):
    eps_IGW_strain_df = pd.read_csv("./scripts/finestructure/method_results/binned_strain_eps.csv")
    dataframe_regression.check(eps_IGW_strain_df, default_tolerance=dict(rtol=1e-05, atol=1e-10))


def test_finestructure_dissipation_derived_data(dataframe_regression):
    eps_IGW_strain_df = pd.read_csv("./derived_data/binned_finestructure_dissipation.csv")
    dataframe_regression.check(eps_IGW_strain_df, default_tolerance=dict(rtol=1e-05, atol=1e-10))


def test_idemix_dissipation_rate_output(dataframe_regression):
    eps_IGW_IDEMIX_df = pd.read_csv("./scripts/IDEMIX_parameterization/method_results/eps_IGW_IDEMIX_results.csv")
    dataframe_regression.check(eps_IGW_IDEMIX_df, default_tolerance=dict(rtol=1e-05, atol=1e-10))


def test_idemix_dissipation_rate_derived_data(dataframe_regression):
    eps_IGW_IDEMIX_df = pd.read_csv("./derived_data/wave_energy_dissipation.csv")
    dataframe_regression.check(eps_IGW_IDEMIX_df, default_tolerance=dict(rtol=1e-05, atol=1e-10))


def test_neutral_density_output(dataframe_regression):
    binned_neutral_density_df = pd.read_csv("./scripts/preprocessing/method_results/binned_gamma_n.csv")
    dataframe_regression.check(binned_neutral_density_df, default_tolerance=dict(rtol=1e-05, atol=1e-08))


def test_neutral_density_derived_data(dataframe_regression):
    binned_neutral_density_df = pd.read_csv("./derived_data/binned_neutral_density.csv")
    dataframe_regression.check(binned_neutral_density_df, default_tolerance=dict(rtol=1e-05, atol=1e-08))

