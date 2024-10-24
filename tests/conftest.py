# conftest.py
# Try to set the output of pytest to use scientific notation
import pandas as pd
import pytest

@pytest.hookimpl(tryfirst=True)
def pytest_configure(config):
    pd.set_option('display.float_format', '{:.3e}'.format)
