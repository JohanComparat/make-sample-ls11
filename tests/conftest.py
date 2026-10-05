import numpy as np
import pytest

from ls11samples.config import load_config
from ls11samples.env import get_paths


@pytest.fixture
def cfg():
    return load_config()


@pytest.fixture
def rng():
    return np.random.default_rng(42)


@pytest.fixture
def first_sweep():
    sweeps = get_paths().sweeps()
    if not sweeps:
        pytest.skip("no DR11 sweep available locally")
    return sweeps[0]
