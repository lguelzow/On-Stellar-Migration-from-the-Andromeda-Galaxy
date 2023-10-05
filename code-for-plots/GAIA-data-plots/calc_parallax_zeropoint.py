import matplotlib
import matplotlib.pyplot as plt
from matplotlib.ticker import PercentFormatter

import numpy as np
import pandas as pd

import astropy
import healpy as hp
import sys
import scipy.stats
from astropy.coordinates import SkyCoord
import astropy.units as u
from mpl_toolkits.mplot3d import Axes3D
from zero_point import zpt

# load coefficient tables
zpt.load_tables()

# data = pd.read_csv('./error33_v 150_zpt.csv')
data = pd.read_csv('./zpt_test.csv')

# calculate zero-point parallax with panda wrapper
zero_point = data.apply(zpt.zpt_wrapper,axis=1)

# add zero-point parallax column to data frame
data['parallax_zpt'] = zero_point

data.drop('phot_g_mean_mag', axis=1, inplace=True)
data.drop('nu_eff_used_in_astrometry', axis=1, inplace=True)
data.drop('pseudocolour', axis=1, inplace=True)
data.drop('ecl_lat', axis=1, inplace=True)
data.drop('astrometric_params_solved', axis=1, inplace=True)

data.to_csv('Gaia_data_with_zpt_parallax.csv', encoding='utf-8')