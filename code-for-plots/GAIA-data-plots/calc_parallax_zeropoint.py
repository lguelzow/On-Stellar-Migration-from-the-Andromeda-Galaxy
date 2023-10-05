import numpy as np
import pandas as pd
import argparse
from zero_point import zpt # install with "pip install gaiadr3-zeropoint"

# start command line argument parser
parser = argparse.ArgumentParser(description='')

# parser asks for path to input file in .csv format
parser.add_argument('path', metavar='PATH', type=str, nargs='*', default=[],
                    help='Choose csv input file with Gaia photometric \
                information in addition to full-phase-space solutions.')

# put command line arguments into convenient variable
args = parser.parse_args()

# load coefficient tables
zpt.load_tables()

# read data from input file
data = pd.read_csv(args.path)

# calculate zero-point parallax with panda wrapper
zero_point = data.apply(zpt.zpt_wrapper,axis=1)

# add zero-point parallax column to data frame
data['parallax_zpt'] = zero_point

# drop unnecessary data for further analysis
data.drop('phot_g_mean_mag', axis=1, inplace=True)
data.drop('nu_eff_used_in_astrometry', axis=1, inplace=True)
data.drop('pseudocolour', axis=1, inplace=True)
data.drop('ecl_lat', axis=1, inplace=True)
data.drop('astrometric_params_solved', axis=1, inplace=True)

# write new data column into new file
data.to_csv('Gaia_data_with_zpt_parallax.csv', encoding='utf-8')