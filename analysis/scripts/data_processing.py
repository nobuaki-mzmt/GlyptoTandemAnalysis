"""
data_preparation.py
N. Mizumoto
This script reads all .h5 results from SLEAP and organize for the further analysis
Gly sat
0: Head 1: Pronotum 2: Tip 3: AntennaR 4: AntennaL 5: BodyCenter 6: Marker
"""

import sys
import os
import glob

import h5py

import numpy as np
from numpy.linalg import norm

import pandas as pd

import scipy
from scipy.interpolate import interp1d

import math
import feather

#------------------------------------------------------------------------------#
# interpolate the data
#------------------------------------------------------------------------------#
def fill_missing(Y, kind="linear"):
    initial_shape = Y.shape
    Y = Y.reshape((initial_shape[0], -1))
    # Interpolate along each slice.
    for i in range(Y.shape[-1]):
        y = Y[:, i]
        # Build interpolant.
        x = np.flatnonzero(~np.isnan(y))
        if len(x) > 3:
          f = interp1d(x, y[x], kind=kind, fill_value=np.nan, bounds_error=False)
          # Fill missing
          xq = np.flatnonzero(np.isnan(y))
          y[xq] = f(xq)
          # Fill leading or trailing NaNs with the nearest non-NaN values
          mask = np.isnan(y)
          y[mask] = np.interp(np.flatnonzero(mask), np.flatnonzero(~mask), y[~mask])
          Y[:, i] = y
          if sum(np.isnan(y)) > 0:
            print("error"+str(i))
            print("error"+(i))
    # Restore to initial shape.
    Y = Y.reshape(initial_shape)
    return Y
#------------------------------------------------------------------------------#

#------------------------------------------------------------------------------#
#------------------------------------------------------------------------------#
def data_preparation(in_dir, out_dir = "data/"):
  files = glob.glob(in_dir + "/*.h5")
  df_ind = pd.DataFrame()
  df_dis = pd.DataFrame()
  
  for f_name in files:
    
    ## load data
    with h5py.File(f_name, "r") as f:
      dset_names = list(f.keys())
      locations = f["tracks"][:].T
      node_names = [n.decode() for n in f["node_names"][:]]
    
    if locations.shape[3] > 2:
      locations = locations[:,:,:,0:2]
    
    print(f_name)
    pair_name = os.path.basename(f_name.replace(".h5", ""))
    
    ## processing locations
    # data filling
    locations = fill_missing(locations)
    
    # filtering
    for i_ind in range(locations.shape[3]):
      for i_coord in range(locations.shape[2]):
        for i_nodes in range(locations.shape[1]):
          locations[:, i_nodes, i_coord, i_ind] = scipy.signal.medfilt( locations[:, i_nodes, i_coord, i_ind], 5)
    
    # obtain meta data
    species = pair_name.split("_")[0]
    colony  = pair_name.split("_")[1]
    
    # measure body length
    body_length = []
    for ind_i in range(locations.shape[3]):
      body_length_temp = ( 
        np.sqrt( (locations[:, 0, 0, ind_i] - locations[:, 1, 0, ind_i])**2 +
              (locations[:, 0, 1, ind_i] - locations[:, 1, 1, ind_i])**2) +
        np.sqrt( (locations[:, 1, 0, ind_i] - locations[:, 2, 0, ind_i])**2 +
              (locations[:, 1, 1, ind_i] - locations[:, 2, 1, ind_i])**2) 
      )
      body_length = np.append(body_length, np.mean(body_length_temp))
    body_length_mean = np.mean(body_length)
    
    
    ## generate dataframes
    ind_name = ["f","m"]
    for i_ind in range(locations.shape[3]):
      df_temp = {
        "frame":  list(range(0,locations.shape[0],1)),
        "head_x": locations[:, 0, 0, i_ind],
        "head_y": locations[:, 0, 1, i_ind],
        "pron_x": locations[:, 1, 0, i_ind],
        "pron_x": locations[:, 1, 1, i_ind],
        "tip_x": locations[:, 2, 0, i_ind],
        "tip_x": locations[:, 2, 1, i_ind],
        "sex": ind_name[i_ind],
        "video": pair_name
        }
      df_temp = pd.DataFrame(df_temp)
      df_ind = pd.concat([df_ind, pd.DataFrame(df_temp)])
    
    # dis between female-male
    fhead_mhead_dis = np.sqrt( (locations[:, 0, 0, 1] - locations[:, 0, 0, 0])**2 + (locations[:, 0, 1, 1] - locations[:, 0, 1, 0])**2 ) / body_length_mean
    fhead_mpron_dis = np.sqrt( (locations[:, 1, 0, 1] - locations[:, 0, 0, 0])**2 + (locations[:, 1, 1, 1] - locations[:, 0, 1, 0])**2 ) / body_length_mean
    fhead_mtip_dis  = np.sqrt( (locations[:, 2, 0, 1] - locations[:, 0, 0, 0])**2 + (locations[:, 2, 1, 1] - locations[:, 0, 1, 0])**2 ) / body_length_mean
    
    fpron_mhead_dis = np.sqrt( (locations[:, 0, 0, 1] - locations[:, 1, 0, 0])**2 + (locations[:, 0, 1, 1] - locations[:, 1, 1, 0])**2 ) / body_length_mean
    fpron_mpron_dis = np.sqrt( (locations[:, 1, 0, 1] - locations[:, 1, 0, 0])**2 + (locations[:, 1, 1, 1] - locations[:, 1, 1, 0])**2 ) / body_length_mean
    fpron_mtip_dis  = np.sqrt( (locations[:, 2, 0, 1] - locations[:, 1, 0, 0])**2 + (locations[:, 2, 1, 1] - locations[:, 1, 1, 0])**2 ) / body_length_mean
      
    ftip_mhead_dis  = np.sqrt( (locations[:, 0, 0, 1] - locations[:, 2, 0, 0])**2 + (locations[:, 0, 1, 1] - locations[:, 2, 1, 0])**2 ) / body_length_mean
    ftip_mpron_dis  = np.sqrt( (locations[:, 1, 0, 1] - locations[:, 2, 0, 0])**2 + (locations[:, 1, 1, 1] - locations[:, 2, 1, 0])**2 ) / body_length_mean
    ftip_mtip_dis   = np.sqrt( (locations[:, 2, 0, 1] - locations[:, 2, 0, 0])**2 + (locations[:, 2, 1, 1] - locations[:, 2, 1, 0])**2 ) / body_length_mean
    
    df_temp = {
        "frame":  list(range(0,locations.shape[0],1)),
        "fhead_mhead_dis": fhead_mhead_dis,
        "fhead_mpron_dis": fhead_mpron_dis,
        "fhead_mtip_dis": fhead_mtip_dis,
        "fpron_mhead_dis": fpron_mhead_dis,
        "fpron_mpron_dis": fpron_mpron_dis,
        "fpron_mtip_dis": fpron_mtip_dis,
        "ftip_mhead_dis": ftip_mhead_dis,
        "ftip_mpron_dis": ftip_mpron_dis,
        "ftip_mtip_dis": ftip_mtip_dis,
        "video": pair_name
        }
    df_temp = pd.DataFrame(df_temp)
    df_dis = pd.concat([df_dis, pd.DataFrame(df_temp)])
    
  filename = os.path.basename(in_dir) + "_df_ind.h5"
  #df_ind.to_hdf(out_dir + filename, key = "df_ind", mode='w',  complevel=9, complib='zlib')
  feather.write_dataframe(df_ind, filename, compression='zstd')

  filename = os.path.basename(in_dir) + "_df_dis.h5"
  #df_dis.to_hdf(out_dir + filename, key = "df_dis", mode='w',  complevel=9, complib='zlib')
  feather.write_dataframe(df_dis, filename, compression='zstd')
#------------------------------------------------------------------------------#

#------------------------------------------------------------------------------#
def main():
  data_place = glob.glob("data/raw/*")
  for data_place_i in data_place:
    print(data_place_i)
    data_preparation(in_dir = data_place_i)
#------------------------------------------------------------------------------#

#------------------------------------------------------------------------------#
if __name__ == "__main__":
    main()
#------------------------------------------------------------------------------#


df_dis = feather.read_dataframe("data/Gly_sat_df_dis.feather")
feather.write_dataframe(df_dis, "data/Gly_sat_df_dis.feather")
