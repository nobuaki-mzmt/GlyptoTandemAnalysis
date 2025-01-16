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
def data_filter(in_dir, dish_size, species):
  files = glob.glob(in_dir + "/*.h5")
  df = pd.DataFrame()
  df_bl_al = pd.DataFrame()
  skip_list = []
  for f_name in files:
    print(f_name)
    ## load data
    with h5py.File(f_name, "r") as f:
      try:
        locations = f['tracks'][:]
        node_names = [n.decode() for n in f["node_names"][:]]
        print("'tracks' exists")
      except KeyError:
        print("'tracks' does not exist")
        skip_list.append(f_name)
        continue;
    
    locations = locations.T
    if 'BodyCenter' in node_names:
        targets = ['Head', 'BodyCenter', 'Tip']
    elif 'Middle' in node_names:
        targets = ['Head', 'Middle', 'Tip']
    indices = [node_names.index(target) for target in targets if target in node_names]
    locations = locations[:,indices,:,:]
    
    pair_name = os.path.basename(f_name.replace(".h5", ""))
    total_frame = locations.shape[0]

    if locations.shape[3] > 2:
      locations = locations[:,:,:,0:2]
  

    ## processing locations
    # data filling
    locations = fill_missing(locations)
    
    # scaling in mm (2000 pixels = dish_size)
    locations[:, :, :, :] = locations[:, :, :, :] / 1200 * dish_size
    
    # filtering
    for i_ind in range(locations.shape[3]):
      for i_coord in range(locations.shape[2]):
        for i_nodes in range(locations.shape[1]):
          locations[:, i_nodes, i_coord, i_ind] = scipy.signal.medfilt( locations[:, i_nodes, i_coord, i_ind], 5)
    
    # measure body length and antenna length
    body_length = []
    for ind_i in range(locations.shape[3]):
      body_length_temp = ( 
        np.sqrt( (locations[:, 0, 0, ind_i] - locations[:, 1, 0, ind_i])**2 +
              (locations[:, 0, 1, ind_i] - locations[:, 1, 1, ind_i])**2) +
        np.sqrt( (locations[:, 1, 0, ind_i] - locations[:, 2, 0, ind_i])**2 +
              (locations[:, 1, 1, ind_i] - locations[:, 2, 1, ind_i])**2)  
      )
      body_length = np.append(body_length, np.mean(body_length_temp))
      
    df_temp = {
        "video":       pair_name,
        "species":     species,
        "total_frame": total_frame,
        "body_length_female": body_length[0],
        "body_length_male": body_length[1]
        }
        
    df_bl_al = pd.concat([df_bl_al, pd.DataFrame([df_temp])])
    
    # summarize data
    mTip_fHead_dis = np.sqrt( (locations[:, 0, 0, 0] - locations[:, 2, 0, 1])**2 + (locations[:, 0, 1, 0] - locations[:, 2, 1, 1])**2 )
    fTip_mHead_dis = np.sqrt( (locations[:, 2, 0, 1] - locations[:, 0, 0, 0])**2 + (locations[:, 2, 1, 1] - locations[:, 0, 1, 0])**2 )
    
    df_temp = {
        "frame": list(range(0,locations.shape[0],1)),
        "fHead_x": locations[:, 0, 0, 0].round(2),
        "fHead_y": locations[:, 0, 1, 0].round(2),
        "mHead_x": locations[:, 0, 0, 1].round(2),
        "mHead_y": locations[:, 0, 1, 1].round(2),
        "fTip_x": locations[:, 2, 0, 0].round(2),
        "fTip_y": locations[:, 2, 1, 0].round(2),
        "mTip_x": locations[:, 2, 0, 1].round(2),
        "mTip_y": locations[:, 2, 1, 1].round(2),
        "fCenter_x": locations[:, 1, 0, 0].round(2),
        "fCenter_y": locations[:, 1, 1, 0].round(2),
        "mCenter_x": locations[:, 1, 0, 1].round(2),
        "mCenter_y": locations[:, 1, 1, 1].round(2),
        "video": pair_name
      }
    df_temp = pd.DataFrame(df_temp)
    
    df = pd.concat([df, pd.DataFrame(df_temp)])
    
  filename = os.path.basename(in_dir) + "_df.feather"
  df.reset_index().to_feather("data_fmt/" + filename)
  print("skipped item: " + str(skip_list))
  return 0
#------------------------------------------------------------------------------#


#------------------------------------------------------------------------------#
def main():
  place = "data_raw/*"
  data_place_species = glob.glob(place)
  print(data_place_species)
  for data_place_species_i in data_place_species:
    print(data_place_species_i)
    species = os.path.basename(data_place_species_i)
    dish_size = 90
    data_filter(in_dir = data_place_species_i, dish_size = dish_size, species = species)
    
#------------------------------------------------------------------------------#

#------------------------------------------------------------------------------#
if __name__ == "__main__":
  # Call the main function without any parameter
  main()
#------------------------------------------------------------------------------#


