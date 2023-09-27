"""
error check for Gly fus
"""
import os
import sys
import glob
import h5py

import math
pi = math.pi

import numpy as np
from numpy.linalg import norm
from scipy.interpolate import interp1d
import pandas as pd

import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib.animation as animation

import tkinter as tk
from tkinter import filedialog

#------------------------------------------------------------------------------#
#------------------------------------------------------------------------------#
def fill_missing(Y, kind="linear"):
    initial_shape = Y.shape
    Y = Y.reshape((initial_shape[0], -1))
    # Interpolate along each slice.
    for i in range(Y.shape[-1]):
        y = Y[:, i]
        # Build interpolant.
        x = np.flatnonzero(~np.isnan(y))
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
def data_preparation(in_dir, dish_size_mm, dish_size_pixel):
  files = glob.glob(in_dir + "/*.h5")
  df = pd.DataFrame()
  for f_name in files:

    with h5py.File(f_name, "r") as f:
      dset_names = list(f.keys())
      locations = f["tracks"][:].T
      node_names = [n.decode() for n in f["node_names"][:]]

    print(f_name)
    
    # data filling
    locations = fill_missing(locations)
    
    # scaling in mm (2000 pixels = 55 mm)
    locations[:, :, :, :] = locations[:, :, :, :] / dish_size_pixel * dish_size_mm
    
    ind_name = ["f", "m"]
    video_name = os.path.basename(f_name).split(".")[0]
    
    error_count = 0
    for ind in range(locations.shape[3]):
      # too high speed movement
      head_speed = np.sqrt( np.diff(locations[:, 0, 0, ind])**2 + np.diff(locations[:, 0, 1, ind])**2 )
      tip_speed = np.sqrt( np.diff(locations[:, 2, 0, ind])**2 + np.diff(locations[:, 2, 1, ind])**2 )
      #plt.figure()
      #plt.hist(tip_speed, bins=200, edgecolor='black')
      #plt.show()
      df_temp = {
        "video": video_name,
        "frame": list(np.where( (tip_speed > 2) & (head_speed > 2))[0] + 2),
        "sex": ind_name[ind],
        "reason": "too fast movement"
        }
      df_temp = pd.DataFrame(df_temp)
      df = pd.concat([df, pd.DataFrame(df_temp)])
      error_count += len(df_temp)
      
      # head tip dis is too small
      head_tip_dis_self = np.sqrt ( (locations[:, 0, 0, ind] - locations[:, 2, 0, ind])**2 +
        (locations[:, 0, 1, ind] - locations[:, 2, 1, ind])**2 )
      #plt.figure()
      #plt.hist(head_tip_dis_self, bins=200, edgecolor='black')
      #plt.show()
      df_temp = {
        "video": video_name,
        "frame": list(np.where( head_tip_dis_self < 2)[0] + 1),
        "sex": ind_name[ind],
        "reason": "abnormal body length"
        }
      df_temp = pd.DataFrame(df_temp)
      df = pd.concat([df, pd.DataFrame(df_temp)])
      error_count += len(df_temp)
      
      # body boundary error
      head_selftip_dis = np.sqrt ( (locations[:, 0, 0, ind] - locations[:, 2, 0, 1-ind])**2 +
        (locations[:, 0, 1, ind] - locations[:, 2, 1, 1-ind])**2 )
      df_temp = {
        "video": video_name,
        "frame": list(np.where( (tip_speed <1 ) & (head_speed > 2) & (head_selftip_dis[1:] < 5))[0] + 2),
        "sex": ind_name[ind],
        "reason": "likely body fusion"
        }
      df_temp = pd.DataFrame(df_temp)
      df = pd.concat([df, pd.DataFrame(df_temp)])
      error_count += len(df_temp)
    
    print(str(error_count) + " error found")
  df.to_csv(in_dir + "/error_frame.csv", index=False)
    

#------------------------------------------------------------------------------#

#------------------------------------------------------------------------------#
def main():
  root = tk.Tk()
  root.withdraw()  # Hide the main window

  directory_path = filedialog.askdirectory(title="Select a Directory")

  if directory_path:
    print("Directory exists")
    data_preparation(directory_path, 90, 1200)
  else:
    print("no Directory")
#------------------------------------------------------------------------------#

#------------------------------------------------------------------------------#
if __name__ == "__main__":
  main()
#------------------------------------------------------------------------------#
