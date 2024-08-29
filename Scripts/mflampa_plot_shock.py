#!/usr/bin/env python3
#  Copyright (C) 2002 Regents of the University of Michigan, 
#  portions used with permission 
#  For more information, see http://csem.engin.umich.edu/tools/swmf
'''
Provide a plotting example for visualizing the shock surface
by MH_data_Shock_* data, written by Weihao Liu
'''

import os
import re
import argparse
import configparser
import numpy as np
import datetime as dt
from scipy.spatial import Delaunay
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d.art3d import Poly3DCollection

# Set default values
default_param_file = 'Param/' + 'PARAM.ini.test.201304.shock'
default_path = 'run_actual13_test_shock2/' + 'SP/IO2/'

# Create an ArgumentParser object
parser = argparse.ArgumentParser(description='Process some paths and parameters.')

# Add an argument for the parameter file
parser.add_argument('--PARAMin', type=str, default=default_param_file,
                    help='Specify the parameter file')
# Add an argument for the path with a default value
parser.add_argument('--PATH', type=str, default=default_path, 
                    help=f'Specify the path (default: {default_path})')

# Parse the arguments
args = parser.parse_args()

# Retrieve the path and parameter file
path = args.PATH
param_file = args.PARAMin

# Get the config file and path from command-line arguments if provided
# param_file = --PARAMin PARAMIN else default_param_file
# path = --PATH PATH else default_path

# Create a ConfigParser object
config = configparser.ConfigParser()
# Load the config file
config.read(param_file)
nFileRead = int(config['MHDATA']['nFileRead'])

# Declare shock coordinate variable number
nShockVar = 7

#==========================================================================================================
def read_line(filename, linenum):
    # Given filename, read No.{linenum} in this file, and return string
    try:
        with open(filename, 'r') as file:
            lines = file.readlines()
            return lines[linenum - 1].strip()
    except FileNotFoundError:
        return f"File {filename} not found."
    except IndexError:
        return f"File {filename} does not have {filename} lines."
    except Exception as e:
        return f"An error occurred: {e}"
#==========================================================================================================
def load_shock_files():
    # Get all the files in the path
    files = sorted(os.listdir(path))
    dic_files_shock = {}
    
    # Filter the files with "MH_data_Shock_e*_*_n*.out"
    files_shock = []; files_time = []
    fmt_file_shock = re.compile(r"^MH_data_Shock_e\d{8}_\d{6}_n\d{6}\.out$")
    for file in files:
        if fmt_file_shock.match(file):
            files_shock.append(file)
            files_time.append(dt.datetime(year=int(file[15:19]), month=int(file[19:21]), \
                day=int(file[21:23]), hour=int(file[24:26]), minute=int(file[26:28])))
    files_shock = np.array(files_shock)
    dic_files_shock['datetime'] = np.array(files_time)
    
    # Read the header
    headerShock = read_line(path+files_shock[0], linenum=5)
    print("Shock file header:", headerShock)
    headerShock = headerShock.split(" ")
    for iHeader in range(len(headerShock)-1, -1, -1):
        if headerShock[iHeader] == "": del headerShock[iHeader]
    dic_files_shock['header'] = headerShock
    
    # Read the files for the shock surface
    print("... Reading in shock files ...")
    for iFile in range(nFileRead):
        CoorShock_II = np.loadtxt(path+files_shock[iFile], skiprows=5, \
            unpack=True, usecols=[i for i in range(nShockVar+1)], dtype=float)
        # Remove those very large values (bugs)
        RShockMedian = np.median(CoorShock_II[5])
        mskEffShockLine = CoorShock_II[5] < 1.8*RShockMedian
        # Get the effective field lines constructing the shock surface
        CoorShockEff_II = CoorShock_II[:, mskEffShockLine]
        dic_files_shock['%d'%(iFile)] = CoorShockEff_II
    
    # Return dictionary
    return dic_files_shock
#==========================================================================================================
def plot_shock_surf(CoorXyz_IB, time, DoSave=False, DoShow=False, figname='Delaunay3d_shock_sphere.pdf'):

    # Perform Delaunay triangulation in 3D
    points = CoorXyz_IB
    mesh = Delaunay(points)

    # Define the tetrahedra
    tetrahedra = mesh.simplices
    # print(tetrahedra)

    # Create a new figure and 3D axis
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    plt.title("3D Shock Surface Triangulation, at %.4d-%.2d-%.2d/%.2d:%.2d" \
        %(time.year, time.month, time.day, time.hour, time.minute))

    # Plot the Delaunay tetrahedra
    for tet in tetrahedra:
        verts = points[tet[:]]
        poly = Poly3DCollection([verts], alpha=0.3, edgecolor='w', linewidths=0.2)
        ax.add_collection3d(poly)

    # Plot the points with colors based on their values
    sc = ax.scatter(points[:, 0], points[:, 1], points[:, 2], marker='o', s=1)  # Reduced size

    # Set labels
    ax.set_xlabel(r'X [$R_\mathrm{sun}$]')
    ax.set_ylabel(r'Y [$R_\mathrm{sun}$]')
    ax.set_zlabel(r'Z [$R_\mathrm{sun}$]')

    # Set equal aspect ratio
    ax.set_box_aspect([1, 1, 1])  # Aspect ratio is 1:1:1

    # Show/save the figure
    if(DoSave): fig.savefig(default_path+"../"+figname)
    if(DoShow): plt.show()
#==========================================================================================================

if __name__ == "__main__":
    # Load the data from simulations
    dic_files_shock = load_shock_files()
    XShock_ = dic_files_shock['header'].index("XShock")
    ZShock_ = dic_files_shock['header'].index("ZShock")

    # Plot and show the shock surface by Delaunay triangulation
    for iFiles in range(60, 61, 1):
        time = dic_files_shock['datetime'][iFiles]
        plot_shock_surf(dic_files_shock['%d'%(iFiles)][XShock_:ZShock_+1, :].T, time, DoShow=True)
