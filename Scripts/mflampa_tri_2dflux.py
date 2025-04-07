#!/usr/bin/env python3
'''
Load and Plot M-FLAMPA 2d spherical triangulated particle flux
Originally Written by Weihao Liu (whliu@umich.edu) in Apr, 2025
'''

import re, os
from io import StringIO
import pandas as pd
import numpy as np
import datetime as dt
import matplotlib.tri as tri
from scipy.ndimage import gaussian_filter

# Read data from some lines in files
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

# Read the MH_data_EChannel.H file
def read_mh_echannel(filepath="./SP/IO2/MH_data_EChannel.H"):
    with open(filepath, 'r') as file:
        lines = file.readlines()

    # Find the line where the header starts
    id_header_line = None
    for i, line in enumerate(lines):
        if re.match(r'\s*FluxChannel\s+ChannelSource', line):
            id_header_line = i
            break
    if id_header_line is None:
        raise ValueError("Table header not found in the file.")

    # Extract only the lines containing the table
    table_lines = lines[id_header_line:]
    table_text = ''.join(table_lines)

    # Use StringIO to treat the text like a file
    table_io = StringIO(table_text)
    # Read the table with whitespace separator
    df = pd.read_csv(table_io, sep=r'\s+', engine='python')

    return df

# Read the MH2D files in the given path
def read_simu_mh2d(pathMH2D="./SP/"):
        
    # mh2d files: first 5 lines are meta data saved as a header
    nLineHeader = 5
    
    # Get all files aligned with the specific mh2d filename format in the given path
    fmt_mh2dfile = re.compile(r"^MH_data_R=\d{4}.\d{2}_e\d{8}_\d{6}_n\d{6}\.out$")
    files_mh2dfile = [ ]
    files = sorted(os.listdir(pathMH2D))
    for j in range(len(files)):
        if fmt_mh2dfile.match(files[j]):
            files_mh2dfile.append(pathMH2D+files[j])
    dic = {}
    
    # Time array
    nIter_MH2D = len(files_mh2dfile)
    time_I = [dt.datetime( \
        year=int(files_mh2dfile[iFile][-27:-23]), month=int(files_mh2dfile[iFile][-23:-21]), \
        day=int(files_mh2dfile[iFile][-21:-19]), hour=int(files_mh2dfile[iFile][-18:-16]), \
        minute=int(files_mh2dfile[iFile][-16:-14]), second=int(files_mh2dfile[iFile][-14:-12])) \
        for iFile in range(nIter_MH2D)]
    dic['Time'] = np.array(time_I)
    
    # Position and flux
    # if there is any mh2d file
    if nIter_MH2D > 0:
        # Get: Numbers of energy channels and field lines
        header_mh2d_I = read_line(files_mh2dfile[0], linenum=nLineHeader).split(" ")
        nEnergymh2d = 0
        for i in range(len(header_mh2d_I)-1, -1, -1):
            if "flux_Channel" in header_mh2d_I[i]: nEnergymh2d += 1
        nLine = np.loadtxt(files_mh2dfile[0], skiprows=nLineHeader, \
            unpack=True, usecols=[0], dtype=float).astype(np.int64)[-1]
        iLon = header_mh2d_I.index("Longitude") # Index of "Longtude"
        iLat = header_mh2d_I.index("Latitude") # Index of "Latitude"
        
        LonLat_III = np.zeros([nIter_MH2D, nLine, 2])
        IFluxOp_III = np.zeros([nIter_MH2D, nLine, nEnergymh2d])
        for iIter in range(nIter_MH2D):
            iLine_I = np.loadtxt(files_mh2dfile[iIter], skiprows=nLineHeader, \
                unpack=True, usecols=[0], dtype=float).astype(np.int64)
            LonLat_II = np.loadtxt(files_mh2dfile[iIter], skiprows=nLineHeader, \
                unpack=True, usecols=[iLon, iLat] + [iLat+1+i for i in range(nEnergymh2d)], dtype=float)
            LonLat_III[iIter, iLine_I-1, 0] = LonLat_II[0]
            LonLat_III[iIter, iLine_I-1, 1] = LonLat_II[1]
            IFluxOp_III[iIter, iLine_I-1, :] = LonLat_II[2:].T # NaN and negative values are included
        dic['Position'] = LonLat_III # deg
        dic['IFluxOp'] = IFluxOp_III # consistent with info in MH_data_EChannel.H read by read_mh_echannel
    
    return dic

# Get the triangulation skeleton
def get_tri_skeleton(datax_I, datay_I):
    ### Parameters:
    ### datax_I, datay_I: 1d array for x and y coordinates
    ### ############
    ### Returns:
    ### tri_skeleton: triangulated skeleton by (datax_I, datay_I)
    ############################################################
    tri_skeleton = tri.Triangulation(datax_I, datay_I)
    return tri_skeleton

# Interpolate the fluxes using the triangulation skeleton
def tri_interp_value(tri_skeleton, dataz_I, dataxinterp, \
        datayinterp, ValueScale='linear', SmoothSigma=1.0):
    ### Parameters:
    ### tri_skeleton: triangulated skeleton
    ### dataz_I: 1d array for the value at (datax_I, datay_I)
    ### mskz_I: msk of the valid values
    ### dataxinterp, datayinterp: interpolated at (x_I, y_I)
    ### ValueScale: interpolation scale, either linear or log
    ### SmoothSigma: if need any smoothing; otherwise use 0.0
    ### ############
    ### Returns:
    ### datazinterp: interpolated values by tri_skeleton
    ############################################################
        
    if ValueScale.lower() in ["linear", "lin"]:
        tri_interpz = tri.LinearTriInterpolator(tri_skeleton, dataz_I)
        datazinterp = tri_interpz(dataxinterp, datayinterp)
    elif ValueScale.lower() in ["log", "logarithm", "logarithmic", "log10", "loge", "ln"]:
        tri_interpz = tri.LinearTriInterpolator(tri_skeleton, np.log(dataz_I))
        datazinterp = np.exp(tri_interpz(dataxinterp, datayinterp))
    else:
        print("Unrecognized ValueScale (%s) for interpolation, please re-enter!"%(ValueScale))
        return

    datazinterp = gaussian_filter(datazinterp, sigma=SmoothSigma)
    return datazinterp


# To run the script
if __name__ == "__main__":
    
    # Test 1: Read the MH_EChannel.H file
    df_mh_echannel = read_mh_echannel(filepath="./Simu/SP/run_1304LatHiPb01/MH_data_EChannel.H")
    print("df_mh_echannel:", df_mh_echannel)
    
    # Test 2: Read the mh2d files
    dic_mh2d = read_simu_mh2d(pathMH2D="./Simu/SP/run_1304LatHiPb01/")
    print("dic_mh2d keywords", dic_mh2d.keys())
