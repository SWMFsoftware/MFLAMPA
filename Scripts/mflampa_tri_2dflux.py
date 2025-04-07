#!/usr/bin/env python3
'''
Load and Plot M-FLAMPA 2d spherical triangulated particle flux
Originally Written by Weihao Liu (whliu@umich.edu) in Apr, 2025
'''

import re, os, argparse
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
            return lines[linenum-1].strip()
    except FileNotFoundError:
        return f"File {filename} not found."
    except IndexError:
        return f"File {filename} does not have {filename} lines."
    except Exception as e:
        return f"An error occurred: {e}"

# Read the command line
def read_args():
    # argparse
    parser = argparse.ArgumentParser( \
        description="Read information for mflampa_tri2dflux")
    parser.add_argument('--dirMH2D', type=str, default='./SP/IO2/', \
        help='Path to M-FLAMPA mh2d output (default: ./SP/IO2/)')
    parser.add_argument('--pathEChannel', type=str, default='./SP/IO2/', \
        help='Path to MH_data_EChannel.H file (default: ./SP/IO2/)')
    parser.add_argument('--nLon', type=int, default=360, \
        help='nLon in the interpolated LonLat map (default: 360)')
    parser.add_argument('--nLat', type=int, default=180, \
        help='nLat in the interpolated LonLat map (default: 180)')
    parser.add_argument('--sepLon', type=float, default=180.0, \
        help='sepLon in the interpolated LonLat map (default: 180.0)')
    parser.add_argument('--ValueScale', type=str, default='log', \
        help='Value scale to use (default: log)')
    parser.add_argument('--SmoothSigma', type=float, default=1.0, \
        help='Smoothing sigma (default: 1.0)')
    
    # return as args
    args = parser.parse_args()
    return args

# Read the MH_data_EChannel.H file
def read_mh_echannel(pathEChannel="./SP/IO2/MH_data_EChannel.H"):
    # Given pathEChannel, read all the mh2d files created by M-FLAMPA
    
    with open(pathEChannel, 'r') as file:
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
def read_simu_mh2d(dirMH2D="./SP/"):
    # Given dirMH2D, read all the mh2d files created by M-FLAMPA
        
    # mh2d files: first 5 lines are meta data saved as a header
    nLineHeader = 5
    
    # Get all files aligned with the specific mh2d filename format in the given path
    fmt_mh2dfile = re.compile(r"^MH_data_R=\d{4}.\d{2}_e\d{8}_\d{6}_n\d{6}\.out$")
    files_mh2dfile = [ ]
    files = sorted(os.listdir(dirMH2D))
    for j in range(len(files)):
        if fmt_mh2dfile.match(files[j]):
            files_mh2dfile.append(dirMH2D+files[j])
    dic = {}
    
    # Time array
    nIter_MH2D = len(files_mh2dfile)
    time_I = [dt.datetime( \
        year=int(files_mh2dfile[iFile][-27:-23]), \
        month=int(files_mh2dfile[iFile][-23:-21]), \
        day=int(files_mh2dfile[iFile][-21:-19]), \
        hour=int(files_mh2dfile[iFile][-18:-16]), \
        minute=int(files_mh2dfile[iFile][-16:-14]), \
        second=int(files_mh2dfile[iFile][-14:-12])) for iFile in range(nIter_MH2D)]
    dic['Time'] = np.array(time_I)
    
    # Position and flux
    # if there is any mh2d file
    if nIter_MH2D > 0:
        # Get: Numbers of energy channels and field lines
        header_mh2d_I = read_line(files_mh2dfile[0], linenum=nLineHeader).split(" ")
        nEnergymh2d = 0
        for i in range(len(header_mh2d_I)-1, -1, -1):
            if header_mh2d_I[i] == "": del header_mh2d_I[i]
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
                unpack=True, usecols=[iLon, iLat]+[iLat+1+i for i in range(nEnergymh2d)], dtype=float)
            LonLat_III[iIter, iLine_I-1, 0] = LonLat_II[0]
            LonLat_III[iIter, iLine_I-1, 1] = LonLat_II[1]
            IFluxOp_III[iIter, iLine_I-1, :] = LonLat_II[2:].T # NaN and negative values are included
        dic['Position'] = LonLat_III # deg
        dic['IFluxOp'] = IFluxOp_III # consistent with information in MH_data_EChannel.H
    else:
        dic['Position'] = None
        dic['IFluxOp'] = None
    
    return dic

def aug_lonlatII(LonOrig_II, LatOrig_II, IFluxOrig_III, sepLon):
    ### Parameters:
    ### LonOrig_II, LatOrig_II: Original LonLat points on field lines
    ### IFluxOrig_III: Original IFlux values on field lines
    ### sepLon: Separation points for the longitude
    ###         [   0, 360] for Lon => sepLon = 180.0
    ###         [-180, 180] for Lon => sepLon = 0.0
    ### ############
    ### Returns:
    ### LonAug_II, LatAug_II: Augmented LonLat points (*2 size)
    ### IFluxAug_III: Augmented IFlux values (*2 size)
    ############################################################
    
    # Get the masks of left and right LonLat map
    MskLonL_II = LonOrig_II <= sepLon # Lon <= sepLon: add 360
    MskLonR_II = LonOrig_II >= sepLon # Lon >= sepLon: subtract 360
    # Handle the Lon array
    LonAugL_II = np.where(MskLonL_II, LonOrig_II+360.0, LonOrig_II) # Add 360 where MskLonL_II=T
    LonAugR_II = np.where(MskLonR_II, LonOrig_II-360.0, LonOrig_II) # Subtract 360 where MskLonR_II=T
    LonAug_II = np.concatenate([LonAugR_II, LonAugL_II], axis=1)
    # Handle the Lat array
    LatAug_II = np.concatenate([LatOrig_II, LatOrig_II], axis=1)
    # Handle the IFlux array
    IFluxAug_III = np.concatenate([IFluxOrig_III, IFluxOrig_III], axis=1)
    
    return LonAug_II, LatAug_II, IFluxAug_III

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
    ### dataxinterp, datayinterp: interpolation points with any dim (1D, 2D, ...)
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
    
    # Step 1: Read arg from the command line
    args = read_args()
    
    # Step 2: Read the MH_data_EChannel.H file
    pathEChannel = os.path.join(args.pathEChannel, 'MH_data_EChannel.H')
    df_mh_echannel = read_mh_echannel(pathEChannel=pathEChannel)
    print("df_mh_echannel:", df_mh_echannel)
    
    # Step 3: Read the mh2d files
    dic_mh2d = read_simu_mh2d(dirMH2D=args.dirMH2D)
    print("dic_mh2d keywords", dic_mh2d.keys())
    
    # Step 4: Check and Prepare the data for triangulation
    if (dic_mh2d['Position'] is None) or (dic_mh2d['IFluxOp'] is None):
        print("Error and Exit: No Useful mh2d Position and Flux Data!")
        exit(0)

    # 4.1: Get original (non-augmented) data
    LonOrig_II = dic_mh2d['Position'][:,:,0].copy() # nIter, nLine
    LatOrig_II = dic_mh2d['Position'][:,:,1].copy() # nIter, nLine
    IFluxOp_III = dic_mh2d['IFluxOp'].copy() # nIter, nLine, nEChannel
    # 4.2: Augmentation of plus/minus 360 degrees
    LonAug_II, LatAug_II, IFluxAug_III = aug_lonlatII( \
        LonOrig_II, LatOrig_II, IFluxOp_III, args.sepLon)
    MskIFluxAug_II = np.all(IFluxAug_III>0.0, axis=2) # nIter, nLine*2

    # 4.3: Preparations before entering the loop for iterations
    LonInterp_I = np.linspace(0.0, 360.0, args.nLon+1)
    LatInterp_I = np.linspace(0.0, 360.0, args.nLat+1)
    LATInterp_II, LONInterp_II = np.meshgrid(LatInterp_I, LonInterp_I)
    IFluxInterp_IV = np.zeros([IFluxAug_III.shape[0], IFluxAug_III.shape[-1], \
        args.nLon+1, args.nLat+1]) # nIter, nEChannel, nLon (remapped), nLat (remapped)
    print(IFluxInterp_IV.shape)
    # In each step, we:
    for iIter in range(len(dic_mh2d['Time'])):
        
        # 4.4: Get the effective and augmented data
        LonEff_I = LonAug_II[iIter][MskIFluxAug_II[iIter]] # deg
        LatEff_I = LatAug_II[iIter][MskIFluxAug_II[iIter]] # deg
        IFluxEff_II = IFluxAug_III[iIter][MskIFluxAug_II[iIter]] # Flux
    
        # Step 5: Get the triangulation skeleton
        tri_mh2d_skeleton = get_tri_skeleton(datax_I=LonEff_I, datay_I=LatEff_I)
        
        # Step 6: Interpolate based on the mh2d triangulation skeleton
        for jChannel in range(IFluxInterp_IV.shape[1]):
            IFluxInterp_IV[iIter, jChannel] = tri_interp_value( \
                tri_skeleton=tri_mh2d_skeleton, dataz_I=IFluxEff_II[:, jChannel], \
                dataxinterp=LONInterp_II, datayinterp=LATInterp_II, \
                ValueScale=args.ValueScale, SmoothSigma=args.SmoothSigma)
        print("iIter = %d"%(iIter), "Time =", dic_mh2d['Time'][iIter], \
            "with Orig LonLat Size =", LonEff_I.shape, LatEff_I.shape, \
            "to the Remapped LonLat Size =", LONInterp_II.shape, LATInterp_II.shape)
    