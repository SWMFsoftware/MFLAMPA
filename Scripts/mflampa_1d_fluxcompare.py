#load all modules and set globle variables
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib
import numpy as np
import os
import glob
import sys

cols=['LagrID', 'X', 'Y', 'Z', 'Rho', 'T', 'Ux', 'Uy', 'Uz', 'Bx', 'By', 'Bz',
       'Wave1', 'Wave2', 'flux_total', 'flux_00005', 'flux_00010',
       'flux_00030', 'flux_00050', 'flux_00060', 'flux_00100', 'eflux']
mu0=1.25663706E-6
rootdir='/Users/whliu/SWMF/SP/MFLAMPA/Doc/'
# rootdir='../Doc/SP-New/IO2/MH_data_001_001_e20120123_040000_n000000.out'

nlon=4
nlat=4
ntime=1

def plot_1d_flux_1line_1time(filename, filenameOld, ax):
    
    def read_outfile(filename):
        data=pd.read_csv(filename,skiprows=16,sep='\s+',names=cols)
        data.columns= data.columns.str.lower()
        data['ux'] /= 1000; data['uy'] /= 1000; data['uz'] /= 1000
        data['r']=(data['x']**2.+data['y']**2.+data['z']**2.)**0.5
        data['u']=(data['ux']**2.+data['uy']**2.+data['uz']**2.)**0.5
        data['b']=(data['bx']**2.+data['by']**2.+data['bz']**2.)**0.5
        data['db2']=(data['wave1']+data['wave2'])*mu0
        
        return data

    dataNew = read_outfile(filename)
    dataOld = read_outfile(filenameOld)

    ax[0,0].plot(dataNew['r'],dataNew['rho']*dataNew['r']**2.,color=colors[ifile])
    ax[1,0].plot(dataNew['r'],dataNew['u'],color=colors[ifile])
    ax[2,0].plot(dataNew['r'],dataNew['db2']/dataNew['b']**2,color=colors[ifile])

    ax[0,1].plot(dataNew['r'],dataNew['flux_00005'],label="Advect Via \nPossion Bracket Alg.",color=colors[ifile])
    ax[1,1].plot(dataNew['r'],dataNew['flux_00010'],color=colors[ifile])
    ax[2,1].plot(dataNew['r'],dataNew['flux_00100'],color=colors[ifile])
    
    ax[0,0].plot(dataOld['r'],dataOld['rho']*dataOld['r']**2., ls='-.', color=colors[-ifile-1])
    ax[1,0].plot(dataOld['r'],dataOld['u'], ls='-.', color=colors[-ifile-1])
    ax[2,0].plot(dataOld['r'],dataOld['db2']/dataOld['b']**2, ls='-.', color=colors[-ifile-1])

    ax[0,1].plot(dataOld['r'],dataOld['flux_00005'],label="Advect Via \nNon-Conservative Alg.",ls='-.', color=colors[-ifile-1])
    ax[1,1].plot(dataOld['r'],dataOld['flux_00010'],ls='-.', color=colors[-ifile-1])
    ax[2,1].plot(dataOld['r'],dataOld['flux_00100'],ls='-.', color=colors[-ifile-1])

    ax[0,0].set_ylabel('$\\rho r^2$')
    ax[1,0].set_ylabel('U [km/s]')
    ax[2,0].set_ylabel('$(db/b)^2$')
    ax[0,0].set_yscale('log')
    ax[2,0].set_yscale('log')

    ax[0,1].set_ylabel('>5 MeV')
    ax[1,1].set_ylabel('>10 MeV')
    ax[2,1].set_ylabel('>100 MeV')

    ax[0,0].set_ylim(1E11,1E14)
    ax[1,0].set_ylim(0,800)
    ax[2,0].set_ylim(1E-4,1E+2)
    
    ax[0,1].set_ylim(9E-3,1.2E-2)
    ax[1,1].set_ylim(9E-3,1.2E-2)
    ax[2,1].set_ylim(9E-3,1.2E-2)
    
    ax[0,1].legend(prop={"size": 9})
    
    for px in ax[:,0]:
        px.set_xlim([0,250])
    for px in ax[:,1]:
        px.set_yscale('log')
        px.set_xlim([0,250])
    for px in ax[2,:]:
        px.set_xlabel('r (Rs)')
    return ax

for ilon in range(1,nlon+1):
    for ilat in range(1,nlat+1):
        
        files_new = sorted(glob.glob(rootdir+'SP-New/IO2/MH_data_'+str(ilon).zfill(3)+'_'+str(ilat).zfill(3)+'*.out'))
        files_old = sorted(glob.glob(rootdir+'SP-Old/IO2/MH_data_'+str(ilon).zfill(3)+'_'+str(ilat).zfill(3)+'*.out'))
        colors = [matplotlib.cm.rainbow(x) for x in np.linspace(0, 1, len(files_new))]
        
        fig, ax = plt.subplots(3,2,figsize=(6,6))
        for ifile, filename in enumerate(files_new):
            if ifile == ntime:
                filenameOld = files_old[ifile]
                timestr = filename[-27:-12]
                plot_1d_flux_1line_1time(filename, filenameOld, ax)
                fig.suptitle('ilon_ilat = '+str(ilon).zfill(3)+'_'+str(ilat).zfill(3)+' at %s'%(timestr))
        plt.subplots_adjust(left=0.13, bottom=0.08, right=0.98, top=0.93, wspace=0.6, hspace=0.2)
        fig.savefig(rootdir+'Fig/Compare/mhd_sep_'+str(ilon).zfill(3)+'_'+str(ilat).zfill(3)+'_%.2d.pdf'%(ntime),dpi=450)
        plt.close(fig) 
        plt.show()
        exit()