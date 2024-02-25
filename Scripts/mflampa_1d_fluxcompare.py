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
rootdir='/Users/whliu/test3/SWMF/SP/MFLAMPA/Doc/'
# rootdir='../Doc/SP-New/IO2/MH_data_001_001_e20120123_040000_n000000.out'

nlon=4
nlat=4
ntime=1

def plot_1d_flux_1line_1time(fileList, ax):
    
    def read_outfile(filename):
        data=pd.read_csv(filename,skiprows=16,sep='\s+',names=cols)
        data.columns= data.columns.str.lower()
        data['ux'] /= 1000; data['uy'] /= 1000; data['uz'] /= 1000
        data['r']=(data['x']**2.+data['y']**2.+data['z']**2.)**0.5
        data['u']=(data['ux']**2.+data['uy']**2.+data['uz']**2.)**0.5
        data['b']=(data['bx']**2.+data['by']**2.+data['bz']**2.)**0.5
        data['db2']=(data['wave1']+data['wave2'])*mu0
        
        return data

    dataTest = [read_outfile(fileList[iFile]) for iFile in range(len(fileList))]
    dataLegend = ['Test: Non-\nConservative Advection', 'Test: Conservative \nPoisson Advection', 'Test: Poisson and GCR \nSpectrum at 1 AU']
    
    for iFile in range(3):
        print(iFile, dataTest[iFile]['flux_00005'][0:10])
        ax[0,0].plot(dataTest[iFile]['r'],dataTest[iFile]['rho']*dataTest[iFile]['r']**2., ls='-.', color=colors[-iFile-2])
        ax[1,0].plot(dataTest[iFile]['r'],dataTest[iFile]['u'], ls='-.', color=colors[-iFile-2])
        ax[2,0].plot(dataTest[iFile]['r'],dataTest[iFile]['db2']/dataTest[iFile]['b']**2, ls='-.', color=colors[-iFile-2])

        ax[0,1].plot(dataTest[iFile]['r'],dataTest[iFile]['flux_00005'],label=dataLegend[iFile],ls='-.', color=colors[-iFile-2])
        ax[1,1].plot(dataTest[iFile]['r'],dataTest[iFile]['flux_00010'],ls='-.', color=colors[-iFile-2])
        ax[2,1].plot(dataTest[iFile]['r'],dataTest[iFile]['flux_00100'],ls='-.', color=colors[-iFile-2])

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
    
    ax[0,1].set_ylim(5E-3,5E-1)
    ax[1,1].set_ylim(5E-3,5E-1)
    ax[2,1].set_ylim(5E-3,5E-1)
    
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
        
        files_mflampa = sorted(glob.glob(rootdir+'run_mflampa/IO2/MH_data_'+str(ilon).zfill(3)+'_'+str(ilat).zfill(3)+'*.out'))
        files_poisson = sorted(glob.glob(rootdir+'run_poisson/IO2/MH_data_'+str(ilon).zfill(3)+'_'+str(ilat).zfill(3)+'*.out'))
        files_poissoncrphi = sorted(glob.glob(rootdir+'run_poissoncrphi/IO2/MH_data_'+str(ilon).zfill(3)+'_'+str(ilat).zfill(3)+'*.out'))
        colors = [matplotlib.cm.rainbow(x) for x in np.linspace(0, 1, len(files_mflampa))]
        
        fig, ax = plt.subplots(3,2,figsize=(6,6))
        for ifile, filename1 in enumerate(files_mflampa):
            if ifile == ntime:
                timestr = filename1[-27:-12]
                filename2 = files_poisson[ifile]
                filename3 = files_poissoncrphi[ifile]
                plot_1d_flux_1line_1time([filename1, filename2, filename3], ax)
                fig.suptitle('ilon_ilat = '+str(ilon).zfill(3)+'_'+str(ilat).zfill(3)+' at %s'%(timestr))
        plt.subplots_adjust(left=0.13, bottom=0.08, right=0.98, top=0.93, wspace=0.6, hspace=0.2)
        fig.savefig(rootdir+'Fig/Compare/mhd_sep_'+str(ilon).zfill(3)+'_'+str(ilat).zfill(3)+'_%.2d_.pdf'%(ntime),dpi=450)
        plt.show(); plt.close(fig)
        exit()