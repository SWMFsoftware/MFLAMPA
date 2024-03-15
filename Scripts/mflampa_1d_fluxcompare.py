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
rootdir='../Doc/'
# rootdir='../Doc/SP-New/IO2/MH_data_001_001_e20120123_040000_n000000.out'

nlon=4
nlat=4
ntime=6

plt.rc('text', usetex=True); plt.rc('font', family='serif')

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
    dataLegend = ['Test: Non-\nConservative Advection', 'Test: Conservative \nPoisson Advection', 
                'Test: Poisson and GCR', 'Test: Steady State']
    lWidth = [1.9, 1.7, 1.6, 1.5]
    
    for iFile in range(4):
        print(iFile, dataTest[iFile]['flux_00005'][0:10])
        ax[0,0].plot(dataTest[iFile]['r'],dataTest[iFile]['rho']*dataTest[iFile]['r']**2., 
                    label=dataLegend[iFile], ls='-.', color=colors[-iFile-2])
        ax[1,0].plot(dataTest[iFile]['r'],dataTest[iFile]['u'], ls='-.', color=colors[-iFile-2])
        ax[2,0].plot(dataTest[iFile]['r'],dataTest[iFile]['db2']/dataTest[iFile]['b']**2, ls='-.', color=colors[-iFile-2])

        ax[0,1].plot(dataTest[iFile]['r'],dataTest[iFile]['flux_00005'],ls='-.', color=colors[-iFile-2], lw=lWidth[iFile])
        ax[1,1].plot(dataTest[iFile]['r'],dataTest[iFile]['flux_00010'],ls='-.', color=colors[-iFile-2], lw=lWidth[iFile])
        ax[2,1].plot(dataTest[iFile]['r'],dataTest[iFile]['flux_00100'],ls='-.', color=colors[-iFile-2], lw=lWidth[iFile])

    ax[0,0].set_ylabel(r'$\rho r^2$ [amu/m$^3 \cdot R_\mathrm{s}^2$]')
    ax[1,0].set_ylabel(r'$U$ [km/s]')
    ax[2,0].set_ylabel(r'$(\mathrm{d}B/B)^2$')
    ax[0,0].set_yscale('log')
    ax[2,0].set_yscale('log')

    ax[0,1].set_ylabel('$>5$ MeV [pfu]')
    ax[1,1].set_ylabel('$>10$ MeV [pfu]')
    ax[2,1].set_ylabel('$>100$ MeV [pfu]')

    ax[0,0].set_ylim(1E11,1E14)
    ax[1,0].set_ylim(0,800)
    ax[2,0].set_ylim(1E-4,1E+2)
    
    ax[0,1].set_ylim(5E-3,5E-1)
    ax[1,1].set_ylim(5E-3,5E-1)
    ax[2,1].set_ylim(5E-3,5E-1)
    
    ax[0,0].legend(prop={"size": 9})
    
    for px in ax[:,0]:
        px.set_xlim([0,250])
    for px in ax[:,1]:
        px.set_yscale('log')
        px.set_xlim([0,250])
    for px in ax[2,:]:
        px.set_xlabel(r'$r$ [$R_\mathrm{s}$]')
    return ax

for ilon in range(1,nlon+1):
    for ilat in range(1,nlat+1):
        
        files_mflampa = sorted(glob.glob(rootdir+'run_mflampa/SP/IO2/MH_data_'+str(ilon).zfill(3)+'_'+str(ilat).zfill(3)+'*.out'))
        files_poisson = sorted(glob.glob(rootdir+'run_poisson/SP/IO2/MH_data_'+str(ilon).zfill(3)+'_'+str(ilat).zfill(3)+'*.out'))
        files_poissoncrphi = sorted(glob.glob(rootdir+'run_poissoncrphi/IO2/MH_data_'+str(ilon).zfill(3)+'_'+str(ilat).zfill(3)+'*.out'))
        files_steady = sorted(glob.glob(rootdir+'run_steady/SP/IO2/MH_data_'+str(ilon).zfill(3)+'_'+str(ilat).zfill(3)+'*n000020.out'))
        colors = [matplotlib.cm.rainbow(x) for x in np.linspace(0, 1, len(files_mflampa))]

        fig, ax = plt.subplots(3,2,figsize=(6,6),sharex=True)
        for ifile, filename1 in enumerate(files_mflampa):
            if ifile == ntime:
                timestr = filename1[-27:-12]
                filename2 = files_poisson[ifile]
                filename3 = files_poissoncrphi[ifile]
                plot_1d_flux_1line_1time([filename1, filename2, filename3, files_steady[0]], ax)
                fig.suptitle('ilon_ilat = '+str(ilon).zfill(3)+'_'+str(ilat).zfill(3)+' at %s'%(timestr))
        plt.subplots_adjust(left=0.11, bottom=0.08, right=0.98, top=0.95, wspace=0.3, hspace=0.15)
        fig.savefig(rootdir+'Fig/Compare/mhd_sep_'+str(ilon).zfill(3)+'_'+str(ilat).zfill(3)+'_%.2d_.pdf'%(ntime),dpi=450)
        plt.show(); plt.close(fig)
        exit()