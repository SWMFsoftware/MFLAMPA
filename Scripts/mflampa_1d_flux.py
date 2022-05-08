#load all modules and set globle variables
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib
import numpy as np
import os
import glob
import sys

cols=['ID','LagrID', 'X', 'Y', 'Z', 'Rho', 'T', 'Ux', 'Uy', 'Uz', 'Bx', 'By', 'Bz',
       'Wave1', 'Wave2', 'flux_total', 'flux_00005', 'flux_00010',
       'flux_00030', 'flux_00050', 'flux_00060', 'flux_00100', 'eflux']
mu0=1.25663706E-6
rootdir='/Users/zhlulu/SWMF/SWMFSOLAR/Results/MFLAMPA/RESULTS_SEP_no_rotate/SP/'

nlon=8
nlat=4

def plot_1d_flux_1line_1time(filename,ax):
    data=pd.read_csv(filename,skiprows=16,sep='\s+',names=cols)
    data.columns= data.columns.str.lower()
    data['r']=(data['x']**2.+data['y']**2.+data['z']**2.)**0.5
    data['u']=(data['ux']**2.+data['uy']**2.+data['uz']**2.)**0.5
    data['b']=(data['bx']**2.+data['by']**2.+data['bz']**2.)**0.5
    data['db2']=(data['wave1']+data['wave2'])*mu0

    ax[0,0].plot(data['r'],data['rho']*data['r']**2.)
    ax[1,0].plot(data['r'],data['u'])
    ax[2,0].plot(data['r'],data['db2']/data['b']**2)

    ax[0,1].plot(data['r'],data['flux_00005'],label=ifile,color=colors[ifile])
    ax[1,1].plot(data['r'],data['flux_00010'],label=ifile,color=colors[ifile])
    ax[2,1].plot(data['r'],data['flux_00100'],label=ifile,color=colors[ifile])

    ax[0,0].set_ylabel('$\\rho r^2$')
    ax[1,0].set_ylabel('U [km/s]')
    ax[2,0].set_ylabel('$(db/b)^2$')
    ax[0,0].set_yscale('log')
    ax[2,0].set_yscale('log')

    ax[0,1].set_ylabel('>5 MeV')
    ax[1,1].set_ylabel('>10 MeV')
    ax[2,1].set_ylabel('>100 MeV')

    ax[0,0].set_ylim(1E10,1E14)
    ax[1,0].set_ylim(0,2E6)
    ax[2,0].set_ylim(1E-4,1)
    
    ax[0,1].set_ylim(1E-3,1E6)
    ax[1,1].set_ylim(1E-3,1E5)
    ax[2,1].set_ylim(1E-3,1E4)
    
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
        files=sorted(glob.glob(rootdir+'MH_data_'+str(ilon).zfill(3)+'_'+str(ilat).zfill(3)+'*.dat'))
        colors = [matplotlib.cm.jet(x) for x in np.linspace(0, 1, len(files))]
        
        fig,ax=plt.subplots(3,2,figsize=(6,6))
        for ifile,filename in enumerate(files):
            plot_1d_flux_1line_1time(filename,ax)
            fig.suptitle('ilon_ilat = '+str(ilon).zfill(3)+'_'+str(ilat).zfill(3))
        plt.subplots_adjust(left=0.13, bottom=0.08, right=0.98, top=0.93, wspace=0.5, hspace=0.2)
        #fig.savefig(rootdir+'fig/1dflux/mhd_sep_'+str(ilon).zfill(3)+'_'+str(ilat).zfill(3)+'.jpg',dpi=300)
        #plt.close(fig)     
        plt.show()
        exit()