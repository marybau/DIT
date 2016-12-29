#!/usr/bin/python

import numpy as np
import matplotlib.pyplot as plt
#import matplotlib as mpl
#import sys
#from matplotlib.font_manager import fontManager, FontProperties
#from pylab import plot, show, title, xlabel, ylabel, subplot, xlim
from scipy import fft, arange, signal
import cStringIO   # To be able to rm ( and ) in the data file

#mpl.rcParams['font.family']='serif'
#mpl.rcParams['font.size']=17
#mpl.rcParams['legend.fontsize']=16
#mpl.rcParams['xtick.labelsize']='13'
#mpl.rcParams['ytick.labelsize']='13'


############################################################################## 
#
#      To get the data point
#
##############################################################################


with open ("U", "r") as myfile:
    data=myfile.read().replace('(', '')  
    data=data.replace(')','')  
    #data=myfile.read().replace('(', '').replace(')','')  #Funky way to replace two things in same line, instead of the two lines above!
    #num_lines = sum(1 for line in myfile) 

field_data=np.genfromtxt(cStringIO.StringIO(data), skip_header=22, skip_footer=29)


u=field_data[:, 0]
v=field_data[:, 1]
w=field_data[:, 2]

#print u

n=128 # number of cells in one direction


u_line=np.reshape(u, (n, n, n)) # To take u in a line of cells to do the fft
v_line=np.reshape(v, (n, n, n)) # NOTE for indices in OF order: u_line(z, y, x)
w_line=np.reshape(w, (n, n, n))

#print u_line
#print np.shape(u_line)

############################################################################## 
#
# Do an average of the velocity in real space in 8 consecutive cells to create a coarser mesh
#
############################################################################## 

uf=np.zeros((n*n*n)/8);   #index organized based on OF mesh ordering, so sweep in x, then y, then z. No need to do an array  like u_line
vf=np.zeros((n*n*n)/8); 
wf=np.zeros((n*n*n)/8); 
cellIndex=0
for k in range (0, n-1, 2):
    for j in range (0, n-1, 2):
        for i in range (0, n-1, 2):
            #print  u_line[k,j,i], u_line[k,j,i+1], u_line[k,j+1,i], u_line[k,j+1,i+1], u_line[k+1,j,i], u_line[k+1,j,i+1],u_line[k+1,j+1,i],u_line[k+1,j+1,i+1]
            uf[cellIndex]=(u_line[k,j,i]+u_line[k,j,i+1]+u_line[k,j+1,i]+u_line[k,j+1,i+1]+u_line[k+1,j,i]+u_line[k+1,j,i+1]+u_line[k+1,j+1,i]+u_line[k+1,j+1,i+1])/8.
            vf[cellIndex]=(v_line[k,j,i]+v_line[k,j,i+1]+v_line[k,j+1,i]+v_line[k,j+1,i+1]+v_line[k+1,j,i]+v_line[k+1,j,i+1]+v_line[k+1,j+1,i]+v_line[k+1,j+1,i+1])/8.
            wf[cellIndex]=(w_line[k,j,i]+w_line[k,j,i+1]+w_line[k,j+1,i]+w_line[k,j+1,i+1]+w_line[k+1,j,i]+w_line[k+1,j,i+1]+w_line[k+1,j+1,i]+w_line[k+1,j+1,i+1])/8.
            #print uf[cellIndex]
            #print i, j, k,  cellIndex
            cellIndex=cellIndex+1


#print uf
#print vf
#print wf


############################################################################## 
#
#      To print filtered field in OF format
#
##############################################################################

with open ("UFiltered", "w") as outFile:
    outFile.write('/*--------------------------------*- C++ -*----------------------------------*'+'\n' 
                  '| =========                 |                                                 |'+'\n'
                  '| \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox           |'+'\n'
                  '|  \\    /   O peration     | Version:  2.1.0                                 |'+'\n'
                  '|   \\  /    A nd           | Web:      www.OpenFOAM.org                      |'+'\n'
                  '|    \\/     M anipulation  |                                                 |'+'\n'
                  '\*---------------------------------------------------------------------------*/'+'\n'
                  "FoamFile"+'\n'
                  "{"+'\n'
                  "    version     2.0;"+'\n'
                  "    format      ascii;"+'\n'
                  "    class       volVectorField;"+'\n'
                  "    location    \"0\";"+'\n'
                  "    object      U;"+'\n'
                  "}"+'\n'
                  "// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //"+'\n'
                  ""+'\n'
                  "dimensions      [0 1 -1 0 0 0 0];"+'\n'
                  ""+'\n'
                  "internalField   nonuniform List<vector>"+'\n'
                  +str((n*n*n/8)) +'\n'
                  "("+'\n')

    for i in range (0, np.size(uf)): # or (0,num_lines-22-29-2):   #rm header lines, footer lines and because it starts at 0.
        outFile.write("(" +str(uf[i])+ "\t\t" +str(vf[i]) + "\t\t" +str(wf[i])+")"+'\n')


    outFile.write(");"+'\n'
                  ""+'\n'
                  "boundaryField"+'\n'
                  "{"+'\n'
                  "    inlet"+'\n'
                  "    {"+'\n'
                  "        type            cyclic;"+'\n'
                  "    }"+'\n'
                  "    outlet"+'\n'
                  "    {"+'\n'
                  "        type            cyclic;"+'\n'
                  "    }"+'\n'
                  "    top"+'\n'
                  "    {"+'\n'
                  "        type            cyclic;"+'\n'
                  "    }"+'\n'
                  "    bottom"+'\n'
                  "    {"+'\n'
                  "        type            cyclic;"+'\n'
                  "    }"+'\n'
                  "    front"+'\n'
                  "    {"+'\n'
                  "        type            cyclic;"+'\n'
                  "    }"+'\n'
                  "    back"+'\n'
                  "    {"+'\n'
                  "        type            cyclic;"+'\n'
                  "    }"+'\n'
                  "}"+'\n'
                  ""+'\n'
                  "// ************************************************************************* //")


############################################################################## 
#
# To plot and compare spectra before and after filtering
#
############################################################################## 

###Before filtering

l=2*np.pi                      #The lenght of the cube
k = np.arange(n)
Fs=l/n                         # sample frequency - delta distance
frq = 2.0*np.pi*(k/(n*Fs))
uHat_line=np.zeros((n, n, n), dtype=complex);
spectrum=np.zeros(n);


#To move to spectral space
for k in range (0, n):
    for j in range (0, n):
        uHat_line[k, j, :]= np.fft.fft(u_line[k,j,:])
        #print uHat_line[k, j, :]
        spectrum=spectrum+(abs((uHat_line[k, j, :])**2.))   ### NOT NORMALIZED!!
        #spectrum=spectrum+(1.0/(n*Fs))*(abs((uHat_line[k, j, :])**2.))*(Fs)**2   ### NOT NORMALIZED!!
        #plt.loglog(frq, abs((uHat_line[k, j, :])**2.), marker='o')

plt.loglog(frq[range(n/2)], spectrum[range(n/2)], marker='o', label='Not filtered') ### NOT NORMALIZED!!


###After filtering

n=n/2
uf_line=np.reshape(uf, (n, n, n)) # To take uf in a line of cells to do the fft

l=2*np.pi                      #The lenght of the cube
k = np.arange(n)
Fs=l/n                         # sample frequency - delta distance
frq = 2.0*np.pi*(k/(n*Fs))
uHat_line=np.zeros((n, n, n), dtype=complex);
spectrum=np.zeros(n);


#To move to spectral space
for k in range (0, n):
    for j in range (0, n):
        uHat_line[k, j, :]= np.fft.fft(uf_line[k,j,:])
        #print uHat_line[k, j, :]
        spectrum=spectrum+(abs((uHat_line[k, j, :])**2.))   ### NOT NORMALIZED!!
        #plt.loglog(frq, abs((uHat_line[k, j, :])**2.), marker='o')

plt.loglog(frq[range(n/2)], spectrum[range(n/2)], marker='*', label='Filtered') ### NOT NORMALIZED!!



legend=plt.legend(loc='lower left',  frameon = 1, ncol=2)
frame = legend.get_frame()
frame.set_facecolor('0.92')

plt.xlabel(r'$\kappa^*_{1}\,\rm{[-]}$', size='17.5')
plt.ylabel(r'$E^*_{11} (\kappa^*_{1}) \,\,\rm{[-]}$', size='20')
#plt.xlim(1, 1e2)
#plt.ylim(1e-6, 1e-1)

#plt.savefig('DITSpectra.pdf', dpi=800)

#plt.grid(False)
plt.show()

