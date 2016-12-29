#!/usr/bin/python

import numpy as np
import matplotlib.pyplot as plt
#import matplotlib as mpl
#import sys
#from matplotlib.font_manager import fontManager, FontProperties
#from pylab import plot, show, title, xlabel, ylabel, subplot, xlim
from scipy import fft, arange, signal
import cStringIO   # To be able to rm ( and ) in the data file
from scipy.interpolate import griddata    # To plot surfaces

#mpl.rcParams['font.family']='serif'
#mpl.rcParams['font.size']=17
#mpl.rcParams['legend.fontsize']=16
#mpl.rcParams['xtick.labelsize']='13'
#mpl.rcParams['ytick.labelsize']='13'

################################################################
"""
TO DO NEXT:


-figure out why the spectrum of the filtered field does not match exactly the unfilter one for the big scales
-figure out why the filter field has higher values of velocity
    -could it be that I need fft of U an not of u_line?


"""
################################################################
# READ DATA FILE
################################################################

with open ("U", "r") as myfile:
    data=myfile.read().replace('(', '')  
    data=data.replace(')','')  
    #data=myfile.read().replace('(', '').replace(')','')  #Funky way to replace two things in same line, instead of the two lines above!
    num_lines = sum(1 for line in myfile) 

field_data=np.genfromtxt(cStringIO.StringIO(data), skip_header=22, skip_footer=29)


u=field_data[:, 0]
v=field_data[:, 1]
w=field_data[:, 2]



n=128 # number of cells in one direction


u_line=np.reshape(u, (n, n, n)) # To take u in a line of cells to do the fft
v_line=np.reshape(v, (n, n, n))
w_line=np.reshape(w, (n, n, n))

#print u_line
#print np.shape(u_line)



################################################################
# CALCULATE 1D SPECTRA AND PLOT STREAMWISE COMPONENT
################################################################
#l=2*np.pi                      #The lenght of the cube
knum1d = np.arange(n)
#Fs=l/n                         # sample frequency - delta distance
#frq = 2.0*np.pi*(k/(n*Fs))
uHat_line=np.zeros((n, n, n), dtype=complex);
spectrum=np.zeros(n);

#To move to spectral space and average 
for k in range (0, n):
    for j in range (0, n):
        uHat_line[k, j, :]= np.fft.fft(u_line[k,j,:])
        #print uHat_line[k, j, :]
        spectrum=spectrum+abs((uHat_line[k, j, :])**2.)   ### NOT NORMALIZED!!
        #plt.loglog(frq, abs((uHat_line[k, j, :])**2.), marker='o')

plt.loglog(knum1d[range(n/2)], spectrum[range(n/2)]/(n*n), marker='o', label='Not filtered - 128^3 cells - 1D fft') ### NOT NORMALIZED!!

"""
##To come back to physical space
uinv_line=np.zeros((n, n, n), dtype=complex);
for k in range (0, n):
    for j in range (0, n):
        uinv_line[k, j, :]= np.fft.ifft(uHat_line[k,j,:])

#print uinv_line.real
"""

################################################################
# CALCULATE SPECTRA WITH 3D FFT AND PLOT STREAMWISE COMPONENT
################################################################

cellNum = np.linspace(0, n-1, n)

kx = np.linspace(0, n-1, n)
ky = np.linspace(0, n-1, n)
kz = np.linspace(0, n-1, n)

kx0=0 #n/2+1; ##??
ky0=0 #n/2+1;
kz0=0 #n/2+1;

knum = np.sqrt((kx-kx0)**2.+(ky-ky0)**2.+(kz-kz0)**2.) 

uHat_line=np.zeros((n, n, n), dtype=complex);
vHat_line=np.zeros((n, n, n), dtype=complex);
wHat_line=np.zeros((n, n, n), dtype=complex);
spectrum=np.zeros(n);

uHat_line=np.fft.fftn(u_line)
vHat_line=np.fft.fftn(v_line)
wHat_line=np.fft.fftn(w_line)

#Average
for k in range (0, n):
    for j in range (0, n):
        spectrum=spectrum+abs((uHat_line[k, j, :])**2.)  ### NOT NORMALIZED!!
        #plt.loglog(knum, abs((uHat_line[k, j, :])**2.), marker='*')

plt.loglog(knum[range(n/2)], spectrum[range(n/2)]/(n*n), marker='o', linewidth=2.0, label='Not filtered - 128^3 cells - 3D fft') ### NOT NORMALIZED!!


"""
% U is a 2D slice of turbulence box
filtsize = 8; % must be adjusted
order = 8; % must be adjusted
nx=size(U,2);ny=size(U,1);
n=min(nx,ny); 
kx0=nx/2+1; 
ky0=ny/2+1;
[ky,kx] = meshgrid(1:ny,1:nx);
k = sqrt((kx-kx0).^2+(ky-ky0).^2);
T = 1./(1+(k/(n/filtsize)).^(order/2));
sp=fftshift(fft2(U));
filtU = real(ifft2(ifftshift(sp.*T)));
"""

################################################################
# FILTER TO A 64^3 MESH WITH CUT-OFF FILTER
################################################################

nFilt=64
uHatFilt_line=np.zeros((nFilt, nFilt, nFilt), dtype=complex);
spectrum=np.zeros(nFilt);

delta=2.0*np.pi/nFilt
kc=np.pi/delta  #k cut off, kmax

uHatFilt_line=uHat_line[0:nFilt, 0:nFilt, 0:nFilt]
vHatFilt_line=vHat_line[0:nFilt, 0:nFilt, 0:nFilt]
wHatFilt_line=wHat_line[0:nFilt, 0:nFilt, 0:nFilt]

#print uHat_line[nFilt-1, nFilt-1, nFilt-1].real
#print np.shape(uHat_line)

#print uHatFilt_line[-1, -1, -1].real
#print np.shape(uHatFilt_line)


#Average
for k in range (0, nFilt):
    for j in range (0, nFilt):
        spectrum=spectrum+abs((uHatFilt_line[k, j, :])**2.)  ### NOT NORMALIZED!!
        #plt.loglog(knum, abs((uHat_line[k, j, :])**2.), marker='*')

plt.loglog(knum[range(nFilt/2)], spectrum[range(nFilt/2)]/(nFilt*nFilt), marker='o', linewidth=2.0, label='Filtered - 64^3 cells - 3D fft') ### NOT NORMALIZED!!

uFilt_line = np.real(np.fft.ifftn(uHatFilt_line));
vFilt_line = np.real(np.fft.ifftn(vHatFilt_line));
wFilt_line = np.real(np.fft.ifftn(wHatFilt_line));


uFilt=np.zeros(nFilt*nFilt*nFilt);
vFilt=np.zeros(nFilt*nFilt*nFilt);
wFilt=np.zeros(nFilt*nFilt*nFilt);

uFilt=np.reshape(uFilt_line, nFilt*nFilt*nFilt);
vFilt=np.reshape(vFilt_line, nFilt*nFilt*nFilt);
wFilt=np.reshape(wFilt_line, nFilt*nFilt*nFilt);

################################################################
# SPECTRA PLOT DETAILS
################################################################

legend=plt.legend(loc=0,  frameon = 1, ncol=1)
frame = legend.get_frame()
frame.set_facecolor('0.92')


#plt.xlabel(r'$\kappa^*_{1}\,\rm{[-]}$', size='17.5')
#plt.ylabel(r'$E^*_{11} (\kappa^*_{1}) \,\,\rm{[-]}$', size='20')
#plt.xlim(1, 1e2)
#plt.ylim(1e-6, 1e-1)

plt.savefig('DITSpectra_1stTry.pdf', dpi=800)

#plt.grid(False)
plt.show()


################################################################
# PLOT A SLICE OF Ux BEFORE AND AFTER FILTERING
################################################################

#Domain size
xmin=0.
xmax=1.
ymin=0.
ymax=1.

########
# 128 cells
########

n=128
shift=(xmax-xmin)/(2.*n) # To plot the center of the first cell
z=u_line[n/2,:,:]

xi = np.linspace(xmin+shift,xmax-shift, n)
yi = np.linspace(ymin,ymax, n)

x, y = np.meshgrid(xi, yi)
levels = np.arange(z.min(), z.max(), 0.2)

fig=plt.figure()
ax = fig.add_subplot(1,1,1)
CS = plt.contourf(x, y, z, levels, extend='both')   #, cmap=plt.cm.gist_gray)

plt.xlabel(r"$x/L\,\,\rm{[-]}$",fontsize=14)
plt.ylabel(r"$y/L\,\,\rm{[-]}$", fontsize=14)

# Make a colorbar for the ContourSet returned by the contourf call.
cbar = plt.colorbar(CS)
cbar.ax.set_ylabel(r"$u\,\,\rm{[-]}$", fontsize=14)

plt.savefig('visualize_u.pdf', dpi=800)
plt.show()


########
# 64 cells    
########

n=64     
shift=(xmax-xmin)/(2.*n) # To plot the center of the first cell
z=uFilt_line[n/2,:,:]  

xi = np.linspace(xmin+shift,xmax-shift, n)
yi = np.linspace(ymin,ymax, n)

x, y = np.meshgrid(xi, yi)
#levels = np.arange(z.min(), z.max(), 0.2)

fig=plt.figure()
ax = fig.add_subplot(1,1,1)
CS = plt.contourf(x, y, z, levels, extend='both')   #, cmap=plt.cm.gist_gray)

plt.xlabel(r"$x/L\,\,\rm{[-]}$",fontsize=14)
plt.ylabel(r"$y/L\,\,\rm{[-]}$", fontsize=14)

# Make a colorbar for the ContourSet returned by the contourf call.
cbar = plt.colorbar(CS)
cbar.ax.set_ylabel(r"$u\,\,\rm{[-]}$", fontsize=14)

plt.savefig('visualize_uFiltered.pdf', dpi=800)
plt.show()

############################################################################## 
#
#      To print filtered field in OF format
#
##############################################################################

with open ("UFiltered", "w") as outFile:
    outFile.write('\*--------------------------------*- C++ -*----------------------------------*/'+'\n' 
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
                  +str(nFilt*nFilt*nFilt)+'\n'
                  "("+'\n')

    for i in range (0, np.size(uFilt)): # or (0,num_lines-22-29-2):   #rm header lines, footer lines and because it starts at 0.
        outFile.write("(" +str(uFilt[i])+ "\t\t" +str(vFilt[i]) + "\t\t" +str(wFilt[i])+")"+'\n')


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


###################################





"""
QUESTIONS


1. fft2 y fftn give similar spectrum but shifted up and left for fftn ? What does that mean? 



"""


