# CALLING XFOIL FROM PYTHON
# Written by: JoshTheEngineer  Modified by: Ben Mertz
# YouTube: www.youtube.com/joshtheengineer
# Website: www.joshtheengineer.com
# Started: 01/01/19
# Updated: 01/01/19 - Started code in MATLAB
#                   - Works as expected
#          02/03/19 - Transferred code from MATLAB to Python
#                   - Works as expected
#          06/26/19 - Modified to use flap deflections and save polars in AeroDyn format

import os
import numpy as np
import matplotlib.pyplot as plt
from Polar import Polar
import weio


# %% CREATE LOADING FILE

# Knowns
AFName       =  "FAAW3211" #"Profile12"  #"NACA63618" "FAAW3211"
AoA_min        = "-9"
AoA_max   = "25"
AoA_inc   = "0.5"
Re        = "9800000"
numNodes   = "250" #number of panels to use
dist_param = "0.6" #TE/LE panel density ratio
IterLimit = "200" #Maximum number of iterations to try and get to convergence
CoordsFlnmAF =  "FAAW3211coords.dat"  #"NACA63618coords.dat" "FAAW3211coords.dat" "Profile12coords.dat"
#saveFlnmAF = "NACA63618_Airfoil.txt"
#saveFlnmPolar = "NACA63618_Polar.txt"
xc_hinge = "0.8" #x/c for flap hinge location
yt_hinge = "0.5" #y/t for flap hinge location
delta_flap = "10" #Flap delfection angle in degrees
c = 5.21 #local chord length [m]
r = 21.3 #local radial location [m]
R = 100.0 #rotor radius [m]
tsr = 10.5 #tip speed ratio
neg_delta_flap = "-" + delta_flap
xfoilFlnm  = 'xfoil_input.txt'

AD_Path = "./AD_Polars/"
Profile_Path = "./Airfoils/Profiles/"
Polars_Path = "./Airfoils/Polars/"

saveFlnmAF = Profile_Path + AFName + "_" + delta_flap + "_Airfoil.txt"
saveFlnmPolar = Polars_Path + AFName + "_" + delta_flap + "_Polar.txt"
saveFlnmAFneg = Profile_Path + AFName + "_" + neg_delta_flap + "_Airfoil.txt"
saveFlnmPolarneg = Polars_Path + AFName + "_" + neg_delta_flap + "_Polar.txt"
saveFlnmAF0 = Profile_Path + AFName + "_0_Airfoil.txt"
saveFlnmPolar0 = Polars_Path + AFName + "_0_Polar.txt"

LoadFlnmAF = Profile_Path + CoordsFlnmAF
# Delete files if they exist
if os.path.exists(saveFlnmAF):
    os.remove(saveFlnmAF)

if os.path.exists(saveFlnmPolar):
    os.remove(saveFlnmPolar)

if os.path.exists(saveFlnmAFneg):
    os.remove(saveFlnmAFneg)

if os.path.exists(saveFlnmPolarneg):
    os.remove(saveFlnmPolarneg)

if os.path.exists(saveFlnmAF0):
    os.remove(saveFlnmAF0)

if os.path.exists(saveFlnmPolar0):
    os.remove(saveFlnmPolar0)
     
# Create the airfoil with flap (pos flap deflection)
fid = open(xfoilFlnm,"w")
fid.write("PLOP \n G \n\n")
fid.write("LOAD \n") 
fid.write( LoadFlnmAF + "\n")
fid.write( AFName + "\n")
fid.write("GDES \n")
fid.write("FLAP \n")
fid.write(xc_hinge + "\n")
fid.write("999\n") #To specify y/t instead of actual distance
fid.write(yt_hinge + "\n")
fid.write(delta_flap + "\n")
fid.write("NAME \n")
fid.write(AFName + "_" + delta_flap + "\n")
fid.write("EXEC \n \n")

# Re-panel with specified number of panes and LE/TE panel density ratio
fid.write("PPAR\n")
fid.write("N \n" )
fid.write(numNodes + "\n")
fid.write("T \n")
fid.write( dist_param + "\n")
fid.write("\n\n")

# Save airfoil coordinates with designation flap deflection
fid.write("PSAV \n")
fid.write( saveFlnmAF + " \n \n") #the extra \n may not be needed

# Set Simulation parameters (Re and max number of iterations)
fid.write("OPER\n")
fid.write("VISC \n")
fid.write( Re + "\n")
fid.write("ITER \n")
fid.write( IterLimit + "\n")

# Run simulations for range of AoA
fid.write("ASEQ 0 " + AoA_min + " -0.5 \n") # The preliminary runs are just to get an initialize airfoil solution at min AoA so that the actual runs will not become unstable
fid.write("PACC\n") #Toggle saving polar on 
fid.write( saveFlnmPolar + " \n \n")
fid.write("ASEQ " + AoA_min + " " + AoA_max + " " + AoA_inc + "\n") #run simulations for desired range of AoA
fid.write("PACC\n\n") #Toggle saving polar off

fid.write("QUIT \n")
fid.close()

# Run the XFoil calling command
os.system("xfoil.exe < xfoil_input.txt")

if os.path.exists(xfoilFlnm):
    os.remove(xfoilFlnm)

# Create the airfoil with flap (neg flap deflection)
fid = open(xfoilFlnm,"w")
fid.write("PLOP \n G \n\n")
fid.write("LOAD \n") 
fid.write( LoadFlnmAF + "\n")
fid.write( AFName + "\n")
fid.write("GDES \n")
fid.write("FLAP \n")
fid.write(xc_hinge + "\n")
fid.write("999\n") #To specify y/t instead of actual distance
fid.write(yt_hinge + "\n")
fid.write(neg_delta_flap + "\n")
fid.write("NAME \n")
fid.write(AFName + "_" + neg_delta_flap + "\n")
fid.write("EXEC \n \n")

# Re-panel with specified number of panes and LE/TE panel density ratio
fid.write("PPAR\n")
fid.write("N \n" )
fid.write(numNodes + "\n")
fid.write("T \n")
fid.write( dist_param + "\n")
fid.write("\n\n")

# Save airfoil coordinates with designation flap deflection
fid.write("PSAV \n")
fid.write( saveFlnmAFneg + " \n \n") #the extra \n may not be needed

# Set Simulation parameters (Re and max number of iterations)
fid.write("OPER\n")
fid.write("VISC \n")
fid.write( Re + "\n")
fid.write("ITER \n")
fid.write( IterLimit + "\n")

# Run simulations for range of AoA
fid.write("ASEQ 0 " + AoA_min + " -1 \n") # The preliminary runs are just to get an initialize airfoil solution at min AoA so that the actual runs will not become unstable
fid.write("PACC\n") #Toggle saving polar on 
fid.write( saveFlnmPolarneg + " \n \n")
fid.write("ASEQ " + AoA_min + " " + AoA_max + " " + AoA_inc + "\n") #run simulations for desired range of AoA
fid.write("PACC\n\n") #Toggle saving polar off

fid.write("QUIT \n")
fid.close()

# Run the XFoil calling command
os.system("xfoil.exe < xfoil_input.txt")

if os.path.exists(xfoilFlnm):
    os.remove(xfoilFlnm)

# Create the airfoil with flap (no deflection)
fid = open(xfoilFlnm,"w")
fid.write("PLOP \n G \n\n")
fid.write("LOAD \n") 
fid.write( LoadFlnmAF + "\n")
fid.write( AFName + "\n")

# Re-panel with specified number of panes and LE/TE panel density ratio
fid.write("PPAR\n")
fid.write("N \n" )
fid.write(numNodes + "\n")
fid.write("T \n")
fid.write( dist_param + "\n")
fid.write("\n\n")

# Save airfoil coordinates with designation flap deflection
fid.write("PSAV \n")
fid.write( saveFlnmAF0 + " \n \n") #the extra \n may not be needed

# Set Simulation parameters (Re and max number of iterations)
fid.write("OPER\n")
fid.write("VISC \n")
fid.write( Re + "\n")
fid.write("ITER \n")
fid.write( IterLimit + "\n")

# Run simulations for range of AoA
fid.write("ASEQ 0 " + AoA_min + " -1 \n") # The preliminary runs are just to get an initialize airfoil solution at min AoA so that the actual runs will not become unstable
fid.write("PACC\n") #Toggle saving polar on 
fid.write( saveFlnmPolar0 + " \n \n")
fid.write("ASEQ " + AoA_min + " " + AoA_max + " " + AoA_inc + "\n") #run simulations for desired range of AoA
fid.write("PACC\n\n") #Toggle saving polar off

fid.write("QUIT \n")
fid.close()

# Run the XFoil calling command
os.system("xfoil.exe < xfoil_input.txt")

if os.path.exists(xfoilFlnm):
    os.remove(xfoilFlnm)


## Convert polars from xfoil into an AeroDyn polar file This part of the code was developed by E. Branlard 
    # --- Reading a existing AD file, just as a template, we'll replace things in it
templatePolFile    = '_TemplateFiles/AeroDyn_Template.dat'
pol = weio.read(templatePolFile)

# --- Creating a Polar object from partial data
# Load the data from the text file
p = np.loadtxt(saveFlnmPolar, skiprows=12)
pneg = np.loadtxt(saveFlnmPolarneg, skiprows=12)
p0 = np.loadtxt(saveFlnmPolar0, skiprows=12)
#p=weio.read(saveFlnmPolar).toDataFrame().values

oldpolar= Polar(float(Re), p[:,0],p[:,1],p[:,2],p[:,4])
oldpolarneg= Polar(float(Re), pneg[:,0],pneg[:,1],pneg[:,2],pneg[:,4])
oldpolar0= Polar(float(Re), p0[:,0],p0[:,1],p0[:,2],p0[:,4])


# Correct for 3D Effects

rR = r/R
cr = c/r
polar3d = oldpolar.correction3D(rR,cr,tsr)
polar3dneg = oldpolarneg.correction3D(rR,cr,tsr)
polar3d0 = oldpolar0.correction3D(rR,cr,tsr)


# Extrapolate to -180 180
cdmax = 1.5
polar = polar3d.extrapolate(cdmax)
polarneg = polar3dneg.extrapolate(cdmax)
polar0 = polar3d0.extrapolate(cdmax) 
#polar = oldpolar.extrapolate(cdmax)
#polarneg = oldpolarneg.extrapolate(cdmax)
#polar0 = oldpolar0.extrapolate(cdmax) 

#plt.plot(oldpolar.alpha,oldpolar.cl,'+',label = 'cl_old'   )
#plt.plot(polar.alpha   ,polar.cl   ,label       = 'cl_')
#plt.legend()
#plt.show()

(alpha0,alpha1,alpha2,cnSlope,cn1,cn2,cd0,cm0)=polar.unsteadyParams()
(alpha0n,alpha1n,alpha2n,cnSlopen,cn1n,cn2n,cd0n,cm0n)=polarneg.unsteadyParams()
(alpha00,alpha10,alpha20,cnSlope0,cn10,cn20,cd00,cm00)=polar0.unsteadyParams()

PolarTable = np.column_stack((polar.alpha,polar.cl,polar.cd,polar.cm))
PolarTableneg = np.column_stack((polarneg.alpha,polarneg.cl,polarneg.cd,polarneg.cm))
PolarTable0 = np.column_stack((polar0.alpha,polar0.cl,polar0.cd,polar0.cm))

#print(cnSlope,cn1,cn2)
# Setting unsteady parameters
pol['Re'] = float(Re)/(1e6) # TODO UNKNOWN
if np.isnan(alpha0):
    pol['alpha0'] = 0
else:
    pol['alpha0'] = np.around(alpha0, 4)
pol['alpha1']    = np.around(alpha1, 4) # TODO approximate
pol['alpha2']    = np.around(alpha2, 4) # TODO approximate
pol['C_nalpha']  = np.around(cnSlope ,4)
pol['Cn1']       = np.around(cn1, 4)    # TODO verify
pol['Cn2']       = np.around(cn2, 4)
pol['Cd0']       = np.around(cd0, 4)
pol['Cm0']       = np.around(cm0, 4)

# Setting
pol['NumAlf'] = polar.cl.shape[0]
pol['AFCoeff'] = np.around(PolarTable, 5)
filename= AD_Path + AFName + "_" + delta_flap + "_AeroDyn15_Polar.dat" #'Polar_out.dat'
pol.write(filename)

#Do it again for neg flap deflection
if np.isnan(alpha0n):
    pol['alpha0'] = 0
else:
    pol['alpha0'] = np.around(alpha0n, 4)
pol['alpha1']    = np.around(alpha1n, 4) # TODO approximate
pol['alpha2']    = np.around(alpha2n, 4) # TODO approximate
pol['C_nalpha']  = np.around(cnSlopen ,4)
pol['Cn1']       = np.around(cn1n, 4)    # TODO verify
pol['Cn2']       = np.around(cn2n, 4)
pol['Cd0']       = np.around(cd0n, 4)
pol['Cm0']       = np.around(cm0n, 4)

# Setting
pol['NumAlf'] = polarneg.cl.shape[0]
pol['AFCoeff'] = np.around(PolarTableneg, 5)
filename= AD_Path + AFName + "_n" + delta_flap + "_AeroDyn15_Polar.dat" #'Polar_out.dat'
pol.write(filename)

# Do it again for no flap
if np.isnan(alpha00):
    pol['alpha0'] = 0
else:
    pol['alpha0'] = np.around(alpha00, 4)
pol['alpha1']    = np.around(alpha10, 4) # TODO approximate
pol['alpha2']    = np.around(alpha20, 4) # TODO approximate
pol['C_nalpha']  = np.around(cnSlope0 ,4)
pol['Cn1']       = np.around(cn10, 4)    # TODO verify
pol['Cn2']       = np.around(cn20, 4)
pol['Cd0']       = np.around(cd00, 4)
pol['Cm0']       = np.around(cm00, 4)

# Setting
pol['NumAlf'] = polar0.cl.shape[0]
pol['AFCoeff'] = np.around(PolarTable0, 5)
filename= AD_Path + AFName + "_0_AeroDyn15_Polar.dat" #'Polar_out.dat'
pol.write(filename)




#print('Writing polar to file:',filename,' thick={}'.format(t))

