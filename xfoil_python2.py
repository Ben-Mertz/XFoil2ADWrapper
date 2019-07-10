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
from xfoil2AD import xfoil2AD

# %% CREATE LOADING FILE

# Knowns

AFName       =  "Profile50" #"Profile12"  #"NACA63618" "FAAW3211"
Re        = 9800000
c = 5.21 #local chord length [m]
r = 21.3 #local radial location [m]
R = 100.0 #rotor radius [m]
tsr = 10.5 #tip speed ratio
rR=r/R
cR=c/R
cdmax=1.5
out_flnm = "./AD_Polars/" + AFName + "_Polar_mulTab.dat"

# Calling class that will be used to run xfoil and create AD polar
x2AD = xfoil2AD(Re,AFName)

flnm_n12=x2AD.runXfoil(-12.0) # Running xfoil for maximum negative flap deflection angle
flnm_0=x2AD.runXfoil(0.0) # Running xfoil for zero flap deflection angle
flnm_12=x2AD.runXfoil(12.0) # Running xfoil for maximum positive flap deflection angle

# Applying 3d corrections, extrapolating from -180 deg to 180 deg, and calculating unsteady parameters (using code based on AirfoilPrep)
(us_params_n12, PolarTab_n12) = x2AD.AFCorrections(flnm_n12,rR,cR,tsr,cdmax) 
(us_params_0, PolarTab_0) = x2AD.AFCorrections(flnm_0,rR,cR,tsr,cdmax)
(us_params_12, PolarTab_12) = x2AD.AFCorrections(flnm_12,rR,cR,tsr,cdmax)

# Saving data in an AeroDyn15 polar file with all three tables (make sure control parameter goes in assendign order)
x2AD.ADwriteFile(out_flnm,-12.0,PolarTab_n12,1,3,True,us_params_n12)
x2AD.ADwriteFile(out_flnm,0.0,PolarTab_0,2,3,True,us_params_0)
x2AD.ADwriteFile(out_flnm,12.0,PolarTab_12,3,3,True,us_params_12)

print('\n\n', out_flnm, ' successfully created')

#print('Writing polar to file:',filename,' thick={}'.format(t))

