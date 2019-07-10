import os
import numpy as np
import matplotlib.pyplot as plt
from Polar import Polar

class xfoil2AD(object):

    def __init__(self, Re, AfName):
        self.Re=str(Re)
        self.AFName=AfName
        

    def runXfoil(self, delta_flap=12.0, xc_hinge=0.8, yt_hinge=0.5, AoA_min=-9, AoA_max=25, AoA_inc=0.5):
        #This function is used to create and run xfoil simulations for a given flap deflection
        Re=self.Re
        AFName=self.AFName
        df=str(delta_flap)
        numNodes   = "250" #number of panels to use
        dist_param = "0.6" #TE/LE panel density ratio
        IterLimit = "200" #Maximum number of iterations to try and get to convergence
        CoordsFlnmAF = AFName + "coords.dat" # Note that I am assuming a filename of the form "AirfoilName"coords.dat
        Profile_Path = "./Airfoils/Profiles/" # Path to the airfoil profile coordinate file
        Polars_Path = "./Airfoils/Polars/"  # Output path where raw xfoil polars will be saved.  Note that these polars will possibly be combined together into a single AD polar.
        saveFlnmAF = Profile_Path + AFName + "_" + df + "_Airfoil.txt"
        saveFlnmPolar = Polars_Path + AFName + "_" + df + "_Polar.txt"
        LoadFlnmAF = Profile_Path + CoordsFlnmAF
        xfoilFlnm  = 'xfoil_input.txt' # Xfoil run script that will be deleted after it is no longer needed

        # Cleaning up old files to prevent replacement issues
        if os.path.exists(saveFlnmAF):
            os.remove(saveFlnmAF)
            
        if os.path.exists(saveFlnmPolar):
            os.remove(saveFlnmPolar)

        # %% Writes the Xfoil run script to read in coordinates, create flap, re-pannel, and create polar
        # Create the airfoil with flap 
        fid = open(xfoilFlnm,"w")
        fid.write("PLOP \n G \n\n") # turn off graphics
        fid.write("LOAD \n") 
        fid.write( LoadFlnmAF + "\n")
        fid.write( AFName + "\n")
        fid.write("GDES \n")
        fid.write("FLAP \n")
        fid.write(str(xc_hinge) + "\n")
        fid.write("999\n") #To specify y/t instead of actual distance
        fid.write(str(yt_hinge) + "\n")
        fid.write(df + "\n")
        fid.write("NAME \n")
        fid.write(AFName + "_" + df + "\n")
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
        fid.write("ASEQ 0 " + str(AoA_min) + " -0.5 \n") # The preliminary runs are just to get an initialize airfoil solution at min AoA so that the actual runs will not become unstable
        fid.write("PACC\n") #Toggle saving polar on 
        fid.write( saveFlnmPolar + " \n \n")
        fid.write("ASEQ " + str(AoA_min) + " " + str(AoA_max) + " " + str(AoA_inc) + "\n") #run simulations for desired range of AoA
        fid.write("PACC\n\n") #Toggle saving polar off
        
        fid.write("QUIT \n")
        fid.close()
        
        # Run the XFoil calling command
        os.system("xfoil.exe < xfoil_input.txt")
        
        # Delete Xfoil run script file
        if os.path.exists(xfoilFlnm):
            os.remove(xfoilFlnm)

        # Return filename with path that can be read in by AFCorrections()
        return saveFlnmPolar



    def AFCorrections(self, polarFlnm, rR, cR, tsr, cdmax ):
        # This code calls Polar.py which is a modified version of AirfoilPrep.py
        p = np.loadtxt(polarFlnm, skiprows=12) # Load in Xfoil polars (note, we are assuming raw Xfoil polars when skipping the first 12 lines)
        oldpolar= Polar(float(self.Re), p[:,0],p[:,1],p[:,2],p[:,4])
        polar3d = oldpolar.correction3D(rR,cR,tsr) # Apply 3D corrections (made sure to change the r/R, c/R, and tsr values appropriately when calling AFcorrections())
        polar = polar3d.extrapolate(cdmax) # Extrapolate polars for alpha between -180 deg and 180 deg
        (alpha0,alpha1,alpha2,cnSlope,cn1,cn2,cd0,cm0)=polar.unsteadyParams() # Calculate usteady parameters needs for AD polar files
        
        # Put outputs into arrays for easy passing between functions
        us_params_out = np.array([alpha0,alpha1,alpha2,cnSlope,cn1,cn2,cd0,cm0])
        PolarTable = np.column_stack((polar.alpha,polar.cl,polar.cd,polar.cm))
        
        # Return arrays of output data
        return (us_params_out, PolarTable)
    
    def ADwriteFile(self, ADfilename, Ctrl, PolarTable ,TabNum, NumTabs, us_flag=False, us_params=[]):
        # This function calls the data created by AFCorrections() and puts it into an AD polar formatted file with the filename (including path) of ADfilename
        # If you want to put multiple polars into the same AD polar file, call this function multiple times with the same filename and use a different TabNum (Table number) each time.  Since this function simply appends polars (for TabNum > 1...it will overwrite file if TabNum == 1), make sure to call function with desired data in the order that you want it in the AD file.
        # The first time you call this file for a given filename, make sure TabNum == 1 so that the header information is at the top of the AD file.
        if TabNum == 1:
            with open(ADfilename,'w') as f:
                f.write('! ------------ AirfoilInfo v1.00.x Input File ----------------------------------.\n')
                f.write('! Airfoil properties to be used with AeroDyn v15\n')
                f.write('! Generated from xfoil2AD which runs xfoil for different ctrl values (i.e. flap deflections) and outputs this AD file with multiple tables\n')
                f.write("! note that this file uses Marshall Buhl's new input file processing; start all comment lines with !\n")
                f.write('{!s:<22}   {:<11} {:}'.format('\"DEFAULT\"', 'InterpOrd', '! Interpolation order to use for quasi-steady table lookup {1=linear; 3=cubic spline; "default"} [default=3]\n'))
                f.write('{:<22.1f}   {:<11} {:}'.format(1.0, 'NonDimArea', '! The non-dimensional area of the airfoil (area/chord^2) (set to 1.0 if unsure or unneeded)\n'))
                f.write('{:<22d}   {:<11} {:}'.format(0, 'NumCoords', '! The number of coordinates in the airfoil shape file.  Set to zero if coordinates not included.\n'))
                f.write('{:<22d}   {:<11} {:}'.format(NumTabs, 'NumTabs', '! Number of airfoil tables in this file.  Each table must have lines for Re and Ctrl.\n'))
                f.write('! ------------------------------------------------------------------------------\n')
                f.write('! data for table 1\n')
                f.write('! ------------------------------------------------------------------------------\n')
                f.write('{:<22f}   {:<11} {:}'.format(float(self.Re)/(1.0e6), 'Re', '! Reynolds number in millions\n'))
                f.write('{:<22f}   {:<11} {:}'.format(Ctrl, 'Ctrl', '! Control setting (must be 0 for current AirfoilInfo)\n'))
                f.write('{!s:<22}   {:<11} {:}'.format(us_flag, 'InclUAdata', '! Is unsteady aerodynamics data included in this table? If TRUE, then include 30 UA coefficients below this line\n'))
                f.write('!........................................\n')
                if us_flag == True: # Populating unsteady parameter portion  
                    alpha0=us_params[0]
                    alpha1=us_params[1]
                    alpha2=us_params[2]
                    cnSlope=us_params[3]
                    cn1=us_params[4]
                    cn2=us_params[5]
                    cd0=us_params[6]
                    cm0=us_params[7]

                    f.write('{:<22f}   {:<11} {:}'.format(alpha0, 'alpha0', '! 0-lift angle of attack, depends on airfoil.\n'))
                    f.write('{:<22f}   {:<11} {:}'.format(alpha1, 'alpha1', '! Angle of attack at f=0.7, (approximately the stall angle) for AOA>alpha0. (deg)\n'))
                    f.write('{:<22f}   {:<11} {:}'.format(alpha2, 'alpha2', '! Angle of attack at f=0.7, (approximately the stall angle) for AOA<alpha0. (deg)\n'))
                    f.write('{:<22d}   {:<11} {:}'.format(1, 'eta_e', '! Recovery factor in the range [0.85 - 0.95] used only for UAMOD=1, it is set to 1 in the code when flookup=True. (-)\n'))
                    f.write('{:<22f}   {:<11} {:}'.format(cnSlope, 'C_nalpha', '! Slope of the 2D normal force coefficient curve. (1/rad)\n'))
                    f.write('{!s:<22}   {:<11} {:}'.format('default', 'T_f0', '! Initial value of the time constant associated with Df in the expression of Df and f''. [default = 3]\n'))
                    f.write('{!s:<22}   {:<11} {:}'.format('default', 'T_V0', '! Initial value of the time constant associated with the vortex lift decay process; it is used in the expression of Cvn. It depends on Re,M, and airfoil class. [default = 6]\n'))
                    f.write('{!s:<22}   {:<11} {:}'.format('default', 'T_p', '! Boundary-layer,leading edge pressure gradient time constant in the expression of Dp. It should be tuned based on airfoil experimental data. [default = 1.7]\n'))
                    f.write('{!s:<22}   {:<11} {:}'.format('default', 'T_VL', '! Initial value of the time constant associated with the vortex advection process; it represents the non-dimensional time in semi-chords, needed for a vortex to travel from LE to trailing edge (TE); it is used in the expression of Cvn. It depends on Re, M (weakly), and airfoil. [valid range = 6 - 13, default = 11]\n'))
                    f.write('{!s:<22}   {:<11} {:}'.format('default', 'b1', '! Constant in the expression of phi_alpha^c and phi_q^c.  This value is relatively insensitive for thin airfoils, but may be different for turbine airfoils. [from experimental results, defaults to 0.14]\n'))
                    f.write('{!s:<22}   {:<11} {:}'.format('default', 'b2', '! Constant in the expression of phi_alpha^c and phi_q^c.  This value is relatively insensitive for thin airfoils, but may be different for turbine airfoils. [from experimental results, defaults to 0.53]\n'))
                    f.write('{!s:<22}   {:<11} {:}'.format('default', 'b5', "! Constant in the expression of K'''_q,Cm_q^nc, and k_m,q.  [from  experimental results, defaults to 5]\n"))
                    f.write('{!s:<22}   {:<11} {:}'.format('default', 'A1', '! Constant in the expression of phi_alpha^c and phi_q^c.  This value is relatively insensitive for thin airfoils, but may be different for turbine airfoils. [from experimental results, defaults to 0.3]\n'))
                    f.write('{!s:<22}   {:<11} {:}'.format('default', 'A2', '! Constant in the expression of phi_alpha^c and phi_q^c.  This value is relatively insensitive for thin airfoils, but may be different for turbine airfoils. [from experimental results, defaults to 0.7]\n'))
                    f.write('{!s:<22}   {:<11} {:}'.format('default', 'A5', "! Constant in the expression of K'''_q,Cm_q^nc, and k_m,q. [from experimental results, defaults to 1]\n"))
                    f.write('{:<22d}   {:<11} {:}'.format(0, 'S1', '! Constant in the f curve best-fit for alpha0<=AOA<=alpha1; by definition it depends on the airfoil. [ignored if UAMod<>1]\n'))
                    f.write('{:<22d}   {:<11} {:}'.format(0, 'S2', '! Constant in the f curve best-fit for         AOA> alpha1; by definition it depends on the airfoil. [ignored if UAMod<>1]\n'))
                    f.write('{:<22d}   {:<11} {:}'.format(0, 'S3', '! Constant in the f curve best-fit for alpha2<=AOA< alpha0; by definition it depends on the airfoil. [ignored if UAMod<>1]\n'))
                    f.write('{:<22d}   {:<11} {:}'.format(0, 'S4', '! Constant in the f curve best-fit for         AOA< alpha2; by definition it depends on the airfoil. [ignored if UAMod<>1]\n'))
                    f.write('{:<22f}   {:<11} {:}'.format(cn1, 'Cn1', '! Critical value of C0n at leading edge separation. It should be extracted from airfoil data at a given Mach and Reynolds number. It can be calculated from the static value of Cn at either the break in the pitching moment or the loss of chord force at the onset of stall. It is close to the condition of maximum lift of the airfoil at low Mach numbers.\n'))
                    f.write('{:<22f}   {:<11} {:}'.format(cn2, 'Cn2', '! As Cn1 for negative AOAs.\n'))
                    f.write('{!s:<22}   {:<11} {:}'.format('default', 'St_sh', "! Strouhal's shedding frequency constant.  [default = 0.19]\n"))
                    f.write('{:<22f}   {:<11} {:}'.format(cd0, 'Cd0', '! 2D drag coefficient value at 0-lift.\n'))
                    f.write('{:<22f}   {:<11} {:}'.format(cm0, 'Cm0', '! 2D pitching moment coefficient about 1/4-chord location, at 0-lift, positive if nose up. [If the aerodynamics coefficients table does not include a column for Cm, this needs to be set to 0.0]\n'))
                    f.write('{:<22d}   {:<11} {:}'.format(0, 'k0', '! Constant in the \\hat(x)_cp curve best-fit; = (\\hat(x)_AC-0.25).  [ignored if UAMod<>1]\n'))
                    f.write('{:<22d}   {:<11} {:}'.format(0, 'k1', '! Constant in the \\hat(x)_cp curve best-fit.  [ignored if UAMod<>1]\n'))
                    f.write('{:<22d}   {:<11} {:}'.format(0, 'k2', '! Constant in the \\hat(x)_cp curve best-fit.  [ignored if UAMod<>1]\n'))
                    f.write('{:<22d}   {:<11} {:}'.format(0, 'k3', '! Constant in the \\hat(x)_cp curve best-fit.  [ignored if UAMod<>1]\n'))
                    f.write('{:<22d}   {:<11} {:}'.format(0, 'k1_hat', '! Constant in the expression of Cc due to leading edge vortex effects.  [ignored if UAMod<>1]\n'))
                    f.write('{!s:<22}   {:<11} {:}'.format('default', 'x_cp_bar', '! Constant in the expression of \\hat(x)_cp^v. [ignored if UAMod<>1, default = 0.2]\n'))
                    f.write('{!s:<22}   {:<11} {:}'.format('default', 'UACutout', '! Angle of attack above which unsteady aerodynamics are disabled (deg). [Specifying the string "Default" sets UACutout to 45 degrees]\n'))
                    f.write('{!s:<22}   {:<11} {:}'.format('default', 'filtCutOff', '! Cut-off frequency (-3 dB corner frequency) for low-pass filtering the AoA input to UA, as well as the 1st and 2nd derivatives (Hz) [default = 20]\n'))

                f.write('!........................................\n')
                f.write('! Table of aerodynamics coefficients\n')
                f.write('{:<22d}   {:<11} {:}'.format(PolarTable.shape[0], 'NumAlf', '! Number of data lines in the following table\n'))
                f.write('!    Alpha      Cl      Cd        Cm\n')
                f.write('!    (deg)      (-)     (-)       (-)\n')
                for row in PolarTable:
                    f.write(' '.join(['{: 2.14e}'.format(val) for val in row])+'\n')
                f.close()

        else: # for TabNum > 1 (note that this section does not start with the header information)
            with open(ADfilename,'a') as f:
                f.write('! ------------------------------------------------------------------------------\n')
                f.write("! data for table %i \n" % (TabNum))
                f.write('! ------------------------------------------------------------------------------\n')
                f.write('{:<22f}   {:<11} {:}'.format(float(self.Re)/(1.0e6), 'Re', '! Reynolds number in millions\n'))
                f.write('{:<22f}   {:<11} {:}'.format(Ctrl, 'Ctrl', '! Control setting (must be 0 for current AirfoilInfo)\n'))
                f.write('{!s:<22}   {:<11} {:}'.format(us_flag, 'InclUAdata', '! Is unsteady aerodynamics data included in this table? If TRUE, then include 30 UA coefficients below this line\n'))
                f.write('!........................................\n')
                if us_flag == True:   # Populating unsteady parameter portion  
                    alpha0=us_params[0]
                    alpha1=us_params[1]
                    alpha2=us_params[2]
                    cnSlope=us_params[3]
                    cn1=us_params[4]
                    cn2=us_params[5]
                    cd0=us_params[6]
                    cm0=us_params[7]

                    f.write('{:<22f}   {:<11} {:}'.format(alpha0, 'alpha0', '! 0-lift angle of attack, depends on airfoil.\n'))
                    f.write('{:<22f}   {:<11} {:}'.format(alpha1, 'alpha1', '! Angle of attack at f=0.7, (approximately the stall angle) for AOA>alpha0. (deg)\n'))
                    f.write('{:<22f}   {:<11} {:}'.format(alpha2, 'alpha2', '! Angle of attack at f=0.7, (approximately the stall angle) for AOA<alpha0. (deg)\n'))
                    f.write('{:<22d}   {:<11} {:}'.format(1, 'eta_e', '! Recovery factor in the range [0.85 - 0.95] used only for UAMOD=1, it is set to 1 in the code when flookup=True. (-)\n'))
                    f.write('{:<22f}   {:<11} {:}'.format(cnSlope, 'C_nalpha', '! Slope of the 2D normal force coefficient curve. (1/rad)\n'))
                    f.write('{!s:<22}   {:<11} {:}'.format('default', 'T_f0', '! Initial value of the time constant associated with Df in the expression of Df and f''. [default = 3]\n'))
                    f.write('{!s:<22}   {:<11} {:}'.format('default', 'T_V0', '! Initial value of the time constant associated with the vortex lift decay process; it is used in the expression of Cvn. It depends on Re,M, and airfoil class. [default = 6]\n'))
                    f.write('{!s:<22}   {:<11} {:}'.format('default', 'T_p', '! Boundary-layer,leading edge pressure gradient time constant in the expression of Dp. It should be tuned based on airfoil experimental data. [default = 1.7]\n'))
                    f.write('{!s:<22}   {:<11} {:}'.format('default', 'T_VL', '! Initial value of the time constant associated with the vortex advection process; it represents the non-dimensional time in semi-chords, needed for a vortex to travel from LE to trailing edge (TE); it is used in the expression of Cvn. It depends on Re, M (weakly), and airfoil. [valid range = 6 - 13, default = 11]\n'))
                    f.write('{!s:<22}   {:<11} {:}'.format('default', 'b1', '! Constant in the expression of phi_alpha^c and phi_q^c.  This value is relatively insensitive for thin airfoils, but may be different for turbine airfoils. [from experimental results, defaults to 0.14]\n'))
                    f.write('{!s:<22}   {:<11} {:}'.format('default', 'b2', '! Constant in the expression of phi_alpha^c and phi_q^c.  This value is relatively insensitive for thin airfoils, but may be different for turbine airfoils. [from experimental results, defaults to 0.53]\n'))
                    f.write('{!s:<22}   {:<11} {:}'.format('default', 'b5', "! Constant in the expression of K'''_q,Cm_q^nc, and k_m,q.  [from  experimental results, defaults to 5]\n"))
                    f.write('{!s:<22}   {:<11} {:}'.format('default', 'A1', '! Constant in the expression of phi_alpha^c and phi_q^c.  This value is relatively insensitive for thin airfoils, but may be different for turbine airfoils. [from experimental results, defaults to 0.3]\n'))
                    f.write('{!s:<22}   {:<11} {:}'.format('default', 'A2', '! Constant in the expression of phi_alpha^c and phi_q^c.  This value is relatively insensitive for thin airfoils, but may be different for turbine airfoils. [from experimental results, defaults to 0.7]\n'))
                    f.write('{!s:<22}   {:<11} {:}'.format('default', 'A5', "! Constant in the expression of K'''_q,Cm_q^nc, and k_m,q. [from experimental results, defaults to 1]\n"))
                    f.write('{:<22d}   {:<11} {:}'.format(0, 'S1', '! Constant in the f curve best-fit for alpha0<=AOA<=alpha1; by definition it depends on the airfoil. [ignored if UAMod<>1]\n'))
                    f.write('{:<22d}   {:<11} {:}'.format(0, 'S2', '! Constant in the f curve best-fit for         AOA> alpha1; by definition it depends on the airfoil. [ignored if UAMod<>1]\n'))
                    f.write('{:<22d}   {:<11} {:}'.format(0, 'S3', '! Constant in the f curve best-fit for alpha2<=AOA< alpha0; by definition it depends on the airfoil. [ignored if UAMod<>1]\n'))
                    f.write('{:<22d}   {:<11} {:}'.format(0, 'S4', '! Constant in the f curve best-fit for         AOA< alpha2; by definition it depends on the airfoil. [ignored if UAMod<>1]\n'))
                    f.write('{:<22f}   {:<11} {:}'.format(cn1, 'Cn1', '! Critical value of C0n at leading edge separation. It should be extracted from airfoil data at a given Mach and Reynolds number. It can be calculated from the static value of Cn at either the break in the pitching moment or the loss of chord force at the onset of stall. It is close to the condition of maximum lift of the airfoil at low Mach numbers.\n'))
                    f.write('{:<22f}   {:<11} {:}'.format(cn2, 'Cn2', '! As Cn1 for negative AOAs.\n'))
                    f.write('{!s:<22}   {:<11} {:}'.format('default', 'St_sh', "! Strouhal's shedding frequency constant.  [default = 0.19]\n"))
                    f.write('{:<22f}   {:<11} {:}'.format(cd0, 'Cd0', '! 2D drag coefficient value at 0-lift.\n'))
                    f.write('{:<22f}   {:<11} {:}'.format(cm0, 'Cm0', '! 2D pitching moment coefficient about 1/4-chord location, at 0-lift, positive if nose up. [If the aerodynamics coefficients table does not include a column for Cm, this needs to be set to 0.0]\n'))
                    f.write('{:<22d}   {:<11} {:}'.format(0, 'k0', '! Constant in the \\hat(x)_cp curve best-fit; = (\\hat(x)_AC-0.25).  [ignored if UAMod<>1]\n'))
                    f.write('{:<22d}   {:<11} {:}'.format(0, 'k1', '! Constant in the \\hat(x)_cp curve best-fit.  [ignored if UAMod<>1]\n'))
                    f.write('{:<22d}   {:<11} {:}'.format(0, 'k2', '! Constant in the \\hat(x)_cp curve best-fit.  [ignored if UAMod<>1]\n'))
                    f.write('{:<22d}   {:<11} {:}'.format(0, 'k3', '! Constant in the \\hat(x)_cp curve best-fit.  [ignored if UAMod<>1]\n'))
                    f.write('{:<22d}   {:<11} {:}'.format(0, 'k1_hat', '! Constant in the expression of Cc due to leading edge vortex effects.  [ignored if UAMod<>1]\n'))
                    f.write('{!s:<22}   {:<11} {:}'.format('default', 'x_cp_bar', '! Constant in the expression of \\hat(x)_cp^v. [ignored if UAMod<>1, default = 0.2]\n'))
                    f.write('{!s:<22}   {:<11} {:}'.format('default', 'UACutout', '! Angle of attack above which unsteady aerodynamics are disabled (deg). [Specifying the string "Default" sets UACutout to 45 degrees]\n'))
                    f.write('{!s:<22}   {:<11} {:}'.format('default', 'filtCutOff', '! Cut-off frequency (-3 dB corner frequency) for low-pass filtering the AoA input to UA, as well as the 1st and 2nd derivatives (Hz) [default = 20]\n'))

                f.write('!........................................\n')
                f.write('! Table of aerodynamics coefficients\n')
                f.write('{:<22d}   {:<11} {:}'.format(PolarTable.shape[0], 'NumAlf', '! Number of data lines in the following table\n'))
                f.write('!    Alpha      Cl      Cd        Cm\n')
                f.write('!    (deg)      (-)     (-)       (-)\n')
                for row in PolarTable:
                    f.write(' '.join(['{: 2.14e}'.format(val) for val in row])+'\n')
        f.close()
       
                
    

if __name__ == '__main__':
    pass




