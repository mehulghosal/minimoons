#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jan 19 10:15:20 2021

@author: rjedicke
"""

from sty import fg
import numpy as np

import math
import subprocess
import utilities
import inspect

from scipy.constants import pi


twopi = 2 * pi



#------------------------------------------------------------------------------------
# fromn https://stackoverflow.com/questions/18425225/getting-the-name-of-a-variable-as-a-string
# and search for 'callers_local_vars'
# NOTE:  seems to only work in the namespace?
def getStringRepresentationOfThisVariable(var):
    callers_local_vars = inspect.currentframe().f_back.f_locals.items()
    return (str([k for k, v in callers_local_vars if v is var][0]))




#------------------------------------------------------------------------------------
def printWarning( szWarningMessage ):
    print( '\n' )
    print( fg.red + '*****************************************' + fg.rs )
    print( fg.red + szWarningMessage + fg.rs )    
    print( fg.red + '*****************************************' + fg.rs )
    print( '\n' )







#------------------------------------------------------------------------------------
# returns a numpy array of random values exponentially distributed in the xRange
def randomExponential( n=1, a=1, x=(0,1) ):
    
    f = np.random.uniform( 0.0, 1.0, n )
    
    return np.log( np.exp(a*x[1]) * f + np.exp(a*x[0]) * (1-f) ) / a 
    
    




#------------------------------------------------------------------------------------
# returns a numpy array of random values exponentially distributed in the xRange
def deltav_Tsiolkovsky( ISP_s=1, m_propellant_kg=2000, m_dry_kg=1000 ):
    
    g0_mps2 = 9.80665
    
    deltav_mps = ISP_s * g0_mps2 * np.log( 1 + m_propellant_kg / m_dry_kg )
    
    return deltav_mps





#------------------------------------------------------------------------------------
# solves law of cosines for length of a side of a triangle given the two other sides and the opposite angle
def lawOfCosines( l1, l2, angle12_deg ):
    return np.sqrt( l1**2 + l2**2 - 2 * l1 * l2 * np.cos( np.radians( angle12_deg ) ) )





#------------------------------------------------------------------------------------
# from https://stackoverflow.com/questions/1628386/normalise-orientation-between-0-and-360
# Normalizes any number to an arbitrary range [min,max)
# by assuming the range wraps around when going below min or above max 
def normalize( value, vmin, vmax ):
    width       = vmax  - vmin 
    offsetValue = value - vmin
    return ( offsetValue - ( math.floor( offsetValue / width ) * width ) ) + vmin





#------------------------------------------------------------------------------------
# rounds and truncates the value to N decimal places
def truncN( val, N ):
    return np.rint( val * 10**N ) / 10**N
 



#------------------------------------------------------------------------------------
# used Marc Paterno's efficiency calculator to determine the efficiency with min/max 
# values within the confidence interval specificed by CI
# i.e. pass in a vector of k successes in n trials with the CI
# returns a vector of efficiencies and min/max efficiencies for the interval
# The appropriate reference for the supporting document is for a
# Fermilab technical memo: FERMILAB-TM-2286-CD.

# returns efficiency, efficiency_min, efficiency_max

def efficiency( k, n, CI=0.68 ):
    
    # write a file containing the k, n vectors
    np.savetxt( 'junk.dat', np.column_stack([k,n]), fmt='%d %d' )    

    # call Paterno's function
    subprocess.check_output( '/Users/rjedicke/Dropbox/src/C++/calceff2/calceff2 junk.dat ' + str( CI ) + '> junk.results', shell=True )

    # format and return the output,  returns efficiency, efficiency_min, efficiency_max
    return np.loadtxt( 'junk.results', unpack=True )
     



#------------------------------------------------------------------------------------
# convert angle in radians to be in range [-pi,+pi)

def inRangeMinusPi2PlusPi( angle_rad ):
                                                                        
    return angle_rad - ( np.floor( ( angle_rad + pi ) / twopi ) ) * twopi
     



# #------------------------------------------------------------------------------------
# # convert angle in degrees to be in range [-180,+180) OR [-pi,+pi)
# # from https://stackoverflow.com/questions/2320986/easy-way-to-keeping-angles-between-179-and-180-degrees/2323034

# def inRange180open( angle ):
                                                                            
#     return angle - ( np.floor( ( angle + 180 ) / 360 ) ) * 360
     



# #------------------------------------------------------------------------------------
# # convert angle in degrees to be in range [-pi,+pi)
# # from https://stackoverflow.com/questions/2320986/easy-way-to-keeping-angles-between-179-and-180-degrees/2323034

# def inRange_pi_open( angle ):
                                                                            
#     return angle - ( np.floor( ( angle + pi ) / twopi ) ) * twopi
     



# #------------------------------------------------------------------------------------
# # convert angle in degrees to be in range [0,2pi)
# # from https://stackoverflow.com/questions/2320986/easy-way-to-keeping-angles-between-179-and-180-degrees/2323034

# def inRange_2pi_open( angle ):
                                                                            
#     return inRange_pi_open( angle ) + pi 



# #------------------------------------------------------------------------------------
# # convert angle in degrees to be in range (-180,+180] OR [-pi,+pi]
# # from https://stackoverflow.com/questions/2320986/easy-way-to-keeping-angles-between-179-and-180-degrees/2323034

# def inRange180closed( angle ):
#     return angle - np.ceil( angle / 360.0 - 0.5) * 360.0
     



#------------------------------------------------------------------------------------
# returns the parameters a, b, c in the quadratic form
# y = a ( x-x2 ) ( x-x3 ) + b ( x-x3 ) ( x-x1 ) + c ( x-x1 ) ( x-x2 ) 
# where the qaudratic goes through the 3 points (x1,y1), (x2,y2), (x3,y3)

def getQuadParams( x1, y1, x2, y2, x3, y3 ):
    
    a = y1 / ( x1 - x2 ) / ( x1 - x3 )
    b = y2 / ( x2 - x3 ) / ( x2 - x1 )
    c = y3 / ( x3 - x1 ) / ( x3 - x2 )

    return a, b, c    