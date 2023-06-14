#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Oct 18 10:24:26 2022

2023-02-02 began modifications to use new REBOUND results per Alessi's README.txt

@author: rjedicke
"""

import numpy as np
import matplotlib.pyplot as pyplot

import os

import HUTILITIES as h
import PLOT_UTILITIES as plot
import utilities as util
# import pyoorb as oo
# import orbit


# oo.pyoorb.oorb_init()

import glob
# ephfile = os.path.join( os.getenv('OORB_DATA'), 'de430.dat' )
# oo.pyoorb.oorb_init( ephfile )


from scipy.constants import pi
from scipy.optimize import curve_fit , newton
from scipy.integrate import quad , dblquad , tplquad

moon_axis_km = 384400 # https://nssdc.gsfc.nasa.gov/planetary/factsheet/moonfact.html

global moon_escape_speed_kps
moon_escape_speed_kps = 2.38

dirSimulation = '../data/'

   
# ejecta initial conditions
global giPar, gLon_rad, gLat_rad, gv0_kps, gAzimuth_rad, gImpactTime_jd, gxv0_ssb_ec, gOrbit0_hel, gOrbit0_geo
global gnParticlesPerCrater
gnParticlesPerCrater = 78

# capture info
global jcrater, jparticle
global t_days, EwrtEarth, EwrtEMBary, xvEjecta_ssb_ec, xvEarth_ssb_ec, xvMoon_ssb_ec, orbit_hel, orbit_geo
global gtBegin_day, gtEnd_day, giParticleID

# escape info
global gEscape_ipar, gEscape_day

# collision info
global gCollide_day, gCollide_ipar, gCollide_body

global gSummaryTable

global gMinDiameter_m
gMinDiameter_m = 1 

global ejecta_speed_max, ejecta_speed_min
ejecta_speed_max = 5000 
ejecta_speed_min = 0

def fraction_captured ( x , a , x_max):
    global moon_escape_speed_kps
    return heavyside(x , moon_escape_speed_kps) * maxwell_boltzmann(x, a , x_max)

def fraction_captured_model ( v , c , a ):
    return c * (v ** a)

def maxwell_boltzmann ( x , a , x_max ) : 
    return  (x-x_max)**2 * np.exp(-((x-x_max)**2) / (2 * a**2) ) / a**3

def heavyside ( x , x_0 ):
    return np.heaviside ( x - x_0 , 1 )
    # return y

def polyfit ( x , y , deg=10 ):
    return np.polyfit ( x , y , deg )

def Gaussian ( x , mean , std , scale ):
    return np.exp ( -.5 * ((x-mean) / std)**2  ) *scale

def lifetime ( x , mean , std , scale ) : 
    global moon_escape_speed_kps
    return Gaussian ( x , mean , std , scale ) * heavyside ( x , moon_escape_speed_kps )

#--------------------------------------------------------------------------------------------------------------------------------------------------------------
def createSummaryTable( nCraters=60 ):

    global gnParticlesPerCrater
    global gSummaryTable
    global giPar, gLon_rad, gLat_rad, gv0_kps, gAzimuth_rad, gImpactTime_jd, gxv0_ssb_ec, gOrbit0_hel, gOrbit0_geo


    loadInitialConditions()

    gSummaryTable = np.full(  (nCraters*gnParticlesPerCrater+1), -1,
                                                    dtype=[ ('craterID',  'int'), 
                                                 ('lon_rad',   'float'),
                                                 ('lat_rad',   'float'),
                                                 ('particleID','int'), 
                                                 ('v0_kps',    'float'),
                                                 ('t0_jd',    'float'),
                                                 ('azimuth_deg',    'float'),
                                                 ('captured' , 'int'),
                                                 ('begin_day' , 'object'),
                                                 ('end_day' , 'object'),
                                                 ] )
    for jCrater in np.arange( nCraters )+1:

        if jCrater==55: continue

        loadCaptureSummary(jCrater)
        global giParticleID, gtBegin_day, gtEnd_day

        # giParticleID_to_particleID = (jCrater-1) *gnParticlesPerCrater + giParticleID

        for iParticle in np.arange( gnParticlesPerCrater )+1:
        
            particleID = ( jCrater - 1 ) * gnParticlesPerCrater + iParticle

            gSummaryTable[particleID][ 'craterID'   ] = jCrater
            gSummaryTable[particleID][ 'lon_rad' ] = gLon_rad[ particleID -1 ]
            gSummaryTable[particleID][ 'lat_rad' ] = gLat_rad[ particleID -1 ]
            gSummaryTable[particleID][ 'particleID' ] = particleID
            gSummaryTable[particleID][ 'v0_kps' ] = gv0_kps[ particleID -1  ]
            gSummaryTable[particleID][ 't0_jd' ] = gImpactTime_jd[ particleID -1  ]
            gSummaryTable[particleID][ 'azimuth_deg' ] = np.rad2deg(gAzimuth_rad[ particleID -1  ])

            captured   = np.where(giParticleID==particleID)

            if len(captured[0])>0  : 
                gSummaryTable[particleID][ 'captured'  ] = 1 
                gSummaryTable[particleID][ 'begin_day' ] = gtBegin_day[captured]
                gSummaryTable[particleID][ 'end_day'   ] = gtEnd_day[captured]

# TODO: prompt captures vs delayed captures
#--------------------------------------------------------------------------------------------------------------------------------------------------------------
def fractionCapturedVSVelocity ( show=False ):
    global gSummaryTable 

    unique_v0_kps   = np.unique(gSummaryTable['v0_kps'] [ np.where (gSummaryTable['v0_kps'] > 0 ) ] )
    total_per_v0    = np.zeros(len(unique_v0_kps))
    captured_per_v0 = np.zeros(len(unique_v0_kps))
    prompt_per_v0   = np.zeros(len(unique_v0_kps))
    delayed_per_v0  = np.zeros(len(unique_v0_kps))

    for ii in range(len(unique_v0_kps)):
    # for ii in range(1,2):
        v_0 = unique_v0_kps[ii]
        v_0_matches = gSummaryTable[np.where(gSummaryTable['v0_kps']==v_0)]
        
        total_per_v0[ii]    = len(v_0_matches)
        captured            = np.where( v_0_matches['captured'] == 1 )

        captured_per_v0[ii] = len(captured[0])
        
        captured_begin_days = v_0_matches[captured] ['begin_day']

        for jj in range ( len(captured_begin_days) ):
            if np.any ( captured_begin_days[jj]<10 ) : 
                prompt_per_v0  [ii] += 1
            if np.any ( captured_begin_days[jj]>10 ) :
                delayed_per_v0 [ii] += 1
        # for jj in range( len(captured[0]) ):
        #     print ( v_0_matches[captured[0][jj]]['begin_day'] )
            # prompt_per_v0  [ii] += len(v_0_matches[captured]['begin_day'][jj]  < 1)
            # delayed_per_v0 [ii] += len(v_0_matches[captured]['begin_day'][jj]  > 1)
        # prompt_per_v0   = len( np.where(v_0_matches[captured]['begin_day']  < 1)[0] )
        # delayed_per_v0  = len( np.where(v_0_matches[captured]['begin_day']  > 1)[0] )

    fraction_per_v0 = captured_per_v0 / total_per_v0
    # prompt_over_delayed = prompt_per_v0 / delayed_per_v0
    param , param_cov = curve_fit ( fraction_captured_model , unique_v0_kps, fraction_per_v0 , p0=[1046 , -8] )

    p_param , p_param_cov = curve_fit ( fraction_captured_model , unique_v0_kps, prompt_per_v0/total_per_v0 , p0=[1046 , -8] )

    d_param , d_param_cov = curve_fit ( fraction_captured_model , unique_v0_kps, delayed_per_v0/total_per_v0 , p0=[1046 , -8] )

    if show: 

        print('unique initial velocities [km/s]: ' , unique_v0_kps)
        print('total per v0: ' , total_per_v0) 
        print('total captured: ' , captured_per_v0) 
        print('delayed capture:' ,delayed_per_v0 )
        print('fraction captured: ' , fraction_per_v0) 

        fig , ax = pyplot.subplots( )

        ax.set_xlabel( 'Ejection speed - lunar escape speed [km/s]' )
        ax.set_ylabel('Fraction captured' )
        ax.set_title('Fraction captured vs ejection speed')
        ax.set_xlim (1e-2 , 5)
        ax.set_ylim (1e-2 , 1.01)
        
        print('fitted parameters : ', param)


        # ax.plot ( unique_v0_kps , fraction_captured(unique_v0_kps , *param))
        ax.plot ( np.arange(0,5,.001) - moon_escape_speed_kps, fraction_captured_model(np.arange(0,5,.001) , *param) , color='red')
        # ax.plot ( np.arange(0,10,.1) , fraction_captured(np.arange(0,10,.1) , 1 , 2.6))

        ax.set_xscale('log')
        ax.set_yscale('log')


        plot.plot2d ( unique_v0_kps-moon_escape_speed_kps , prompt_per_v0/total_per_v0  ,  markersize=5 , xrange=(1e-2,2) , yrange=(1e-2,1.01) , title='prompt Fraction captured vs ejection speed' , xlabel='Ejection speed - lunar escape [km/s]' , ylabel='prompt Fraction captured' )
        plot.plot2d ( unique_v0_kps-moon_escape_speed_kps , delayed_per_v0/total_per_v0 ,  markersize=5  , xrange=(1e-2,2) , yrange=(1e-2,1.01) , title='delayed Fraction captured vs ejection speed' , xlabel='Ejection speed- lunar escape [km/s]' , ylabel='delayed Fraction captured' )
        
        ax.scatter ( unique_v0_kps-moon_escape_speed_kps , prompt_per_v0/total_per_v0  , label='prompt Fraction captured '  , color='green' )
        ax.scatter ( unique_v0_kps-moon_escape_speed_kps , delayed_per_v0/total_per_v0 , label='delayed Fraction captured ' ,color='blue')


        ax.plot ( np.arange(0,5,.001) - moon_escape_speed_kps, fraction_captured_model(np.arange(0,5,.001) , *p_param) , color='green')
        ax.plot ( np.arange(0,5,.001) - moon_escape_speed_kps, fraction_captured_model(np.arange(0,5,.001) , *d_param) , color='blue')
        
        ax.scatter ( unique_v0_kps-moon_escape_speed_kps , fraction_per_v0 ,label='total fraction captured' , color='red')

        ax.legend()

    return unique_v0_kps , fraction_per_v0 , param

# TODO prompt and delayed captures
#--------------------------------------------------------------------------------------------------------------------------------------------------------------
def captureLifetimeVSVelocity ( show=True ) :
    global gSummaryTable

    unique_v0_kps   = np.unique(gSummaryTable['v0_kps'][ np.where (gSummaryTable['v0_kps'] > 0 ) ])
    captured_per_v0 = np.zeros(len(unique_v0_kps))

    lifetime_per_v0 = np.zeros(len(unique_v0_kps))

    prompt_per_v0   = np.zeros(len(unique_v0_kps))
    delayed_per_v0  = np.zeros(len(unique_v0_kps))

    prompt_lifetime  = np.zeros(len(unique_v0_kps))
    delayed_lifetime = np.zeros(len(unique_v0_kps))

    for ii in range(len(unique_v0_kps)):
        v_0 = unique_v0_kps[ii]
        v_0_matches = gSummaryTable[np.where(gSummaryTable['v0_kps']==v_0)]

        captured  = np.where( v_0_matches['captured'] == 1 )
        captured_per_v0[ii] = len(captured[0])

        begin_day = v_0_matches[captured]['begin_day']
        end_day   = v_0_matches[captured]['end_day'  ]

        for jj in range(len(captured[0])):
            prompt  = np.where(begin_day [jj] < 10)
            delayed = np.where(begin_day [jj] > 10)

            lifetime_per_v0[ii] += np.sum ( end_day[jj] - begin_day[jj] )

            prompt_lifetime  [ii] += np.sum ( end_day[jj][prompt ] - begin_day[jj][prompt ] ) 
            delayed_lifetime [ii] += np.sum ( end_day[jj][delayed] - begin_day[jj][delayed] ) 

            prompt_per_v0 [ii]  += len(prompt[0])
            delayed_per_v0 [ii] += len(delayed[0])


    avg_capture_lifetime = lifetime_per_v0 / captured_per_v0 
    param , param_cov = curve_fit ( lifetime , unique_v0_kps[np.where(lifetime_per_v0>0)] , avg_capture_lifetime[np.where(lifetime_per_v0>0)] )


    if show: 
        print('unique initial velocities [km/s]: ' , unique_v0_kps)
        print('lifetime per v0 [days] : ' ,         lifetime_per_v0 ) 

        print(captured_per_v0)


        print('prompt lifetime per v0 [days] : ' , prompt_lifetime  )
        print(prompt_per_v0)

        print('delayed lifetime per v0 [days] : ' , delayed_lifetime)
        print(delayed_per_v0)


        # plot.plot2d ( unique_v0_kps , lifetime_per_v0 / captured_per_v0 , markersize=5 , logx=True , logy=True , xrange=(2,5) , yrange=(.1, 10000), title='Average Capture lifetime vs ejection speed' , xlabel='Ejection speed [km/s]' , ylabel='Average Capture lifetime [days]' )
    
    

        fig , ax = pyplot.subplots( )

        ax.scatter ( unique_v0_kps , avg_capture_lifetime )
        ax.set_xlabel('Ejection speed [km/s]' )
        ax.set_ylabel('Average Capture lifetime [days]')
        ax.set_title('Average Capture lifetime vs ejection speed')
        ax.set_xlim (1 , 25)
        ax.set_ylim (1e-1 , 10e4)
    

        ax.plot ( np.arange(0,10,.001) , lifetime(np.arange(0,10,.001) , *param) )

        ax.set_xscale('log')
        ax.set_yscale('log')




        plot.plot2d ( unique_v0_kps , prompt_lifetime / prompt_per_v0   , markersize=5 , logx=True , logy=True , xrange=(2,5) , title='Average prompt Capture lifetime vs ejection speed' , xlabel='Ejection speed [km/s]' , ylabel='Average prompt Capture lifetime [days]' )
        plot.plot2d ( unique_v0_kps , delayed_lifetime / delayed_per_v0  , bShow=False , markersize=5 , logx=True , logy=True , xrange=(2,5) , title='Average delayed Capture lifetime vs ejection speed' , xlabel='Ejection speed [km/s]' , ylabel='Average delayed Capture lifetime [days]' )

        print()


    return unique_v0_kps , avg_capture_lifetime , param

def interpolate_polynomials ( small_impactor_param , large_impactor_param , impactor_diameter_km ) :
    interp = np.zeros ( (len(impactor_diameter_km) , len(small_impactor_param)) )
    for ii in range ( len(impactor_diameter_km) ) :
        for jj in range (  len(small_impactor_param) ) :
            interp[ii ,jj] = small_impactor_param[jj] + (large_impactor_param[jj]-small_impactor_param[jj]) * ( impactor_diameter_km[ii] - .100 ) / ( 72 - .100 ) 
    return interp

#--------------------------------------------------------------------------------------------------------------------------------------------------------------
def cumulativeLunarImpactRate_Myr (   impactor_speed_kps=1 , diameter_km=1 , param=0 , p=2.5 , A=1737.4/6378 ):
    # cumulative rate of impactors larger than 1km per 1 million years
    rateOnEarth_1km_1Myr = 1
    rateOnMoon_1km_1Myr  = rateOnEarth_1km_1Myr*(A)**2
    # print(param)

    # dividing by 45 makes it diferential in v... ? but still in D
    integrated_prob = quad(param , 0 , 45)
    probability = param(impactor_speed_kps) / integrated_prob[0]
    if probability < 1e-3: probability = 0
    # probability = 1
    # NOTE 2/23
    # equation 10: N = c * d **(-p)
    # p = 2.5
    return probability * rateOnMoon_1km_1Myr * diameter_km ** (-p)

def returnIntegrated_impactRate ( impactor_diameter_km , param=(0) , p=2.5 , A=1737.4/6378 ):
    I = quad ( cumulativeLunarImpactRate_Myr , 0 , 50 , args=(impactor_diameter_km , param , p , A) )
    return I[0]

#--------------------------------------------------------------------------------------------------------------------------------------------------------------
def diffLunarImpactRate_Myr (   impactor_speed_kps=1 , diameter_km=1 , param=0 , p=2.5 , A=1737.4/6378 ):
    # cumulative rate of impactors larger than 1km per 1 million years
    rateOnEarth_1km_1Myr = 1
    rateOnMoon_1km_1Myr  = rateOnEarth_1km_1Myr*(A)**2
    # print(param)

    # dividing by 45 makes it diferential in v... ? but still in D
    integrated_prob = quad(param , 0 , 45)
    probability = param(impactor_speed_kps) / integrated_prob[0]
    if probability < 1e-3: probability = 0
    # probability = 1
    # NOTE 2/23
    # equation 10: N = c * d **(-p)
    # p = 2.5
    return p*probability * rateOnMoon_1km_1Myr * diameter_km ** (-p-1)

def return_diffIntegrated_impactRate ( impactor_diameter_km , param=(0) , p=2.5 , A=1737.4/6378 ):
    I = quad ( diffLunarImpactRate_Myr , 0 , 50 , args=(impactor_diameter_km , param , p , A) )
    return I[0]


def nEjectaDensity_ ( impactor_diameter_km , ejecta_diameter_km , impactor_speed_kps , show=False , p=2.5):
    # 1km impactor creates 10km crater
    # ratio of crater diameter to depth = 10/1
    # assuming hemispherical bowl : volume of ejecta
    # volume of ejecta has same size distribution as impactors : size freq of ejecta
    global moon_escape_speed_kps
    R = .2

    ejecta_d_min_m = 1e-6 

    # ejecta diameter vs speed: assume random and uniform
    # crater_diameter_m = x * impactor_diameter_km * 1000

    impactor_speed_mps = impactor_speed_kps * 1000

    moon_diameter_m = 3_474_200
    moon_mass_kg    = 7.34767309e22

    rho_impactor_kg_m = 1700
    rho_moon_kg_m     = 2550 
    impactor_KE       = 1/12 * np.pi * rho_impactor_kg_m * (impactor_speed_mps ** 2) * ((impactor_diameter_km*1000)**3)
    # impactor_KE       = 1/12* np.pi * rho_impactor_kg_m * (impactor_speed_mps ** 2) * ((impactor_diameter_km*1000)**3)

    impactor_radius_m = 1000*impactor_diameter_km/2

    g_impactor = 6.67e-11 * (rho_impactor_kg_m * 4/3 * np.pi * impactor_radius_m**3) * impactor_radius_m**-2
    g_moon     = 6.67e-11 * (rho_moon_kg_m * 4/3 * np.pi * (moon_diameter_m/2)**3) * (moon_diameter_m/2)**-2

    sin45 = (2**.5)/2
    crater_diameter_m = (2.7e-2) *  (rho_impactor_kg_m ** (1/6)) * (rho_moon_kg_m**(-.5)) * (impactor_KE**.28) * ( 1 - .0095 * ( 1 - sin45 )) * (( g_moon / g_impactor ) ** (3/16))

    crater_depth_m    = R * crater_diameter_m 
    # volume of hemispherical bowl
    # crater_volume_m   = np.pi/3 * (crater_diameter_m/2 - crater_depth_m) * crater_depth_m**2

    crater_volume_m = R * np.pi/6 * ( .75 + R**2 ) * (crater_diameter_m**3)

    x = (3-p)/p * R * (.75 + R**2)

    LHS = lambda c : c * np.abs(c ** (3/p - 1) - ejecta_d_min_m ** (3-p) ) 
    RHS = x * crater_diameter_m**3

    lambda_C = lambda c : LHS(c) - RHS


    r = np.linspace(5000 , RHS , 500)


    C_first_guess = RHS ** (1/1.2)
    
    # C = r[np.where ( lambda_C(r) < .1 )][0]
    # C = r[np.where ( np.abs(lambda_LHS(r) - RHS) < 1000 )][0]
    C = newton ( lambda_C , C_first_guess ) 
    print(f'C={C}   impactorD_km={impactor_diameter_km}   impactorV_kps={impactor_speed_kps}')
    
    ejecta_d_max_m = C ** (1/p)
    ejecta_volume_m = np.pi * C * p / (6 * (3-p)) * np.abs ( ejecta_d_max_m**(3-p) - ejecta_d_min_m**(3-p) )
    print ( f'crater_volume={crater_volume_m}  ejecta_volume={ejecta_volume_m}    crater-ejecta_V={np.abs(crater_volume_m-ejecta_volume_m)}  ' )
    print()
    # N = C * p * (ejecta_diameter_km*1000) ** (-p-1)
    N = C * (ejecta_diameter_km*1000) ** (-p)


    ejecta_speed_kps = velocity_kps_vs_impactorD_ejectaD ( ejecta_diameter_km/impactor_diameter_km )

    try: 
        N[ejecta_speed_kps < moon_escape_speed_kps] = 0
    except Exception as e: 
        if ejecta_speed_kps < moon_escape_speed_kps: N = 0


    return N

def returnIntegrated_n () : pass


# FUNCTION OF impactor diameter , ejecta diameter , ejecta speed
#--------------------------------------------------------------------------------------------------------------------------------------------------------------
def craterDiameter_km ( impactor_diameter_km , ejecta_diameter_km , impactor_speed_kps , show=False , p=2.5):
    # 1km impactor creates 10km crater
    # ratio of crater diameter to depth = 10/1
    # assuming hemispherical bowl : volume of ejecta
    # volume of ejecta has same size distribution as impactors : size freq of ejecta
    global moon_escape_speed_kps
    R = .2


    # ejecta diameter vs speed: assume random and uniform
    # crater_diameter_m = x * impactor_diameter_km * 1000

    # impactor_speed_kps = np.array([5 , 10 , 20 , 30 , 40])
    # impactor_speed_mps = impactor_speed_kps * 1000

    moon_diameter_m = 3_474_200
    moon_mass_kg    = 7.34767309e22

    rho_impactor_kg_m = 1700
    rho_moon_kg_m     = 2550 
    # impactor_KE       = np.pi * rho_impactor_kg_m * (impactor_speed_kps ** 2) * ((impactor_diameter_km*1000)**3)
    impactor_KE       = 1/12* np.pi * rho_impactor_kg_m * ((impactor_speed_kps*1000) ** 2) * ((impactor_diameter_km*1000)**3)


    g_impactor = 6.67e-11 * (rho_impactor_kg_m * 4/3 * np.pi * (1000*impactor_diameter_km/2)**3) * (1000*impactor_diameter_km/2)**-2
    g_moon     = 6.67e-11 * (rho_moon_kg_m * 4/3 * np.pi * (moon_diameter_m/2)**3) * (moon_diameter_m/2)**-2

    crater_diameter_m = (2.7e-2) *  (rho_impactor_kg_m ** (1/6)) * (rho_moon_kg_m**(-.5)) * (impactor_KE**.28) * ( 1 - .0095 * ( 1 - (2**.5)/2 )) * (( g_moon / g_impactor ) ** (3/16))

    return crater_diameter_m / 1000

# returns polynomial coefficients for marchi2009 and yue2013 figs 
def digitize_plots ( poly_order = 8 , show=False ):
    marchi2009_small = np.loadtxt ( '../data/digitized-figs/plot-marchi2009-small.csv' , skiprows=1 , delimiter=',' )
    marchi2009_large = np.loadtxt ( '../data/digitized-figs/plot-marchi2009-large.csv' , skiprows=2 , delimiter=',' )
    yue2013          = np.loadtxt ( '../data/digitized-figs/plot-yue2013.csv' , skiprows=1 , delimiter  =',' )

    marchi_small_param = polyfit ( marchi2009_small[:,0] , marchi2009_small[:,1] , poly_order )
    
    marchi_large_param = polyfit ( marchi2009_large[:,0] , marchi2009_large[:,1] , poly_order )
    
    yue_param = polyfit ( yue2013[:,0] , yue2013[:,1] , poly_order )
        
    if show: 
        fig_marchi_sm , ax_marchi_sm = pyplot.subplots()
        fig_marchi_lg , ax_marchi_lg = pyplot.subplots()
    
        ax_marchi_sm.scatter ( marchi2009_small[:,0] , marchi2009_small[:,1] , label='marchi2009 100 meter impactor' )
        ax_marchi_lg.scatter ( marchi2009_large[:,0] , marchi2009_large[:,1] , label='marchi2009 72 km impactor' , color='red' )
        ax_marchi_sm.set_xlabel('impactor speed [km/s]')
        ax_marchi_sm.set_ylabel('probability')
        ax_marchi_lg.set_xlabel('impactor speed [km/s]')
        ax_marchi_lg.set_ylabel('probability')
        print( 'marchi2009 small impactor params:', marchi_small_param)

        print( 'marchi2009 large impactor params:', marchi_large_param)
        impactor_speed_kps = np.linspace ( 0,45 , 100 )
        ax_marchi_lg.plot ( impactor_speed_kps , np.polyval ( marchi_large_param , impactor_speed_kps  ) , label='72km best fit' )
        ax_marchi_sm.plot ( impactor_speed_kps , np.polyval ( marchi_small_param , impactor_speed_kps  ) , label='100m best fit' )
        ax_marchi_lg.legend()
        ax_marchi_sm.legend()

        fig_yue , ax_yue = pyplot.subplots()
        ax_yue.scatter ( yue2013[:,0] , yue2013[:,1] , label='Digitized data' )

        ax_yue.plot ( yue2013[:,0] , np.polyval ( yue_param , yue2013[:,0] ) , label='best fit' )
        ax_yue.set_xlabel('impactor speed [km/s]')
        ax_yue.set_ylabel('encounters')
        ax_yue.legend()
        # print('yue2013 params: ' , yue_param)

    return marchi_small_param , marchi_large_param , yue_param

def velocity_kps_vs_impactorD_ejectaD ( ratio_ejectaD_impactorD ) :
    return 10 ** ( (2.34 - .52 * np.log10(ratio_ejectaD_impactorD) ) ) / 1000

def integrand ( impactor_D_km , ejecta_speed_kps , impactor_speed_kps , param_1 , param_2 , param_3 , ejecta_D_km ) : 

    global ejecta_speed_max, ejecta_speed_min

    f = fraction_captured ( ejecta_speed_kps , *param_1)
    l = lifetime          ( ejecta_speed_kps , *param_2)

    F = cumulativeLunarImpactRate_Myr ( impactor_D_km, impactor_speed_kps , param_3 ) * (365 * 1e6 )
    N = nEjectaDensity_ (impactor_D_km , ejecta_D_km , impactor_speed_kps) / (ejecta_speed_max - ejecta_speed_min) / 45 

    return f * l * F * N

# integrate previous functions to calculate stead state population 
#--------------------------------------------------------------------------------------------------------------------------------------------------------------
def minimoons( ejecta_diameter_km ):
    global gMinDiameter_m

    unique_v0_kps , fraction_per_v0 , param_1 = fractionCapturedVSVelocity(show=False)
    marchi_small_param , marchi_large_param , yue_param = digitize_plots (  )


    _, avg_lifetime_per_v0 , param_2 = captureLifetimeVSVelocity(show=False)

    # integral = dblquad ( integrand , 0 , np.inf , 0 , np.inf , args=(param_1 , param_2 , ejecta_diameter_km) )
    integral = tplquad ( integrand , 0 , 10 , 0 , 10 , 0 , 45 ,  args=(param_1 , param_2 , marchi_small_param , ejecta_diameter_km) )

    return integral
#--------------------------------------------------------------------------------------------------------------------------------------------------------------
def loadCollisions( iMaxCraterID=60 ):

    print( '\n' )
    print( 'Loading collisions...' )    
    
    #1-MERCURY,2-VENUS, 3-EARTH, 4-MARS, 5-JUPITER, 6-SATURN, 7-URANUS, 8-NEPTUNE, 10-MOON, 11-SUN
    global gCollide_day, gCollide_ipar, gCollide_body
    
    gCollide_day  = np.empty( 0, dtype='float' )
    gCollide_ipar = np.empty( 0, dtype='int' )
    gCollide_body = np.empty( 0, dtype='int' )

    for jCrater in np.arange(iMaxCraterID)+1:
        
        craterDir = dirSimulation + 'crater' + str(jCrater)
        
        if( os.path.exists(craterDir) ):
            
            collideFilename = craterDir+'/collision_'+str(jCrater)+'.out'
            
            if( os.path.getsize(collideFilename) > 0 ):
                
                tday, jpar, ibody = np.loadtxt( collideFilename, dtype="float,int,int", unpack=True )
                
                ipar = ( jCrater-1 )*gnParticlesPerCrater + jpar
                
                gCollide_day  = np.append( gCollide_day,   tday )
                gCollide_ipar = np.append( gCollide_ipar,  ipar )
                gCollide_body = np.append( gCollide_body,  ibody )

    print( 'found ', len(gCollide_ipar), ' collisions' )

    h.plot( np.log10( gCollide_day/365.25/100 ), xlabel='$\log_{10}$( time of collision / [century] )' )

    fig, ax = pyplot.subplots( figsize=(8,8), dpi=200 )

    nPerBin, binEdges, patches = ax.hist( gCollide_body, bins=11, histtype='step', range=(0.5,11.5) )
    
    pyplot.xticks( ticks=np.linspace(1,11,11), labels=('MERCURY','VENUS','EARTH','MARS','JUPITER','SATURN','URANUS','NEPTUNE','PLUTO','MOON','SUN'), rotation=30 )
    pyplot.xlabel( 'body' ) 
    pyplot.ylabel( 'number' ) 
    pyplot.yscale( "log" )

    # pyplot.show()





#--------------------------------------------------------------------------------------------------------------------------------------------------------------
def loadEscapes( iMaxCraterID=60 ):

    print( '\n' )
    print( 'Loading escapes...' )    

    global gEscape_ipar, gEscape_day
    
    gEscape_ipar  = np.empty( 0, dtype='int' )
    gEscape_day = np.empty( 0, dtype='float' )

    for jCrater in np.arange(iMaxCraterID)+1:
        
        craterDir = dirSimulation + 'crater' + str(jCrater)
        
        if( os.path.exists(craterDir) ):
            
            escapeFilename = craterDir+'/escape_'+str(jCrater)+'.out'
            
            if( os.path.getsize(escapeFilename) > 0 ):
                
                jpar, t = np.loadtxt( escapeFilename, dtype="int,float", unpack=True )
                
                ipar = ( jCrater-1 )*gnParticlesPerCrater + jpar
                
                gEscape_ipar  = np.append( gEscape_ipar,  ipar )
                gEscape_day = np.append( gEscape_day, t )

    print( 'found ', len(gEscape_ipar), ' escapes' )

    h.plot( np.log10( gEscape_day/365.25/100 ), xlabel='$\log_{10}$( time of escape / [century] )' )





#--------------------------------------------------------------------------------------------------------------------------------------------------------------
def loadInitialConditions():

    # index of the fragment, crater longitude, latitude, initial ejecta speed from the Moon (km/s), azimuth, time (JD), 
    # SSB-centric particle state x,y,z,vx,vy,v,
    # heliocentric particle orbital elements a,e,inc,Omega,omega, 
    # geocentric particle orbital elements Ea, Ee, Einc, EOmega, Eomega
    # where here ID of the fragment goes from (icrater-1)*78+1 to (icrater-1)*78+78

    print( '\n' )
    print( 'Loading all initial conditions...' )
    
    szInitialConditionsFile = dirSimulation+'/initialconditions_all.out'

    global giPar, gLon_rad, gLat_rad, gv0_kps, gAzimuth_rad, gImpactTime_jd, gxv0_ssb_ec, gOrbit0_hel, gOrbit0_geo
    global gnParticlesPerCrater

    giPar, gLon_rad, gLat_rad, gv0_kps, gAzimuth_rad, gImpactTime_jd = np.loadtxt( szInitialConditionsFile, usecols=(0,1,2,3,4,5), dtype="int,float,float,float,float,float", unpack=True )

    gLon_rad = util.inRangeMinusPi2PlusPi( gLon_rad )
    
    x,y,z,vx,vy,vz = np.loadtxt( szInitialConditionsFile, usecols=( 6,7,8,9,10,11), unpack=True )
    
    gxv0_ssb_ec = np.column_stack( ( x,y,z,vx,vy,vz ) )
    
    a, e, i, O, w = np.loadtxt( szInitialConditionsFile, usecols=(12,13,14,15,16), unpack=True )
    
    gOrbit0_hel = np.column_stack( ( a, e, np.rad2deg(i), np.rad2deg(O), np.rad2deg(w) ) )
    
    del a, e, i, O, w
    
    a, e, i, O, w = np.loadtxt( szInitialConditionsFile, usecols=(17,18,19,20,21), unpack=True )
    
    gOrbit0_geo = np.column_stack( ( a, e, np.rad2deg(i), np.rad2deg(O), np.rad2deg(w) ) )
    
    del a, e, i, O, w
    
    gnParticles = len(giPar)
    
    print( 'loaded ', gnParticles, ' particles from ', gnParticles/gnParticlesPerCrater, ' craters' )
    
    h.plot( np.log10(gv0_kps),            xlabel='$\log_{10}$( ejection speed / [km/s] )' )
    h.plot( np.rad2deg(gAzimuth_rad),     xlabel='azimuth [deg]'  )
    
    h.plot(           gOrbit0_geo[:,0],   xlabel='geocentric semi-major axis [au]', logy=True )
    h.plot( np.log10( gOrbit0_geo[:,1] ), xlabel='$\log_{10}$( geocentric eccentricity )' )
    h.plot(           gOrbit0_geo[:,2],   xlabel='geocentric inclination [deg]', xrange=(0,180) )
    
    plot.plot2d( gOrbit0_geo[:,0], np.log10( gOrbit0_geo[:,1] ), xlabel='geocentric semi-major axis [au]', ylabel='$\log_{10}$( geocentric eccentricity )' )

    plot.plot2d( gv0_kps, np.log10( gOrbit0_geo[:,1] ), xlabel='ejection speed [kps]', xrange=(0,25), ylabel='$\log_{10}$( geocentric eccentricity )' )

    pyplot.figure( figsize=(8,4), dpi=200 )
    pyplot.subplot( projection="aitoff" )
    pyplot.title( "crater locations on moon" )
    pyplot.grid( True )
    pyplot.scatter( gLon_rad, gLat_rad, marker = '*', color = 'red', s = 2 )
  
    




#--------------------------------------------------------------------------------------------------------------------------------------------------------------
def loadCaptureSummary( icrater=53, size=(8,8), dpi=200, markerstyle='.', markersize=1 ):
    
    print( 'loading capture summary from crater ', icrater )
    
    captureSummaryFile = dirSimulation+'crater' + str(icrater) + '/capture_' + str(icrater) + '.out'

    global giParticleID, gtBegin_day, gtEnd_day
    iEjecta, gtBegin_day, gtEnd_day = np.loadtxt( captureSummaryFile, unpack=True )

    giParticleID = ( icrater - 1 ) * gnParticlesPerCrater + iEjecta


    h.plot( np.log10(gtEnd_day-gtBegin_day), xlabel='$\log_{10}$( capture duration / days )', nbins=20, xrange=(-1,4), title='crater '+str(icrater) )
        
    # capture start and stop times for each ejecta
    pyplot.figure( figsize=size, dpi=dpi )
    pyplot.ylim( (0,80) )
    pyplot.xscale( 'log' )
    for j in np.arange( len(giParticleID) ):
        x = [ gtBegin_day[j], gtEnd_day[j] ]
        y = [    giParticleID[j],  giParticleID[j] ]
        pyplot.plot( x, y, 'b', linestyle='-' )
    pyplot.xlabel( 'days from impact' )
    pyplot.ylabel( 'ejecta particle ID' )
    pyplot.title( 'crater ' + str(icrater) )
    # pyplot.show()
    




#--------------------------------------------------------------------------------------------------------------------------------------------------------------
def loadCapture( icrater=53, iparticle=4057 ):

    print( 'loading capture of particle ', iparticle, ' from crater ', icrater )

    global jcrater, jparticle
    global t_days, EwrtEarth, EwrtEMBary, xvEjecta_ssb_ec, xvEarth_ssb_ec, xvMoon_ssb_ec, orbit_hel, orbit_geo

    jcrater   = icrater
    jparticle = iparticle
    ipar      = 100000000 + iparticle

    captureDataFile = dirSimulation+'crater' + str(icrater) + '/capture_' + str(ipar) + '.out'

    t_days, EwrtEarth, EwrtEMBary = np.loadtxt( captureDataFile, usecols=(0,1,2), unpack=True )
    
    x,y,z,vx,vy,vz = np.loadtxt( captureDataFile, usecols=(3,4,5,6,7,8), unpack=True )
    
    xvEjecta_ssb_ec = np.column_stack( ( x,y,z,vx,vy,vz ) )
    
    del x,y,z,vx,vy,vz
     
    x,y,z,vx,vy,vz = np.loadtxt( captureDataFile, usecols=(9,10,11,12,13,14), unpack=True )
    
    xvEarth_ssb_ec = np.column_stack( ( x,y,z,vx,vy,vz ) )
    
    del x,y,z,vx,vy,vz
     
    x,y,z,vx,vy,vz = np.loadtxt( captureDataFile, usecols=(15,16,17,18,19,20), unpack=True )
    
    xvMoon_ssb_ec = np.column_stack( ( x,y,z,vx,vy,vz ) )
    
    del x,y,z,vx,vy,vz
    
    a, e, i, O, w = np.loadtxt( captureDataFile, usecols=(21,22,23,24,25), unpack=True )
    
    orbit_hel = np.column_stack( ( a, e, np.rad2deg(i), np.rad2deg(O), np.rad2deg(w) ) )
    
    del a, e, i, O, w
    
    a, e, i, O, w = np.loadtxt( captureDataFile, usecols=(26,27,28,29,30), unpack=True )
    
    orbit_geo = np.column_stack( ( a, e, np.rad2deg(i), np.rad2deg(O), np.rad2deg(w) ) )
    
    del a, e, i, O, w
    
    
    h.plot( t_days,     xlabel='time while captured [days]' )
    h.plot( EwrtEMBary, xlabel='energy wrt EM barycenter in crazy units' )

    util.printWarning( 'CURRENTLY DISPLAY *ALL* CAPTURES, NOT A SPECIFIC ONE' )
    
    plot.plot2d( t_days, EwrtEMBary,     xlabel='day', ylabel='energy wrt EM barycenter [crazy]', xrange=(t_days[0],t_days[-1]) )
    plot.plot2d( t_days, orbit_geo[:,0], xlabel='day', ylabel='geocentric semi-major axis [au]',  xrange=(t_days[0],t_days[-1]), yrange=(0,0.01) )
    plot.plot2d( t_days, orbit_geo[:,1], xlabel='day', ylabel='geocentric eccentricity',          xrange=(t_days[0],t_days[-1]), yrange=(0,1.00) )
    plot.plot2d( t_days, orbit_geo[:,2], xlabel='day', ylabel='geocentric inclination [deg]',     xrange=(t_days[0],t_days[-1]), yrange=(0,180.) )
    
    x_wrt_Moon = xvEjecta_ssb_ec[:,0] - xvMoon_ssb_ec[:,0]
    y_wrt_Moon = xvEjecta_ssb_ec[:,1] - xvMoon_ssb_ec[:,1]
    z_wrt_Moon = xvEjecta_ssb_ec[:,2] - xvMoon_ssb_ec[:,2]
    
    dMoon_au = np.sqrt( x_wrt_Moon**2 + y_wrt_Moon**2 + z_wrt_Moon**2 )
    
    plot.plot2d( t_days, dMoon_au,       xlabel='day', ylabel='lunacentric distance [au]',       xrange=(t_days[0],t_days[-1]), yrange=(0,0.01) )

