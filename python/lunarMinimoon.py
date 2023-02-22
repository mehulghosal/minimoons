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

import pyoorb as oo

import HUTILITIES as h
import PLOT_UTILITIES as plot
import utilities as util
import pyoorb as oo
import orbit


oo.pyoorb.oorb_init()

import glob
# ephfile = os.path.join( os.getenv('OORB_DATA'), 'de430.dat' )
# oo.pyoorb.oorb_init( ephfile )


from scipy.constants import pi

moon_axis_km = 384400 # https://nssdc.gsfc.nasa.gov/planetary/factsheet/moonfact.html



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
def fractionCapturedVSVelocity ( show=True ):
    global gSummaryTable 

    unique_v0_kps   = np.unique(gSummaryTable['v0_kps'])
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
        for jj in range(len(captured[0])):
            prompt_per_v0  [ii] += len( np.where(v_0_matches[captured]['begin_day'] [jj] < 1)[0] )
            delayed_per_v0 [ii] += len( np.where(v_0_matches[captured]['begin_day'] [jj] > 1)[0] )

    fraction_per_v0 = captured_per_v0 / total_per_v0
    prompt_over_delayed = prompt_per_v0 / delayed_per_v0

    print('unique initial velocities [km/s]: ' , unique_v0_kps)
    print('total per v0: ' , total_per_v0) 
    print('total captured: ' , captured_per_v0) 
    print('fraction captured: ' , fraction_per_v0) 
    print('ratio of prompt over delayed captures: ', prompt_over_delayed)
    print()

    plot.plot2d ( unique_v0_kps , fraction_per_v0 ,  markersize=5 , xrange=(1,25) , yrange=(-.01,1) , title='Fraction captured vs ejection speed' , xlabel='Ejection speed [km/s]' , ylabel='Fraction captured' )
    plot.plot2d ( unique_v0_kps , fraction_per_v0 ,  markersize=5 , logx=True , logy=True , xrange=(1,25) , yrange=(0.001,1) , title='Fraction captured vs ejection speed' , xlabel='Ejection speed [km/s]' , ylabel='Fraction captured' )

    plot.plot2d ( unique_v0_kps , prompt_over_delayed , bShow=show ,  markersize=5 , logx=True , logy=True , xrange=(1,25) , yrange=(0.001,1) , title='Fraction captured vs ejection speed' , xlabel='Ejection speed [km/s]' , ylabel='Fraction captured' )

    return unique_v0_kps , fraction_per_v0

# TODO prompt and delayed captures
#--------------------------------------------------------------------------------------------------------------------------------------------------------------
def captureLifetimeVSVelocity ( show=True ) :
    global gSummaryTable

    unique_v0_kps   = np.unique(gSummaryTable['v0_kps'])
    lifetime_per_v0 = np.zeros(len(unique_v0_kps))

    prompt_per_v0   = np.zeros(len(unique_v0_kps))
    delayed_per_v0  = np.zeros(len(unique_v0_kps))

    for ii in range(len(unique_v0_kps)):
        v_0 = unique_v0_kps[ii]
        v_0_matches = gSummaryTable[np.where(gSummaryTable['v0_kps']==v_0)]

        captured  = np.where( v_0_matches['captured'] == 1 )

        begin_day = v_0_matches[captured]['begin_day']
        end_day   = v_0_matches[captured]['end_day'  ]

        for jj in range(len(captured[0])):
            prompt  = np.where(v_0_matches[captured]['begin_day'] [jj] < 1)
            delayed = np.where(v_0_matches[captured]['begin_day'] [jj] > 1)

            lifetime_per_v0[ii] += np.sum ( end_day[jj] - begin_day[jj] )

            prompt_per_v0  [ii] += np.sum ( end_day[jj][prompt ] - begin_day[jj][prompt ] )
            delayed_per_v0 [ii] += np.sum ( end_day[jj][delayed] - begin_day[jj][delayed] )


    print('unique initial velocities [km/s]: ' , unique_v0_kps)
    print('lifetime per v0 [days] : ' , lifetime_per_v0) 
    print('prompt lifetime per v0 [days] : ' , prompt_per_v0)
    print('delayed lifetime per v0 [days] : ' , delayed_per_v0)
    print()

    plot.plot2d ( unique_v0_kps , lifetime_per_v0 , markersize=5 , xrange=(1,25)  , title='Capture lifetime vs ejection speed' , xlabel='Ejection speed [km/s]' , ylabel='Capture lifetime [days]' )
    plot.plot2d ( unique_v0_kps , lifetime_per_v0 , markersize=5 , logx=True , logy=True , xrange=(1,25) , title='Capture lifetime vs ejection speed' , xlabel='Ejection speed [km/s]' , ylabel='Capture lifetime [days]' )
    
    plot.plot2d ( unique_v0_kps , prompt_per_v0  , markersize=5 , logx=True , logy=True , xrange=(1,25) , title='prompt Capture lifetime vs ejection speed' , xlabel='Ejection speed [km/s]' , ylabel='Capture lifetime [days]' )
    plot.plot2d ( unique_v0_kps , delayed_per_v0 , bShow=show , markersize=5 , logx=True , logy=True , xrange=(1,25) , title='delayed Capture lifetime vs ejection speed' , xlabel='Ejection speed [km/s]' , ylabel='Capture lifetime [days]' )

    return unique_v0_kps , lifetime_per_v0

# TODO
#--------------------------------------------------------------------------------------------------------------------------------------------------------------
def fluxDensity ( diameter=1 ):
    # cumulative rate of impactors larger than 1km per 1 million years
    rateOnEarth_1km_1Myr = 1
    rateOnMoon_1km_1Myr  = rateOnEarth_1km_1Myr*(1737.4/6378)**2

    # size distribution of objects function of diamter - read paper


# TODO
# FUNCTION OF impactor diameter , ejecta diameter , ejecta speed
#--------------------------------------------------------------------------------------------------------------------------------------------------------------
def nEjectaDensity ( diameter , ejecta_diameter , ejecta_speed ):
    # 1km impactor creates 10km crater
    # ratio of crater diameter to depth = 10/1
    # assuming hemispherical bowl : volume of ejecta
    # volume of ejecta has same size distribution as impactors : size freq of ejecta
    
    # ejecta diameter vs speed: assume random and uniform
    crater_diameter = 10 * diameter
    crater_depth    = crater_diameter / 10
    crater_volume   = np.pi/3 * (crater_diameter/2 - crater_depth) * crater_depth**2



    pass

# integrate previous functions to calculate stead state population 
#--------------------------------------------------------------------------------------------------------------------------------------------------------------
def minimoons(  ):
    global gMinDiameter_m
    pass

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

