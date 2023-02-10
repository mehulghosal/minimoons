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
import oorb
import orbit


oo.pyoorb.oorb_init()

import glob
ephfile = os.path.join( os.getenv('OORB_DATA'), 'de430.dat' )
oo.pyoorb.oorb_init( ephfile )


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

# escape info
global gEscape_ipar, gEscape_day

# collision info
global gCollide_day, gCollide_ipar, gCollide_body

global gSummaryTable





#--------------------------------------------------------------------------------------------------------------------------------------------------------------
def createSummaryTable( nCraters=60 ):

    global gnParticlesPerCrater    
    global gSummaryTable


    loadInitialConditions()

   #summaryTable = np.array( (nCraters*gnParticlesPerCrater+1,5), 
    gSummaryTable = np.full(  (nCraters*gnParticlesPerCrater+1,7), -1,
                                                    dtype=[ ('craterID',  'int'), 
                                                 ('lon_rad',   'float'),
                                                 ('lat_rad',   'float'),
                                                 ('particleID','int'), 
                                                 ('v0_kps',    'float'),
                                                 ('t0_jd',    'float'),
                                                 ('azimuth_deg',    'float'),
                                                 ] )

    for jCrater in np.arange( nCraters )+1:
        
        for iParticle in np.arange( gnParticlesPerCrater )+1:
        
            particleID = ( jCrater - 1 ) * gnParticlesPerCrater + iParticle

            gSummaryTable[ 'craterID'   ][particleID] = jCrater
            gSummaryTable[ 'lon_rad' ][particleID] = gLon_rad[ giPar==particleID ]
            gSummaryTable[ 'lat_rad' ][particleID] = gLat_rad[ giPar==particleID ]
            gSummaryTable[ 'particleID' ][particleID] = particleID
            gSummaryTable[ 'v0_kps' ][particleID] = gv0_kps[ giPar==particleID ]
          #gSummaryTable[ '' ][particleID] = 

        






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

    pyplot.show()





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

    giEjecta, gtBegin_day, gtEnd_day = np.loadtxt( captureSummaryFile, unpack=True )
    
    h.plot( np.log10(gtEnd_day-gtBegin_day), xlabel='$\log_{10}$( capture duration / days )', nbins=20, xrange=(-1,4), title='crater '+str(icrater) )
        
    # capture start and stop times for each ejecta
    pyplot.figure( figsize=size, dpi=dpi )
    pyplot.ylim( (0,80) )
    pyplot.xscale( 'log' )
    for j in np.arange( len(giEjecta) ):
        x = [ gtBegin_day[j], gtEnd_day[j] ]
        y = [    giEjecta[j],  giEjecta[j] ]
        pyplot.plot( x, y, 'b', linestyle='-' )
    pyplot.xlabel( 'days from impact' )
    pyplot.ylabel( 'ejecta particle ID' )
    pyplot.title( 'crater ' + str(icrater) )
    pyplot.show()
    




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

