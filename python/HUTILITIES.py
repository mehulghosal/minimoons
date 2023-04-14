#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jul 17 10:06:21 2020

@author: rjedicke
"""

import numpy as np
import matplotlib.pyplot as pyplot
import math
from scipy.optimize import curve_fit
from scipy.ndimage.interpolation import shift
from scipy.constants import pi
# from lmfit import Parameters
# from lmfit.models import GaussianModel

import utilities as util

global nPerBin
global nValues
global bins
global patches
global mode
global binCenters
global binEdges
global mean
global median
global rms
global sem
global labelx

root2pi = math.sqrt( 2 * pi )



# utilities for analyzing and plotting histograms

median  = 0.0
mean    = 0.0
mode    = 0.0
sem     = 0.0  # standard error of the mean
rms     = 0.0
nValues = 0   # number of values in the array

nPerBin    = np.zeros( 1 )  # entries in each bin
binCenters = np.zeros( 1 )





#-----------------------------------------------------------------------------------
# basic histogram of an array with basic statistics
def plot( npx, title='', xlabel='', ylabel='number', nbins=100, xrange=(), yrange=(), doNormalize=False,
         info=True, logx=False, logy=False, saveto="", 
          xsize=8, ysize=8, dpi=200, doReturn=False, doShow=False,
          bShowErrorBars=False, bPoissonErrors=False, bUseYErrorVectors=False, yErrorPos=[], yErrorNeg=[] ):
    
    global nPerBin
    global nValues
    global bins
    global patches
    global mode
    global binEdges
    global binCenters
    global labelx

    nparray = np.ndarray.flatten( npx ) # converts the dimensionality of the input array into a column array
    labelx = xlabel
    fig, ax = pyplot.subplots( figsize=(xsize, ysize), dpi=dpi )
    
    if( doNormalize ):
        nTotal = np.sum( nparray[ ( xrange[0] <= nparray ) * ( nparray < xrange[1] ) ] )
    else:
        nTotal = 1

    weights = np.ones_like( nparray ) / nTotal

    if( xrange ):
        pyplot.xlim( xrange )
        if ( logx ):
            xmin = xrange[0]
            xmax = xrange[1]
            logbins = 10 ** np.linspace( np.log10(xmin), np.log10(xmax) )
            nPerBin, binEdges, patches = ax.hist( nparray, bins=logbins, weights=weights, histtype='step', range=xrange )
        else:
            nPerBin, binEdges, patches = ax.hist( nparray, bins=nbins,   weights=weights, histtype='step', range=xrange )
    else:
        nPerBin, binEdges, patches = ax.hist( nparray, bins=nbins, weights=weights, histtype='step' )
    
    binCenters = (shift(binEdges,-1)[:-1]+binEdges[:-1])/2  #should use np.roll NOT shift! here and elsewhere

    if( info ):
        stats( nparray, xrange )
        textstr = '\n'.join((
            r'$N=%5.3e$'         % (nValues, ),
            r'$\bar{x}=%5.3e$'   % (mean, ),
            r'$\tilde{x}=%5.3e$' % (median, ),
            r'$\sigma=%5.3e$'    % (rms, )))
        ax.text(0.02, 0.98, textstr, transform=ax.transAxes, verticalalignment='top', fontsize=10 )
        
    if( bShowErrorBars ):
        if( bPoissonErrors ):  yError = np.sqrt( nPerBin );  yErrorPos = yError;  yErrorNeg = yError
        pyplot.errorbar( binCenters, nPerBin, yerr=[ yErrorNeg, yErrorPos ], fmt='none' )
        
    modeBin = np.argmax( nPerBin )
    mode = ( binEdges[modeBin] + binEdges[modeBin+1] ) / 2

    pyplot.xlabel( xlabel ) 
    pyplot.ylabel( ylabel ) 

    if( title != '' ):  pyplot.title( title )
    
    if ( yrange ):  pyplot.ylim( yrange )
    if ( logy   ):  pyplot.yscale( "log" )
    if ( logx   ):  pyplot.xscale( "log" )
    if ( saveto ):  pyplot.savefig( saveto, bbox_inches='tight' )
    if ( doShow ):  pyplot.show()
    else:           pyplot.close()
    
    if ( doReturn ):  return nPerBin, binEdges, binCenters
    


    
#-----------------------------------------------------------------------------------
# basic array statistics
def stats( nparray, xrange=() ):
    global nPerBin
    global nValues
    global median
    global mean
    global rms
    global sem
    if ( xrange ):
        nparrayInRange = nparray[ (nparray>=xrange[0]) & (nparray<xrange[1]) ]
    else:
        nparrayInRange = nparray
    nValues =      len( nparrayInRange ) 
    median = np.median( nparrayInRange )
    mean   = np.mean(   nparrayInRange )
    rms    = np.std(    nparrayInRange ) 
    sem    = rms / math.sqrt( nValues )




#-----------------------------------------------------------------------------------
# print basic array statistics
def stats_print( nparray, xrange=(), szName="" ):
    print( "\n" )
    print( szName )
    stats( nparray, xrange )
    print( "n      = ", nValues )
    print( "median = ", median )
    print( "mean   = ", mean )
    print( "sem    = ", sem )
    print( "std    = ", rms )
    print( "\n" )
    
    
    
    
#-----------------------------------------------------------------------------------
# basic histogram of an array with basic statistics
def plotCumulative( nparray, xlabel="", nBins = 100 ):
    pyplot.figure( figsize=(8, 8), dpi=80 )
    pyplot.hist( nparray, nBins, cumulative=True, histtype="step" )
    pyplot.xlabel( xlabel )
    pyplot.ylabel( "number" )
    pyplot.show()
    stats( nparray )
    
    
    
    
    
#-----------------------------------------------------------------------------------
# basic histogram of an array with basic statistics
def plotCumulativeNormalized( nparray, xlabel="", nBins = 100 ):
    pyplot.figure( figsize=(8, 8), dpi=80 )
    pyplot.hist( nparray, nBins, cumulative=True, density=True, histtype="step" )
    pyplot.xlabel( xlabel )
    pyplot.ylabel( "number" )
    pyplot.show()
    stats( nparray )
 



#-----------------------------------------------------------------------------------
# basic histogram of an array with basic statistics
# def plot( npx, xlabel="", ylabel='number', nbins=100, xrange=(), info=True, logx=False, logy=False, saveto="", 
#          xsize=8, ysize=8, dpi=200, doReturn=False ):
#     global nPerBin
#     global nValues
#     global bins
#     global patches
#     global mode
#     global binEdges
#     global binCenters
#     global labelx
#     nparray = np.ndarray.flatten( npx ) # converts the dimensionality of the input array into a column array
#     labelx = xlabel
#     stats( nparray, xrange ) 
#     fig, ax = pyplot.subplots( figsize=(xsize, ysize), dpi=dpi )
#     if ( xrange ):
#         pyplot.xlim( xrange )
#         if ( logx ):
#             xmin = xrange[0]
#             xmax = xrange[1]
#             logbins = 10 ** np.linspace( np.log10(xmin), np.log10(xmax) )
#             nPerBin, binEdges, patches = ax.hist( nparray, bins=logbins, histtype='step', range=xrange )
#         else:
#             nPerBin, binEdges, patches = ax.hist( nparray, bins=nbins,   histtype='step', range=xrange )
#     else:
#         nPerBin, binEdges, patches = ax.hist( nparray, nbins, histtype='step' )
#     modeBin = np.argmax( nPerBin )
#     mode = ( binEdges[modeBin] + binEdges[modeBin+1] ) / 2
#     #print( "mode = ", mode, "( bin", modeBin, ")" )
#     pyplot.xlabel( xlabel )
#     pyplot.ylabel( ylabel )
#     if ( info ):
#         textstr = '\n'.join((
#             r'$N=%5.3e$'         % (nValues, ),
#             r'$\bar{x}=%5.3e$'   % (mean, ),
#             r'$\tilde{x}=%5.3e$' % (median, ),
#             r'$\sigma=%5.3e$'    % (rms, )))
#         ax.text(0.02, 0.98, textstr, transform=ax.transAxes, verticalalignment='top', fontsize=10 )
#     if ( logy ):
#         pyplot.yscale( "log" )
#     if ( logx ):
#         pyplot.xscale( "log" )
#     if ( saveto ):
#         pyplot.savefig( saveto, bbox_inches='tight' )
#     pyplot.show()
    
#     binCenters = (shift(binEdges,-1)[:-1]+binEdges[:-1])/2
    
#     if ( doReturn ):
#         return nPerBin, binEdges, binCenters






#-----------------------------------------------------------------------------------
# basic weighted histogram of an array with basic statistics
def plot_weighted( npx, weights, xlabel="", ylabel='number', nbins=100, xrange=(), logx=False, yrange=(), logy=False, info=True, saveto="", xsize=8, ysize=8, dpi=200, doReturn=False, doShow=True ):
    global nPerBin
    global nValues
    global bins
    global patches
    global mode
    global binEdges
    global binCenters
    global labelx
    
    nparray = np.ndarray.flatten( npx ) # converts the dimensionality of the input array into a column array
    labelx = xlabel
    stats( nparray, xrange ) 
    fig, ax = pyplot.subplots( figsize=(xsize, ysize), dpi=dpi )
    
    if ( xrange ):
        pyplot.xlim( xrange )
        if ( logx ):
            xmin = xrange[0]
            xmax = xrange[1]
            logbins = 10 ** np.linspace( np.log10(xmin), np.log10(xmax) )
            nPerBin, binEdges, patches = ax.hist( nparray, weights=weights, bins=logbins, histtype='step', range=xrange )
        else:
            nPerBin, binEdges, patches = ax.hist( nparray, weights=weights, bins=nbins,   histtype='step', range=xrange )
    else:
        nPerBin,     binEdges, patches = ax.hist( nparray, weights=weights, bins=nbins,   histtype='step' )
        
    modeBin = np.argmax( nPerBin )
    mode = ( binEdges[modeBin] + binEdges[modeBin+1] ) / 2
    
    binCenters = (shift(binEdges,-1)[:-1]+binEdges[:-1])/2

    pyplot.xlabel( xlabel )
    pyplot.ylabel( ylabel )
    
    if ( info ):
        textstr = '\n'.join((
            r'$N=%5.3e$'         % (nValues, ),
            r'$\bar{x}=%5.3e$'   % (mean, ),
            r'$\tilde{x}=%5.3e$' % (median, ),
            r'$\sigma=%5.3e$'    % (rms, )))
        ax.text(0.02, 0.98, textstr, transform=ax.transAxes, verticalalignment='top', fontsize=10 )
        
    if ( logx   ):  pyplot.xscale( "log" )
    if ( yrange ):  pyplot.ylim( yrange )
    if ( logy   ):  pyplot.yscale( "log" )
    if ( saveto ):  pyplot.savefig( saveto, bbox_inches='tight' )
    if ( doShow ):  pyplot.show()
    else:           pyplot.close()
    
    if ( doReturn ):  return nPerBin, binEdges
       
    
    
    
    
#-----------------------------------------------------------------------------------
# plots  npn vs npx as if the data is a histogram instead of histogramming the array
# assumes that npx is the bin edges returned by np.hist
def plotAsHistogram( npBinEdges, npValues, title='', xrange=(), xlabel="", xticks=(), yrange=(), ylabel='value', logx=False, logy=False, doShow=False,
                    bShowErrorBars=False, bPoissonErrors=False, bUseYErrorVectors=False, yErrorNeg=[], yErrorPos=[],
                    plotFunc=True, xfunc=(), yfunc=(), figsize=(8,8), dpi=200 ):
    
    pyplot.figure( figsize=figsize, dpi=dpi )
    
    if ( xrange ): pyplot.xlim( xrange )
    if ( yrange ): pyplot.ylim( yrange )
        
    npValuesNoNaN = np.nan_to_num( npValues )
    
    pyplot.hist( npBinEdges[:-1], npBinEdges, weights=npValuesNoNaN, histtype="step" )
    
    if ( logy ): pyplot.yscale( "log" )
        
    if ( logx ): pyplot.xscale( "log" )
    
    binCenters = (shift(npBinEdges,-1)[:-1]+npBinEdges[:-1])/2
        
    if( bShowErrorBars ):
        if( bPoissonErrors ):  yError = np.sqrt( npValues );  yErrorPos = yError;  yErrorNeg = yError         
        pyplot.errorbar( binCenters, npValues, yerr=[ yErrorNeg, yErrorPos ], fmt='none' )

    pyplot.xlabel( xlabel )
    pyplot.ylabel( ylabel )

    if( plotFunc ):  pyplot.plot( xfunc, yfunc, 'g-' )
    
    if( xticks ):  pyplot.xticks( xticks[0], xticks[1] )

    if( title != '' ):  pyplot.title( title )
    
    if ( doShow ):  pyplot.show()
    # else:           pyplot.close()
    
    
    
    
#-----------------------------------------------------------------------------------
# basic histogram of an array with basic statistics
def plotWithXLimits( nparray, xlabel="", xMin=0.0, xMax=1.0,  nBins=100 ):
    pyplot.figure( figsize=(8, 8), dpi=80 )
    pyplot.hist( nparray, nBins, range=(xMin,xMax), histtype="step" )
    pyplot.xlabel( xlabel )
    pyplot.ylabel( "number" )
    pyplot.show()
    stats( nparray )




#-----------------------------------------------------------------------------------
# fits the current histogram to the specified function
def fit_g( xrange=(), gNorm=-1, gMean=-1, gSigma=-1 ):
    global binEdges
    global binCenters
    global nPerBin
    global nValues
    global mean
    global rms
    global labelx
    
    if ( -1 == gNorm ):
        print('nValues, mean, rms = ', nValues, mean, rms )
        popt, pcov = curve_fit( gauss, binCenters, nPerBin, p0=[nValues,3.0,1.0] )
    else:
        print('nValues, mean, rms = ', gNorm, gMean, gSigma )
        popt, pcov = curve_fit( gauss, binCenters, nPerBin, p0=[gNorm,gMean,gSigma] )
    print( 'popt = ', popt )
    
    pyplot.figure( figsize=(8, 8), dpi=80 )
    if ( xrange ):
        pyplot.xlim( xrange )
    npValuesNoNaN = np.nan_to_num( nPerBin )
    pyplot.hist( binEdges[:-1], binEdges, weights=npValuesNoNaN, histtype="step" )
    pyplot.xlabel( labelx )
    pyplot.ylabel( 'number' )
    xmin = binEdges[0]
    xmax = binEdges[len(binEdges)-1]
    xstep = ( xmax - xmin ) / 1000
    x = np.arange( xmin, xmax, xstep )
    pyplot.plot( x, gauss(x,*popt), 'r' )
    pyplot.show()
    
    return popt, pcov
    

    
    
#-----------------------------------------------------------------------------------
def gaussgauss( x, norm1, mean1, sigma1, norm2, mean2, sigma2 ):
    return norm1 * np.exp( - 0.5 * ( (x-mean1) / sigma1 )**2 ) / sigma1 / root2pi \
        +  norm2 * np.exp( - 0.5 * ( (x-mean2) / sigma2 )**2 ) / sigma2 / root2pi





#-----------------------------------------------------------------------------------
# fits the current histogram to the specified function
def fit_gg( gNorm1, gMean1, gSigma1, gNorm2, gMean2, gSigma2, xrange=() ):
    global binEdges
    global binCenters
    global nPerBin
    global nValues
    global mean
    global rms
    global labelx
    
    popt, pcov = curve_fit( gaussgauss, binCenters, nPerBin, sigma=np.sqrt(nPerBin), absolute_sigma=True, p0=[gNorm1,gMean1,gSigma1,gNorm2,gMean2,gSigma2] )
    #popt=[gNorm1,gMean1,gSigma1,gNorm2,gMean2,gSigma2]
    
    pyplot.figure( figsize=(8, 8), dpi=80 )
    if ( xrange ):
        pyplot.xlim( xrange )
    npValuesNoNaN = np.nan_to_num( nPerBin )
    pyplot.hist( binEdges[:-1], binEdges, weights=npValuesNoNaN, histtype="step" )
    pyplot.xlabel( labelx )
    pyplot.ylabel( 'number' )
    xmin = binEdges[0]
    xmax = binEdges[len(binEdges)-1]
    pyplot.xlim( xmin, xmax )
    xstep = ( xmax - xmin ) / 1000
    x = np.arange( xmin, xmax, xstep )
    pyplot.plot( x, gaussgauss(x,*popt)*0.1, 'r' )
    pyplot.show()
    
    return popt, pcov





#-----------------------------------------------------------------------------------
# fits the current histogram to a double gaussian
def fitDoubleGaussian( xFitRange=(), \
                    g1='g1', g1center=0, g1sigma=1, g1amplitude=1, g1centerVARY=True, g1sigmaVARY=True, g1amplitudeVARY=True, \
                    g2='g2', g2center=0, g2sigma=1, g2amplitude=1, g2centerVARY=True, g2sigmaVARY=True, g2amplitudeVARY=True, \
                    showFit=True, showg1=False, showg2=False, showInit=False ):

    global binEdges
    global binCenters
    global nPerBin    
    global labelx
    
    if ( xFitRange ):
        data_x     = binCenters[ (binCenters >= xFitRange[0]) & (binCenters <= xFitRange[1]) ]
        data_y     =    nPerBin[ (binCenters >= xFitRange[0]) & (binCenters <= xFitRange[1]) ]
        data_edges =   binEdges[ (binEdges   >= xFitRange[0]) & (binEdges   <= xFitRange[1]) ]
    else:
        data_x     = binCenters
        data_y     = nPerBin
        data_edges = binEdges
    
    parameters = Parameters()
    
    gauss1 = GaussianModel( prefix=g1 )    
    parameters.update( gauss1.make_params( ))
    
    parameters[g1+'center'   ].set( value=g1center,    vary=g1centerVARY )
    parameters[g1+'sigma'    ].set( value=g1sigma,     vary=g1sigmaVARY )
    parameters[g1+'amplitude'].set( value=g1amplitude, vary=g1amplitudeVARY )
    
    gauss2 = GaussianModel( prefix=g2 )
    parameters.update( gauss2.make_params( ))
    
    parameters[g2+'center'   ].set( value=g2center,    vary=g2centerVARY )
    parameters[g2+'sigma'    ].set( value=g2sigma,     vary=g2sigmaVARY )
    parameters[g2+'amplitude'].set( value=g2amplitude, vary=g2amplitudeVARY )
    
    parameters_init = parameters
    
    function = gauss1 + gauss2
    
    #init = function.eval(        parameters, x=data_x )
    out  = function.fit( data_y, parameters, x=data_x )
    
    parameters_final = out.params
    
    print( out.fit_report( min_correl=0.5) )    
    
    pyplot.figure( figsize=(8, 8), dpi=80 )
    
    #pyplot.plot( data_x, data_y,       'b' )
    
    # plot the histogram
    pyplot.hist( binEdges[:-1], binEdges, weights=nPerBin, histtype="step" )
    pyplot.xlabel( labelx )
    pyplot.ylabel( 'number' )
    pyplot.xlim( binEdges[0], binEdges[len(binEdges)-1] )
    
    # plot the function(s) in the fitting range
    xmin = data_edges[0]
    xmax = data_edges[len(data_edges)-1]
    xstep = ( xmax - xmin ) / 1000
    xHighRes = np.arange( xmin, xmax, xstep )

    if ( showInit ):
        pyplot.plot( xHighRes, function.eval( parameters_init,  x=xHighRes ), 'k--', label='initial fit' ) # initial function
        
    if ( showFit  ):
        pyplot.plot( xHighRes, function.eval( parameters_final, x=xHighRes ), 'r-',  label='best fit' ) # final function
    
    comps = out.eval_components( x=xHighRes )
    if ( showg1 ):
        pyplot.plot( xHighRes, comps[g1], 'g--', label=g1 )   # first  gaussian only
        
    if ( showg2 ):
        pyplot.plot( xHighRes, comps[g2], 'm--', label=g2 )   # second gaussian only
    
    pyplot.legend( loc='best' )
  
    pyplot.show()




#-----------------------------------------------------------------------------------
def logarithmicalBins( xrange ):

    #create bin edges aligned with the log scale
    nBins = 9 * (math.log10(xrange[1]) - math.log10(xrange[0]) )
    nEdges = nBins + 1        
    binEdges = np.zeros( int(nEdges) )       
    print( 'nBins, nEdges = ', nBins, nEdges )
    iBin = 0        
    for jlog in np.arange( 1, 5 ):            
        fj = 10**jlog            
        for i in np.arange( 1, 10 ):            
            binEdges[iBin] = fj * i
            iBin += 1                
    binEdges[iBin] = 10**math.log10(xrange[1])

