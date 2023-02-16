#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jul 22 17:14:40 2021

@author: rjedicke
"""

import matplotlib as mpl
import matplotlib.pyplot as pyplot

    
    


#-----------------------------------------------------------------------------------------------
# used to be plot2d
def plot2d( x, y, xlabel='', xrange=(), ylabel='', title='', yrange=(), markerstyle='.', 
           markersize=1, size=(8,8), dpi=200, bShow=False, xfunc=(), yfunc=(),
           logx=False, logy=False, saveto='' ):
    
    pyplot.figure( figsize=size, dpi=dpi )
    
    pyplot.plot( x, y, markerstyle, markersize=markersize )
    
    pyplot.xlabel( xlabel )    
    pyplot.ylabel( ylabel )
    pyplot.title(  title )
    
    if ( logy   ):  pyplot.yscale( "log" )
    if ( logx   ):  pyplot.xscale( "log" )
    
    if ( xrange ):  pyplot.xlim( xrange )
    if ( yrange ):  pyplot.ylim( yrange )
    
    # plot an additional function if defined
    if ( len(xfunc) != 0 ): pyplot.plot( xfunc, yfunc, 'g-' )
    
    if ( saveto ):  pyplot.savefig( saveto, bbox_inches='tight' )
    
    if ( bShow ):  pyplot.show()
    



#-----------------------------------------------------------------------------------------------
# used to be plot3d
def plot3d( x, y, z, datalabel='data', xlabel='', xrange=(), ylabel='', yrange=(), zlabel='', zrange=() ):
      
    mpl.rcParams['legend.fontsize'] = 10
    
    fig = pyplot.figure( figsize=(5,5), dpi=200 )
    
    ax = fig.gca( projection='3d' )
    
    ax.plot( x, y, z, label=datalabel )

    ax.legend()
    
    pyplot.show()




#-----------------------------------------------------------------------------------------------
def plot2horiz( x1, y1, x2, y2, x1label='', x1range=(), y1label='', y1range=(),
                                x2label='', x2range=(), y2label='', y2range=()):
    
    plot2d( x1, y1, x1label, x1range, y1label, y1range, size=(8,4), dpi=200, bShow=False )
    plot2d( x2, y2, x2label, x2range, y2label, y2range, size=(8,4), dpi=200, bShow=True  )
