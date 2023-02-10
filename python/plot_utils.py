import matplotlib.pyplot as plt
from scipy.ndimage.interpolation import shift
import numpy as np
from mpl_toolkits.axes_grid1 import make_axes_locatable

def plot ( x , y , xlabel='' , xlim=() , ylim=() , ylabel=() , markerstyle='.' , markersize=1 , 
		   size=(6,6) , dpi=100 , show=False , logx=False , logy=False , alpha=1 , title='', 
		   color='black' , linewidth=0 , linestyle='-' , label='' ) : 

	
	fig, ax = plt.subplots( figsize=size , dpi=dpi )
	ax.plot ( x , y , markerstyle , markersize=markersize , linestyle=linestyle , 
		linewidth=linewidth , color=color , alpha=alpha , label=label )
	# ax.legend()
	ax.set_title  ( title  )
	ax.set_xlabel ( xlabel )
	ax.set_ylabel ( ylabel )
	if xlim: ax.set_xlim   ( xlim   )
	if ylim: ax.set_ylim   ( ylim   )


	if logx : ax.set_xscale ( 'log' )
	if logy : ax.set_yscale ( 'log' )

	if show: fig.show( )
	return fig, ax


def plot_existing ( fig , ax , x , y , show=False , markerstyle='.' , markersize=1 , alpha=1 ,
					linewidth=0 , linestyle='-' , color='black' , label='' , leg=False , zorder=2 ,
					markerfacecolor='black' , markeredgecolor='black' ):
	
	ax.plot ( x , y , markerstyle , markersize=markersize , linestyle=linestyle , 
		linewidth=linewidth , color=color , alpha=alpha , label=label , zorder=2, 
		markerfacecolor=markerfacecolor , markeredgecolor=markerfacecolor )
	# if leg : ax.legend()

	if show: fig.show()
	return fig, ax

def show ( fig , leg_loc=(0,0)): 
	fig.legend(bbox_to_anchor=leg_loc)
	plt.tight_layout()
	fig.show()
def show ( leg_loc=(0,0) ): 
	# plt.legend(bbox_to_anchor=leg_loc)
	plt.tight_layout()
	plt.show()

def hist2d ( x , y , title='' , xlabel='' , ylabel='number' , binsx=100 , binsy=100 , xlim=() , dopdf=False ,
	 	   logx=False , logy=False , size=(6,6), dpi=100 , show=False , ylim=() , norm=False , 
	 	   cmap='hot' , alpha=1 , label='' ) :

	fig, ax = plt.subplots ( figsize=size , dpi=dpi )
	h,xbin,ybin,img = ax.hist2d ( x , y , bins=[binsx,binsy] , cmap=cmap , range=[ xlim , ylim ]  )
	divider = make_axes_locatable(ax)
	cax = divider.append_axes('right', size='5%', pad=0.05)
	fig.colorbar ( img , cax=cax )
	if xlim : ax.set_xlim  ( xlim )
	if ylim : ax.set_ylim  ( ylim )
	if logx : ax.set_xscale( 'log' )
	if logy : ax.set_yscale( 'log' )
	ax.set_xlabel ( xlabel )
	ax.set_ylabel ( ylabel )
	ax.set_title  ( title  )

	if show: fig.show ( )
	
	return fig, ax


def hist ( x , title='' , xlabel='' , ylabel='number' , bins=100 , xlim=() , dopdf=False ,
	 	   logx=False , logy=False , size=(6,6), dpi=100 , show=False , ylim=() , norm=False , 
	 	   color='black' , ret_hist=False , cumulative=False, label='' ) :

	nPerBin, binEdges, patches = 0 , 0 , 0
	nparray = np.ndarray.flatten( x )

	if( norm ):
		# nTotal = np.sum( nparray )
		nTotal = len(nparray)
	else:
		nTotal = 1
	
	fig , ax = plt.subplots ( figsize=size , dpi=dpi )
	
	if xlim : 
		ax.set_xlim ( xlim )
		if logx : 
			xmin , xmax = xlim[0] , xlim[1]
			logbins = 10 ** np.linspace( np.log10(xmin), np.log10(xmax) )
			nPerBin, binEdges, patches = ax.hist( nparray/nTotal, bins=logbins, density=dopdf, histtype='step', range=xlim , color=color , cumulative=cumulative )
		else:
			nPerBin, binEdges, patches = ax.hist( nparray/nTotal, bins=bins,   density=dopdf, histtype='step', range=xlim  , color=color , cumulative=cumulative, label=label)
	else:
		nPerBin, binEdges, patches = ax.hist( nparray/nTotal, bins, density=dopdf, histtype='step' , color=color, cumulative=cumulative, label=label )
	binCenters = (shift(binEdges,-1)[:-1]+binEdges[:-1])/2  #should use np.roll NOT shift! here and elsewhere
	
	ax.set_xlabel ( xlabel )
	ax.set_ylabel ( ylabel )
	ax.set_title  ( title  )
	if ylim : ax.set_ylim   ( ylim   )
	if logx : ax.set_xscale ( 'log' )
	if logy : ax.set_yscale ( 'log' )
	# if not label == '' : ax.legend()
	if show : fig.show( )

	ret = fig, ax
	if ret_hist: ret = fig, ax , nPerBin , binEdges , patches

	return ret

def hist_existing ( fig , ax , x, show=False , label='', color='black', bins=100, ret_hist=False , cumulative=False): 
	nPerBin, binEdges, patches = ax.hist ( x , label=label , color=color, bins=bins, cumulative=cumulative)
	# ax.legend()
	if show : fig.show()
	ret = fig, ax
	if ret_hist: ret = fig, ax, nPerBin, binEdges, patches
	return ret