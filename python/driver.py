import lunarMinimoon as mm

# mm.loadCapture()

# mm.loadCaptureSummary()

# mm.loadInitialConditions()

# mm.loadCollisions()
# mm.loadEscapes() 

import matplotlib.pyplot as plt
import numpy as np
import PLOT_UTILITIES as plot
import HUTILITIES as hist


# if True: plt.show()

mm.createSummaryTable()

print(mm.gSummaryTable[2])
print()

# for i in mm.gSummaryTable: print(i)
plt.close('all')


unique_v0_kps , fraction_per_v0 , param_f = mm.fractionCapturedVSVelocity(show=True)


print(np.min(unique_v0_kps))
print('Fraction captured (boltzmann dist) parameters a, x_0' , param_f)

_ , avg_lifetime_per_v0 , param_l = mm.captureLifetimeVSVelocity()
print('capture lifetime (Gaussian dist) parameters mean, std , scale , t' , param_l)

# if True: plt.show()

impactor_diameter_km = np.logspace( -3 , 0 , 101 , base=10)
# impactor_diameter_km = np.arange( 0 , 5 , .01)
# cumulative_impact_rate = mm.cumulativeLunarImpactRate_Myr(impactor_diameter_km)
# plot.plot2d ( impactor_diameter_km , cumulative_impact_rate , title='cumulative impact rate vs impactor diameter' , xlabel='Diameter [km]' , ylabel='Cumulative impact rate per Myr' , logy=True )


# ejecta_diameter_km = np.arange(0.0001 , .1 , .001)
ejecta_diameter_km = np.logspace( -3 , 0 , 101 , base=10)
impactor_speed_kps = np.logspace( 0 , 2 , 101 , base=10)

# N = mm.nEjectaDensity ( 1 ,ejecta_diameter_km ,  )

# plot.plot2d ( impactor_speed_kps , N , xlabel='impactor_speed_kps [km]' , ylabel=' N (D, d , v)' , title='N(1km impactor, 1m ejecta, v)' , logx=True , logy=True )

ratio_ejectaD_impactorD = np.arange(0 , .1 , .0001)
ejecta_velocity_kps = mm.velocity_kps_vs_impactorD_ejectaD ( ratio_ejectaD_impactorD )
plot.plot2d ( ratio_ejectaD_impactorD , ejecta_velocity_kps , xlabel='D_ejecta/D_impactor' , ylabel='v ejecta [km/s]' , title='V ejecta vs D_ej/D_imp' , logx=True , logy=True )
# if True: plt.show()

# fig_VvsD , ax_VvsD = plt.subplots ( )
# for i in np.arange( .1 , 1.2 , .2 ) :
# # for i in range(0,len (impactor_diameter_km), 10): 
# 	# D_impactor = impactor_diameter_km[i]
# 	D_impactor = i
# 	ejecta_vel_kps = mm.velocity_kps_vs_impactorD_ejectaD ( ejecta_diameter_km/D_impactor )
# 	if np.any(ejecta_vel_kps > 2.4):
# 		ax_VvsD.scatter ( ejecta_diameter_km , ejecta_vel_kps , label=f'D_impactor={D_impactor} km' )
# ax_VvsD.axhline ( 2.4 , label='Lunar escape velocity' )
# ax_VvsD.set_xlabel ( 'D_ejecta[km]' )
# ax_VvsD.set_ylabel ( 'v_ejecta[km/s]' )
# ax_VvsD.set_title  ( 'v_ejecta vs D_ejecta' )
# ax_VvsD.legend     ()
# ax_VvsD.set_xscale ('log')
# ax_VvsD.set_yscale ('log')

minimoons = []
for ii in range( len( ejecta_diameter_km ) - 1): 
	minimoons.append(mm.minimoons ( ejecta_diameter_km[ii] )[0] * (ejecta_diameter_km[ii+1] - ejecta_diameter_km[ii])  )
print(minimoons)



# # plot.plot2d ( ejecta_diameter_km[:-1] , minimoons , xlabel='ejecta diameter [km]' , ylabel='n' , title=' n minimoons vs diameter_km' , logx=True , logy=True )
hist.plotAsHistogram ( ejecta_diameter_km[:] , minimoons , xlabel='ejecta diameter [km]' , ylabel='n' , title='incremental n minimoons vs diameter_km' , logx=True , logy=True )

plt.show()