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
marchi_small_param , marchi_large_param , yue_param = mm.digitize_plots ( show=True )
impactor_diameter_km = np.logspace( -3 , 0 , 101 , base=10)
impactor_speed_kps = np.logspace( 0 , 2 , 101 , base=10)

# impactor_diameter_km = np.arange( 0 , 5 , .01)
cumulative_impact_rate = mm.cumulativeLunarImpactRate_Myr(1 , impactor_speed_kps , marchi_small_param)
plot.plot2d ( impactor_speed_kps , cumulative_impact_rate , title='differential impact rate vs impactor velocity of impactors smaller than 1km' , xlabel='impactor speed [km/s]' , ylabel='differential impact rate per Myr' , logy=True )

# if True: plt.show()
# ejecta_diameter_km = np.arange(0.0001 , .1 , .001)
ejecta_diameter_km = np.logspace( -3 , -1 , 50 , base=10)


# ratio_ejectaD_impactorD = ejecta_diameter_km/impactor_diameter_km
# ejecta_velocity_kps = mm.velocity_kps_vs_impactorD_ejectaD ( ratio_ejectaD_impactorD )

# N_1 = mm.nEjectaDensity_ ( 1 ,ejecta_diameter_km , impactor_speed_kps )
# N_2 = mm.nEjectaDensity_ ( impactor_diameter_km ,.001 , impactor_speed_kps )
# N_3 = mm.nEjectaDensity ( impactor_diameter_km ,ejec/ta_diameter_km , 12 )
# plot.plot2d ( impactor_speed_kps , N_1 , xlabel='impactor_speed_kps [km]' , ylabel=' N (1km, d , v)' , title='N(1km impactor, d, v)' , logx=True , logy=True )
# plot.plot2d ( ejecta_diameter_km , N_1 , xlabel='ejecta_diameter_km [km]' , ylabel=' N (1km, d , v)' , title='N(1km impactor, d, v)' , logx=True , logy=True )

# plot.plot2d ( impactor_diameter_km , N_2 , xlabel='impactor_diameter_km [km]' , ylabel=' N (D, 1m , v)' , title=' N (D, 1m , v)' , logx=True , logy=True )
# plot.plot2d ( impactor_speed_kps , N_2 , xlabel='impactor_speed_kps [km]' , ylabel=' N (D, 1m , v)' , title=' N (D, 1m , v)' , logx=True , logy=True )

# hist.plotAsHistogram ( impactor_speed_kps , N_1[:-1] , xlabel='impactor_speed_kps [kms]' , ylabel=' N (1km, d , v)' , title='N(1km impactor, d, v)' , logx=True , logy=True )
# hist.plotAsHistogram ( ejecta_diameter_km , N_1[:-1] , xlabel='ejecta_diameter_km [km]' , ylabel=' N (1km, d , v)' , title='N(1km impactor, d, v)' , logx=True , logy=True )

# hist.plotAsHistogram ( impactor_diameter_km , N_2[:-1] , xlabel='impactor_diameter_km [km]' , ylabel=' N (D, 1m , v)' , title=' N (D, 1m , v)' , logx=True , logy=True )
# hist.plotAsHistogram ( impactor_speed_kps , N_2[:-1] , xlabel='impactor_speed_kps [kms]' , ylabel=' N (D, 1m , v)' , title=' N (D, 1m , v)' , logx=True , logy=True )


# plot.plot2d ( ratio_ejectaD_impactorD , ejecta_velocity_kps , xlabel='D_ejecta/D_impactor' , ylabel='v ejecta [km/s]' , title='V ejecta vs D_ej/D_imp' , logx=True , logy=True )
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
	print(ii)
	minimoons.append(mm.minimoons ( ejecta_diameter_km[ii] )[0] * (ejecta_diameter_km[ii+1] - ejecta_diameter_km[ii])  )
print(minimoons)


output_file_name = '../data/outputs/minimoons.dat'
np.savetxt ( output_file_name , np.hstack ( ejecta_diameter_km , minimoons) )

# # plot.plot2d ( ejecta_diameter_km[:-1] , minimoons , xlabel='ejecta diameter [km]' , ylabel='n' , title=' n minimoons vs diameter_km' , logx=True , logy=True )
hist.plotAsHistogram ( ejecta_diameter_km[:] , minimoons , xlabel='ejecta diameter [km]' , ylabel='n' , title='incremental n minimoons vs diameter_km' , logx=True , logy=True )

plt.show()