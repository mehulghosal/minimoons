import lunarMinimoon as mm

# import sys
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
mm.captureLifetimeVSVelocity()


print(param_f)
# if True: plt.show()
if False:
	fig_A5 , ax_A5 = plt.subplots ( )
	ejecta_diameter_km = np.array ( [1e-3 , 1e-2 , 1e-1] )
	impactor_diameter_km = np.linspace ( 0 , 10 , 1000 )
	lunar_escape_kps = 2.38
	for ejectaD_km in ejecta_diameter_km:
		ratio_ejectaD_impactorD = ejectaD_km / impactor_diameter_km
		ax_A5.plot ( impactor_diameter_km , mm.velocity_kps_vs_impactorD_ejectaD (ratio_ejectaD_impactorD) , label=f'{ejectaD_km*1000} m ejecta' )

	ax_A5.axhline (lunar_escape_kps , label='lunar escape speed')
	ax_A5.set_xlabel('impactor diameter [km]')
	ax_A5.set_ylabel('ejecta speed [km/s] ')
	ax_A5.set_ylim  (0 , 10)
	ax_A5.set_xlim  (0 , 10)
	ax_A5.legend()



	print(np.min(unique_v0_kps))
	print('Fraction captured (boltzmann dist) parameters a, x_0' , param_f)


	# if True: plt.show()
	marchi_small_param , marchi_large_param , yue_param = mm.digitize_plots ( show=True , poly_order=10 )
	impactor_diameter_km = np.logspace( -3 , 1 , 101 , base=10)
	# impactor_diameter_km = np.linspace(.001 , 1.5 , 101)
	impactor_speed_kps = np.logspace( 0 , 2 , 101 , base=10)
	# impactor_speed_kps = 10

	interpolated_polynomials = mm.interpolate_polynomials ( marchi_small_param , marchi_large_param , impactor_diameter_km )

	from scipy.integrate import quad , dblquad , tplquad
	incremental_impact_rate = []
	integrated_cum_impact_rate = []
	for i in range( len(impactor_diameter_km) -1 ):

		marchi_poly = np.poly1d(interpolated_polynomials[i])

		integrated_cum_impact_rate.append( mm.returnIntegrated_impactRate      ( impactor_diameter_km[i  ] , marchi_poly , ) )
		incremental_impact_rate   .append( mm.return_diffIntegrated_impactRate ( impactor_diameter_km[i  ] , marchi_poly , ) * (impactor_diameter_km[i+1]-impactor_diameter_km[i]) )


		# incremental_impact_rate.append( np.average([differential_impact_rate[i], differential_impact_rate[i+1]]) * ((impactor_diameter_km[i+1]-impactor_diameter_km[i])) )


	plot.plot2d (impactor_diameter_km[:-1] , integrated_cum_impact_rate , title='integrated cumulative impact rate vs impactor diameter' , xlabel='impactor diameter [km]' , ylabel='cumulative impact rate per Myr' , logy=True , logx=True , xrange=(0.1,10) , yrange=(1e-5 , 100) )


	hist.plotAsHistogram ( impactor_diameter_km , incremental_impact_rate,  title='inc impact rate vs impactor diameter' , xlabel='impactor diameter [km]' , ylabel='inc impact rate per Myr' , logy=True , logx=True, xrange=(0.1,10) , yrange=(1e-5 , 100))

	another_impactor_D_km = [.1  , 5 , 10 , 20 , 50 , 72 ]
	another_impactor_V_km = np.linspace (0,45,101)

	another_interpolation = mm.interpolate_polynomials ( marchi_small_param , marchi_large_param , another_impactor_D_km )
	fig_poly , ax_poly = plt.subplots (  )
	print(another_interpolation[0])

	for i in range ( len(another_impactor_D_km) ):
	# for i in range(1):
		marchi_poly = np.poly1d(another_interpolation[i])

		integrated_prob = quad(marchi_poly , 0 , 45)

		P = marchi_poly(another_impactor_V_km) / integrated_prob[0]
		cutoff = P>1e-3
		ax_poly.plot ( another_impactor_V_km[cutoff] , P[cutoff] , label=f'{another_impactor_D_km[i]} km impactor' )


	ax_poly.set_xlabel('impactor speed [km/s]')
	ax_poly.set_ylabel('impactor probability ')
	ax_poly.set_ylim  (0 , .07)
	ax_poly.legend()

# ejecta_diameter_km = np.arange(0.0001 , .1 , .001)


ejecta_diameter_km = np.logspace( -3 , 0 , 501 , base=10)
impactor_diameter_km = np.logspace( -3 , 1.7 , 501 , base=10)

# ejecta_diameter_km = np.linspace ( .001 , 1 , 501 )
# ejecta_diameter_km = np.linspace ( .001 , 1 , 501 )

impactor_speed_kps = np.array ( [5 , 10 , 20 , 30 , 40] )
fig_CraterDiameter_vs_impactorDiam , ax_CraterDiameter_vs_impactorDiam = plt.subplots ( )

for speed in impactor_speed_kps:
	ax_CraterDiameter_vs_impactorDiam. plot ( impactor_diameter_km , mm.craterDiameter_km(impactor_diameter_km  , speed) , label=f'{speed}km/s impactor')

ax_CraterDiameter_vs_impactorDiam.set_xlabel('impactor diameter [km]')
ax_CraterDiameter_vs_impactorDiam.set_ylabel('crater diameter [km] ')

ax_CraterDiameter_vs_impactorDiam.legend()


impactor_diameter_km = np.array ( [.1 , .2 , .5 , 1 , 2 , 5 , 10 , 20 ] )
# impactor_diameter_km = [.1]
fig_A4_left , ax_A4_left = plt.subplots ( )
for impactorD in impactor_diameter_km : 
	ax_A4_left . plot ( ejecta_diameter_km , mm.nEjectaDensity_ ( impactorD , ejecta_diameter_km , 12 ) , label=f'{impactorD} km impactor' )
	# plt.show()
ax_A4_left.set_xlabel('ejecta diameter [km]')
ax_A4_left.set_ylabel('cumulative N ejecta > d ')

ax_A4_left.set_ylim  (1 , 1e11)
ax_A4_left.set_xlim  (1e-3 , 1)
ax_A4_left.set_xscale('log')
ax_A4_left.set_yscale('log')
ax_A4_left.legend()

impactor_speed_kps   = np.array ( [ 8 , 12 , 20 , 35 ])
fig_A4_right , ax_A4_right = plt.subplots ( )
for impactorV in impactor_speed_kps : 
	ax_A4_right . plot ( ejecta_diameter_km , mm.nEjectaDensity_ ( 10 , ejecta_diameter_km , impactorV ) , label=f'{impactorV} kps impactor' )
	# plt.show()
ax_A4_right.set_xlabel('ejecta diameter [km]')
ax_A4_right.set_ylabel('cumulative N ejecta > d ')
ax_A4_right.set_ylim  (1 , 1e11)
ax_A4_right.set_xlim  (1e-3 , 1)
ax_A4_right.set_xscale('log')
ax_A4_right.set_yscale('log')
ax_A4_right.legend()


print()

plt.show()



