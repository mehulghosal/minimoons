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


plt.show()

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

impactor_diameter_km = np.arange(0 , 5 , .01)
cumulative_impact_rate = mm.cumulativeLunarImpactRate_Myr(impactor_diameter_km)
plot.plot2d ( impactor_diameter_km , cumulative_impact_rate , title='cumulative impact rate vs impactor diameter' , xlabel='Diameter [km]' , ylabel='Cumulative impact rate per Myr' , logy=True )


# ejecta_diameter_km = np.arange(0.0001 , .1 , .001)
ejecta_diameter_km = np.logspace( -3 , 0 , 101 , base=10)
N = mm.nEjectaDensity ( 1 , ejecta_diameter_km , unique_v0_kps[3] )
# # print( N )

plot.plot2d ( ejecta_diameter_km , N , xlabel='ejecta diameter [km]' , ylabel=' N (d, D , v)' , title='cumulative N ejecta vs ejecta_diameter_km for 1km impactor' , logx=True , logy=True )



minimoons = []
for ii in range( len( ejecta_diameter_km ) - 1): 
	minimoons.append(mm.minimoons ( ejecta_diameter_km[ii] )[0] * (ejecta_diameter_km[ii+1] - ejecta_diameter_km[ii])  )
print(minimoons)



# plot.plot2d ( ejecta_diameter_km[:-1] , minimoons , xlabel='ejecta diameter [km]' , ylabel='n' , title=' n minimoons vs diameter_km' , logx=True , logy=True )
hist.plotAsHistogram ( ejecta_diameter_km[:] , minimoons , xlabel='ejecta diameter [km]' , ylabel='n' , title='incremental n minimoons vs diameter_km' , logx=True , logy=True )

plt.show()