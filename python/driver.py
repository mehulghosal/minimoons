import lunarMinimoon as mm

# mm.loadCapture()

# mm.loadCaptureSummary()

# mm.loadInitialConditions()

# mm.loadCollisions()
# mm.loadEscapes() 

import matplotlib.pyplot as plt
import numpy as np
import PLOT_UTILITIES as plot


plt.show()

mm.createSummaryTable()

print(mm.gSummaryTable[2])
print()

# for i in mm.gSummaryTable: print(i)
plt.close('all')


unique_v0_kps , fraction_per_v0 , param_f = mm.fractionCapturedVSVelocity(show=False)

_ , avg_lifetime_per_v0 , param_l = mm.captureLifetimeVSVelocity()


impactor_diameter_km = np.arange(0 , 5 , .01)
cumulative_impact_rate = mm.cumulativeLunarImpactRate_Myr(impactor_diameter_km)
plot.plot2d ( impactor_diameter_km , cumulative_impact_rate , title='cumulative impact rate vs impactor diameter' , xlabel='Diameter [km]' , ylabel='Cumulative impact rate per Myr' , logy=True )

# 
# N = mm.nEjectaDensity ( diameter_km , 1 , unique_v0_kps )
# print( N )

# plot.plot2d ( diameter_km , N , xlabel='impactor diameter [km]' , ylabel='N (d, D , v)' , title='N ejecta for ejecta smaller than 1km' )


ejecta_diameter_km = np.arange(0.0001 , .1 , .001)
N = mm.nEjectaDensity ( 1 , ejecta_diameter_km , unique_v0_kps )
# print( N )

plot.plot2d ( ejecta_diameter_km , np.log10(N) , xlabel='ejecta diameter [km]' , ylabel='log10 N (d, D , v)' , title='log10 N ejecta vs ejecta_diameter_km for 1km impactor' , logx=True )



# integral = mm.minimoons( ejecta_diameter_km , 1 )
minimoons = []
for ii in range( len( ejecta_diameter_km ) ): 
	minimoons.append(mm.minimoons ( ejecta_diameter_km[ii] ))
print(minimoons)

plot.plot2d ( ejecta_diameter_km , minimoons , xlabel='ejecta diameter [km]' , ylabel='n' , title=' n minimoons vs diameter_km' , logx=True , logy=True )

plt.show()