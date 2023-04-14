import numpy as np
import lunarMinimoon_noPlot as mm

ejecta_diameter_km = np.logspace( -3 , -1 , 100 , base=10)

mm.createSummaryTable()

minimoons = []
for ii in range( len( ejecta_diameter_km ) - 1): 
	print(ii)
	minimoons.append(mm.minimoons ( ejecta_diameter_km[ii] )[0] * (ejecta_diameter_km[ii+1] - ejecta_diameter_km[ii])  )
print(minimoons)


output_file_name = '../data/outputs/minimoons.dat'
np.savetxt ( output_file_name , np.hstack ( ejecta_diameter_km , minimoons) )
