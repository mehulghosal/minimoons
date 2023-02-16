import lunarMinimoon as mm

# mm.loadCapture()

# mm.loadCaptureSummary()

# mm.loadInitialConditions()

# mm.loadCollisions()
# mm.loadEscapes() 

import matplotlib.pyplot as plt

plt.show()

mm.createSummaryTable()

print(mm.gSummaryTable[0])

for i in mm.gSummaryTable: print(i)