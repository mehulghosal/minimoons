import lunarMinimoon as mm

# mm.loadCapture()

# mm.loadCaptureSummary()

# mm.loadInitialConditions()

# mm.loadCollisions()
# mm.loadEscapes() 

import matplotlib.pyplot as plt

plt.show()

mm.createSummaryTable()

print(mm.gSummaryTable[2])
print()

# for i in mm.gSummaryTable: print(i)
plt.close('all')


unique_v0_kps , fraction_per_v0 = mm.fractionCapturedVSVelocity(show=False)

mm.captureLifetimeVSVelocity()
