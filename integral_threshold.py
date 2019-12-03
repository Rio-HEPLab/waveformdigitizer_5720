import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import EventAnalyze as analyze

maxEvents = 10000
df = pd.read_hdf('output.h5','df')

integral_vals = []
x_leading_func_vals = []

for i in range( df.shape[0] ):
    if maxEvents >= 0 and i >= maxEvents: break

    event = df.loc[i,'Vals']
    
    result = analyze.fit_event(event)
    if result != None:
        integral_vals.append( result.integral[0] )
        x_leading_func_vals.append( result.threshold )


fig, axes = plt.subplots(2, figsize=(20,10))
axes[0].hist( integral_vals, bins=50 )
axes[1].hist( x_leading_func_vals, bins=50)
plt.show()
