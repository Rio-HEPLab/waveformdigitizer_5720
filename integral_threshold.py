import numpy as np
import h5py
import matplotlib.pyplot as plt
import EventAnalyze as analyze

with h5py.File('output-wave-06-11-2019.h5','r') as f:

    integral_vals = []
    x_leading_func_vals = []
    
    i = 0
    for event in f['Vals']:
        print(i)
        i += 1

        result = analyze.fit_event(event)
        if result != None:
            integral_vals.append( result.integral[0] )
            x_leading_func_vals.append( result.threshold )


    fig, axes = plt.subplots(2, figsize=(20,10))
    axes[0].hist( integral_vals, bins=50 )
    axes[1].hist( x_leading_func_vals, bins=50)
    plt.show()
