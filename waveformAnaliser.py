#!/bin/env/python

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

def main():
    df = pd.read_hdf('output.h5','df')
    #events = df['Vals']

    invert = True

    signals = []

    df['Bsl'] = 0
    df['Itg'] = 0
    df['Pk'] = 0
    #for i, event in enumerate(events):
    for i in range( df.shape[0] ):
        event = df.loc[i,'Vals']

        #df['Itg'] = integrate(event)
        baseline = GetBaseLine(event, 200)
        #df['Bsl'] = baseline
        event -= baseline
        if invert: event *= -1

        df.loc[i,'Bsl'] = baseline
        df.loc[i,'Itg'] = integrate(event)
        df.loc[i,'Pk']  = GetMax(event)
        #df.loc[i,'Vals'] = event

    print(df)

    fig, axes = plt.subplots(1,2) 
    for i in range(5):
        axes[0].plot( range(0, df.loc[i,'Vals'].size),df.loc[i,'Vals'] )

    axes[1].hist( df['Itg'], bins=50 )

    plt.show()

def integrate(event):
    return np.sum(event)

def GetBaseLine(event, _range):
    baseline = 0
    for i in range(0, _range - 1):
        baseline += event[i]

    return baseline / _range

def GetMax(event):
    return np.nanmax(event)


if __name__ == '__main__':
    main()

