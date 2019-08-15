#!/bin/env/python

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

def main():
    df = pd.read_hdf('output.h5','df')
    events = df['Vals']

    signals = []

    for i, event in enumerate(events):
        df['Itg'] = integrate(event)
        baseline = GetBaseLine(event, 200)
        df['Bsl'] = baseline
        event -= baseline
        event *= -1
        df['Vals'][i] = event
        df['Pk'] = GetMax(event)

    print(df)
    plt.plot(range(0, df['Vals'][0].size),df['Vals'][0])
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
