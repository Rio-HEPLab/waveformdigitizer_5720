#!/bin/env/python

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import ipdb

XMIN = 450
XMAX = 550

def main():
    #df = pd.read_hdf('output.h5','df')
    df = pd.read_hdf('output-wave1-SingleTrig-V2-10k.h5','df')
    #df = pd.read_hdf('output-wave1-DoubleTrig.h5','df')

    invert = True
    xmin=400
    xmax=600
    conv = 1.0

    #ipdb.set_trace()

    signals = []

    df['Baseline'] = 0
    df['Integral'] = 0
    df['Mean'] = 0
    df['Max'] = 0
    for i in range( df.shape[0] ):
        event = df.loc[i,'Vals']

        baseline = GetBaseLine(event, 200)
        event -= baseline
        if invert: event *= -1

        #xmin, xmax = findPulseWidth(event, 1/3)

        df.loc[i,'Baseline'] = baseline
        df.loc[i,'Integral'] = integrate(event,xmin,xmax,conv)
        df.loc[i,'Mean']  = mean(event,xmin,xmax,conv)
        df.loc[i,'Max']  = GetMax(event)

    event = df.loc[0,'Vals']
    xmin, xmax = findPulseWidth(event, 1/3)
    xdata = np.arange(xmin, (xmax if xmax < event.size else event.size) ) * conv
    event_r = event[xmin:xmax]
    popt, pcov = curve_fit(Gaus, xdata, event_r, p0=(1000.,500.,50.))
    print (popt)
    print (pcov)

    popt_WithC, pcov_WithC = curve_fit(GausWithConstant, xdata, event_r, p0=(popt[0],popt[1],popt[2],0.))
    print (popt_WithC)
    print (pcov_WithC)

    fig, axes = plt.subplots(2,2)
    for i in range(5):
        axes[0,0].plot( range(0, df.loc[i,'Vals'].size),df.loc[i,'Vals'] )

    axes[0,1].hist( df['Mean'], bins=50 )

    axes[1,0].plot( xdata, event_r, 'ko' )
    x_plot = np.linspace(xmin,xmax,200)
    axes[1,0].plot( x_plot, Gaus( x_plot, A=popt[0], mean=popt[1], sigma=popt[2] ), 'r' )
    axes[1,0].plot( x_plot, GausWithConstant( x_plot, A=popt_WithC[0], mean=popt_WithC[1], sigma=popt_WithC[2], C=popt_WithC[3] ), 'b' )

    axes[1,1].hist( df['Integral'], bins=200, range=(0,10000) )

    plt.show()

def findPulseWidth(event, threshold):
    if (threshold > 1):
        return 0

    center = int(np.argmax(event))
    thresholdvalue = threshold * np.max(event)
    enumerateEvent = list(zip(range(event[1:center].size), event[1:center]))
    firstThreshold = bisearchLeft( enumerateEvent, thresholdvalue)
    if firstThreshold == None:
        firstThreshold = XMIN

    enumerateEvent = list(zip(range(event[center : event.size].size), event[center : event.size]))
    secondThreshold = bisearchRight(enumerateEvent, thresholdvalue) + center
    if secondThreshold == None:
        secondThreshold = XMAX

    return (firstThreshold - 10, secondThreshold + 10)

def bisearchLeft(array, value):
    x = int(len(array) / 2)

    if array[x][1] > value and array[x - 1][1] < value:
        return array[x][0]
    elif array[x][1] > value :
        return bisearchLeft(array[0 : x], value)
    elif array[x][1] < value :
        return bisearchLeft(array[x : len(array)], value)

def bisearchRight(array, value):
    x = int(len(array) / 2)

    if array[x][1] < value and array[x - 1][1] > value:
        return array[x][0]
    elif array[x][1] > value :
        return bisearchRight(array[x : len(array)], value)
    elif array[x][1] < value :
        return bisearchRight(array[0 : x], value)

def integrate(event, xmin, xmax, conv=1.0):
    #x = np.arange(xmin, (xmax if xmax < event.size else event.size) ) * conv
    #x = np.arange(0, event.size) * conv
    event_r = event[xmin:xmax]
    return np.sum( event_r )

def mean(event, xmin, xmax, conv=1.0):
    x = np.arange(xmin, (xmax if xmax < event.size else event.size) ) * conv
    event_r = event[xmin:xmax]
    mean = np.sum( event_r * x ) / np.sum( event_r )
    return mean

def Gaus(x, A, mean, sigma):
   return A * np.exp( -0.5*( (x - mean)/sigma )**2 )

def GausWithConstant(x, A, mean, sigma, C):
   return A * np.exp( -0.5*( (x - mean)/sigma )**2 ) + C

def GetBaseLine(event, _range):
    baseline = 0
    for i in range(0, _range - 1):
        baseline += event[i]

    return baseline / _range

def GetMax(event):
    return np.nanmax(event)


if __name__ == '__main__':
    main()
