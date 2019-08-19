#!/bin/env/python

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

def main():
    df = pd.read_hdf('output.h5','df')
    #events = df['Vals']

    invert = True
    xmin=450
    xmax=550
    conv = 1.0
 
    signals = []

    df['Baseline'] = 0
    df['Integral'] = 0
    df['Mean'] = 0
    df['Max'] = 0
    #for i, event in enumerate(events):
    for i in range( df.shape[0] ):
    #for i in range( 10 ):
        event = df.loc[i,'Vals']

        #df['Itg'] = integrate(event)
        baseline = GetBaseLine(event, 200)
        #df['Bsl'] = baseline
        event -= baseline
        if invert: event *= -1

        df.loc[i,'Baseline'] = baseline
        df.loc[i,'Integral'] = integrate(event,xmin,xmax,conv)
        df.loc[i,'Mean']  = mean(event,xmin,xmax,conv)
        df.loc[i,'Max']  = GetMax(event)
        #df.loc[i,'Vals'] = event

    print(df)

    event = df.loc[0,'Vals']
    xdata = np.arange(xmin, (xmax if xmax < event.size else event.size) ) * conv
    event_r = event[xmin:xmax]
    popt, pcov = curve_fit(Gaus, xdata, event_r, p0=(1000.,500.,50.))   
    print (popt)
    print (pcov)

    #baseline = df.loc[0,'Baseline']
    popt_WithC, pcov_WithC = curve_fit(GausWithConstant, xdata, event_r, p0=(popt[0],popt[1],popt[2],0.))   
    print (popt_WithC)
    print (pcov_WithC)

    fig, axes = plt.subplots(2,2) 
    #for i in (0,100,200,300):
    for i in range(5):
        axes[0,0].plot( range(0, df.loc[i,'Vals'].size),df.loc[i,'Vals'] )

    #axes[1].hist( df['Itg'], bins=50 )
    axes[0,1].hist( df['Mean'], bins=50 )

    axes[1,0].plot( xdata, event_r, 'ko' )
    x_plot = np.linspace(xmin,xmax,200)
    axes[1,0].plot( x_plot, Gaus( x_plot, A=popt[0], mean=popt[1], sigma=popt[2] ), 'r' )
    axes[1,0].plot( x_plot, GausWithConstant( x_plot, A=popt_WithC[0], mean=popt_WithC[1], sigma=popt_WithC[2], C=popt_WithC[3] ), 'b' )

    plt.show()

def integrate(event, xmin, xmax, conv=1.0):
    #x = np.arange(xmin, (xmax if xmax < event.size else event.size) ) * conv
    #x = np.arange(0, event.size) * conv
    event_r = event[xmin:xmax]
    return np.sum( event_r )

def mean(event, xmin, xmax, conv=1.0):
    x = np.arange(xmin, (xmax if xmax < event.size else event.size) ) * conv
    event_r = event[xmin:xmax]
    mean = np.sum( event_r * x ) / np.sum( event_r )
    #print( x )
    #print( event_r )
    #print( np.sum( event_r ), np.sum( event_r * x ), event_r * x )
    return mean

#def Gaus(x, A, mean, sigma, C):
def Gaus(x, A, mean, sigma):
   #return A * np.exp( -0.5*( (x - mean)/sigma )**2 ) + C
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

