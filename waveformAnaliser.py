#!/bin/env python

import argparse
import numpy as np
import pandas as pd
import h5py
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from scipy.integrate import simps
from scipy.special import erfc

import argparse
import ipdb

XMIN = 450
XMAX = 550

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
    secondThreshold = bisearchRight(enumerateEvent, thresholdvalue)
    if secondThreshold == None:
        secondThreshold = XMAX
    else:
        secondThreshold += center

    #assimetrico por causa do formato do pico ter maior tempo de queda que de subida
    return (firstThreshold - 10, secondThreshold + 150)

def bisearchLeft(array, value, it = 0):
    it += 1
    if (len(array) <= 0 or array == None or it > 100):
        return None

    x = int(len(array) / 2)

    if array[x][1] > value and array[x - 1][1] < value:
        return array[x][0]
    elif array[x][1] > value :
        return bisearchLeft(array[0 : x], value, it)
    elif array[x][1] < value :
        return bisearchLeft(array[x : len(array)], value, it)

def bisearchRight(array, value, it = 0):
    it += 1
    if (len(array) <= 0 or array == None or it > 100):
        return 0

    x = int(len(array) / 2)

    if array[x][1] < value and array[x - 1][1] > value:
        return array[x][0]
    elif array[x][1] > value :
        return bisearchRight(array[x : len(array)], value, it)
    elif array[x][1] < value :
        return bisearchRight(array[0 : x], value, it)

def integrate(event, xmin, xmax, conv=1.0):
    if (xmin <= 0 or xmin == None or xmax <= 0 or xmax == None or xmax - xmin < 3):
        return 0

    event_r = event[xmin:xmax]
    xdata = np.arange(xmin, (xmax if xmax < event.size else event.size) ) * conv
    return simps(event_r, xdata)

def mean(event, xmin, xmax, conv=1.0):
    xmin_corr = ( xmin if xmin >= 0 else 0 )
    xmax_corr = ( xmax if xmax < event.size else event.size )
    x = np.arange( xmin_corr, xmax_corr ) * conv
    event_r = event[xmin_corr:xmax_corr]
    mean = np.sum( event_r * x ) / np.sum( event_r )
    return mean

def GetBaseLine(event, _range):
    baseline = 0
    for i in range(0, _range - 1):
        baseline += event[i]

    return baseline / _range

def baseline(vals, n_low, n_high):
    vals_low = vals[:n_low]
    vals_high = vals[(vals.size - n_high):-1]
    vals_mean = ( np.sum(vals_low) + np.sum(vals_high) ) / ( vals_low.size + vals_high.size )

    return vals_mean

def GetMax(event):
    return np.nanmax(event)

def model_exp_gaus_res(x, *pars):
    C = pars[0]
    t0 = pars[1]
    sigma = pars[2]
    A = pars[3]
    tau = pars[4]
    
    val = A * ( 1./(2*tau) ) * np.exp( (sigma**2)/(2*(tau**2)) - (x - t0)/tau ) * erfc( sigma/(np.sqrt(2.)*tau) - (x - t0)/(np.sqrt(2.)*sigma) ) 
    val += C
    
    return val

class ModelExpRCRes:
    #def __init__(self, sigma=None):
    #    self.sigma = sigma
    def __call__(self, x, *pars):
        C = pars[0]
        t0 = pars[1]
        sigma = pars[2]
        A = pars[3]
        tau_0 = pars[4]
        a_1 = pars[5]
        tau_1 = pars[6]
    
        a   = np.array( ( (1. - a_1), a_1) )
        tau = np.array( (tau_0, tau_1) )
        
        val = ( ( a[0]/(2*tau[0]) ) * np.exp( (sigma**2)/(2*(tau[0]**2)) - (x - t0)/tau[0] ) * erfc( sigma/(np.sqrt(2.)*tau[0]) - (x - t0)/(np.sqrt(2.)*sigma) ) + 
            ( a[1]/(2*tau[1]) ) * np.exp( (sigma**2)/(2*(tau[1]**2)) - (x - t0)/tau[1] ) * erfc( sigma/(np.sqrt(2.)*tau[1]) - (x - t0)/(np.sqrt(2.)*sigma) ) )                
             
        val *= A
        val += C
    
        return val

def main():
    parser = argparse.ArgumentParser(description = 'Programa que recebe waveforms e extrai suas informações')
    #parser.add_argument('-df', action = 'store', dest = 'waveforms', required = True, help = 'Arquivo waveform do pandas' )
    parser.add_argument('--h5', action = 'store', dest = 'h5_file', required = True, help = 'Arquivo waveforms.' )
    args = parser.parse_args()

    df_out = pd.DataFrame(columns=['Baseline','Integral','Mean','Max'])

    #df = pd.read_hdf(arguments.waveforms,'df')
    with h5py.File( args.h5_file, 'r') as f:

        dset = f['Vals']

        invert = True
        xmin=400
        xmax=600
        conv = 1.0
    
        #ipdb.set_trace()
    
        #df['Baseline'] = 0
        #df['Integral'] = 0
        #df['Mean'] = 0
        #df['Max'] = 0
        #for i in range( df.shape[0] ):
        for i_evt in range( dset.shape[0] ):
            #event = df.loc[i,'Vals']
            event = dset[i_evt].copy()
            print ( "Event {:d}.format( i_evt)" )
            print ( event )   
 
            baseline = GetBaseLine(event, 200)
            event -= baseline
            if invert: event *= -1
            print ( event )
    
            xmin, xmax = findPulseWidth(event, 1/3)
            print ( "X(min), X(max) = {:d}, {:d}".format( xmin, xmax ) )    

            #df.loc[i,'Baseline'] = baseline
            #df.loc[i,'Integral'] = integrate(event,xmin,xmax,conv)
            #df.loc[i,'Mean']  = mean(event,xmin,xmax,conv)
            #df.loc[i,'Max']  = GetMax(event)
            df_out = df_out.append(
                {'Baseline' : baseline, 
                 'Integral' : integrate(event,xmin,xmax,conv), 
                 'Mean' : mean(event,xmin,xmax,conv), 
                 'Max': GetMax(event)
                }, 
                ignore_index=True 
                ) 
    
        #fig, axes = plt.subplots(2)
    
        print(df['Integral'])
    
        #axes[0].hist( df['Mean'], bins=50 )
        plt.hist( df['Integral'], bins=200 )
    
        plt.show(block = True)

if __name__ == '__main__':
    main()
