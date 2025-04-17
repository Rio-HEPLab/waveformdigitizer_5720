#!/usr/bin/env python

import argparse
import numpy as np
import pandas as pd
import h5py
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from scipy.integrate import simpson
from scipy.special import erfc

import argparse
#import ipdb

from fit_event import *

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
    return simpson(event_r, xdata)

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


def main():
    parser = argparse.ArgumentParser(description = 'Programa que recebe waveforms e extrai suas informações')
    #parser.add_argument('-df', action = 'store', dest = 'waveforms', required = True, help = 'Arquivo waveform do pandas' )
    parser.add_argument('--h5', '--h5_0', action = 'store', dest = 'h5_file_0', required = True, help = 'Arquivo waveforms.' )
    parser.add_argument('--h5_1', action = 'store', dest = 'h5_file_1', required = False, help = 'Arquivo waveforms.' )
    parser.add_argument('-n', '--events', dest = 'events', type = int, required = False, default = 10000, help = 'Numero de eventos' )
    args = parser.parse_args()

    fit_selected_events = False
    show_channel = 1

    #df = pd.read_hdf(arguments.waveforms,'df')

    #with h5py.File( args.h5_file, 'r') as f:

    run_double_channel = False 
    if args.h5_file_1:
        run_double_channel = True

    if run_double_channel:
        print ( "Running on channels 0/1" )
    else:
        print ( "Running on channel 0" )

    if not run_double_channel: show_channel = 0

    f0 = h5py.File( args.h5_file_0, 'r')
    f1 = None
    if run_double_channel: f1 = h5py.File( args.h5_file_1, 'r')
    
    #df_out = pd.DataFrame(columns=['Baseline','Integral','Mean','Max'])
    df_out = None
    if run_double_channel:
        df_out = pd.DataFrame( 
            columns=[ 'Baseline_ch0','Integral_ch0','Mean_ch0','Max_ch0','FitIntegral_ch0','BinnedLE_ch0','FitLE_ch0','FitChi2_ch0','FitSigma_ch0','FitTau0_ch0','FitTau1_ch0','FitA1_ch0',
                      'Baseline_ch1','Integral_ch1','Mean_ch1','Max_ch1','FitIntegral_ch1','BinnedLE_ch1','FitLE_ch1','FitChi2_ch1','FitSigma_ch1','FitTau0_ch1','FitTau1_ch1','FitA1_ch1' ] 
            )
    else: 
        df_out = pd.DataFrame( columns=['Baseline','Integral','Mean','Max','FitIntegral','BinnedLE','FitLE','FitChi2','FitSigma','FitTau0','FitTau1','FitA1'] )

    if f0: 
        #dset = f['Vals']
        dset0 = f0['Waveform']
        dset_metadata0 = f0['Metadata']
        print ( "Number of events: {:d}".format( dset0.shape[0] ) )
    
        dset1 = None
        dset_metadata1 = None
        if run_double_channel and f1:
            dset1 = f1['Waveform']
            dset_metadata1 = f1['Metadata']
            print ( "Number of events: {:d}".format( dset1.shape[0] ) )
        
        df0 = pd.DataFrame( columns=('Event','Channel','Waveform') )
        df0['Event']   = dset_metadata0[:,0]
        df0['Channel'] = dset_metadata0[:,1]
        # df0['Waveform'] = dset0
        for i in range( dset0.shape[0] ):
            df0[ 'Waveform' ].iloc[ i ] = dset0[ i ]
        # df0 = df0.set_index( 'Event' )

        df1 = None
        if run_double_channel and dset1:
            df1 = pd.DataFrame( columns=('Event','Channel','Waveform') )
            df1['Event']   = dset_metadata1[:,0]
            df1['Channel'] = dset_metadata1[:,1]
            # df1['Waveform'] = dset1
            for i in range( dset1.shape[0] ):
                df1[ 'Waveform' ].iloc[ i ] = dset1[ i ]
            # df1 = df1.set_index( 'Event' )

        df = None
        if run_double_channel:
            print ( "Merging datasets Ch0/1" )
            df = df0.merge( df1, left_on='Event', right_on='Event', suffixes=('_ch0','_ch1') )
        else:
            df = df0

        print ( df ) 

        invert = True
        xmin=400
        xmax=600
        conv = 4.0
 
        rows = 4
        cols = 4
        #event_numbers_rand = np.random.randint( 0, dset.shape[0], ( rows*cols ) )
        event_numbers_rand = np.random.randint( 0, df.shape[0], ( rows*cols ) )
        
        fig, axes = None, None
        if fit_selected_events:
            fig, axes = plt.subplots( rows, cols, figsize=(20,20) )   

        #ipdb.set_trace()

        #df['Baseline'] = 0
        #df['Integral'] = 0
        #df['Mean'] = 0
        #df['Max'] = 0
        #for i in range( df.shape[0] )
       
        #event_numbers = event_numbers_rand if fit_selected_events else range( dset.shape[0] )
        event_numbers = event_numbers_rand if fit_selected_events else range( df.shape[0] )
        i_row = 0
        i_col = 0
        for i_evt in event_numbers:
            if not fit_selected_events and args.events >= 0 and i_evt >= args.events: break

            print ( "Event {:d}".format( i_evt ) )

            #event = df.loc[i,'Vals']
            #event = dset[i_evt].copy()
            event0 = np.array( df.loc[ i_evt, ('Waveform_ch0' if run_double_channel else 'Waveform') ] )
            print ( event0 ) 

            event1 = None
            if run_double_channel:
                event1 = np.array( df.loc[ i_evt, 'Waveform_ch1' ] )
                print ( event1 ) 
              
            #baseline = GetBaseLine(event, 200)
            val_baseline0 = baseline(event0, 200, 200)
            event0 -= val_baseline0
            if invert: event0 *= -1
            print ( event0 )
    
            if run_double_channel:
                val_baseline1 = baseline(event1, 200, 200)
                event1 -= val_baseline1
                if invert: event1 *= -1
                print ( event1 )
             
            xmin0, xmax0 = findPulseWidth(event0, 1/3)
            print ( "X(min), X(max) (Ch0) = {:d}, {:d}".format( xmin0, xmax0 ) )    

            xmin1, xmax1 = None, None
            if run_double_channel:
                xmin1, xmax1 = findPulseWidth(event1, 1/3)
                print ( "X(min), X(max) (Ch1) = {:d}, {:d}".format( xmin1, xmax1 ) )    

            #result0 = fit_event(event0, False, 400, 800, conv)
            debug_fit = False
            if fit_selected_events: debug_fit = True
            result0 = fit_event(event0, debug_fit, 400, 800, conv, 15., (0., 1950., 20., 20000., 250.))
            chi2_fit0, integral_fit_result0, x_leading_binned0, x_leading_fit0, popt_fit0, pcov_fit0, x_data_range0, event_range0 = result0
            print ( "Fit resuts Ch0" )
            print ( popt_fit0 )
            print ( pcov_fit0 ) 

            chi2_fit1, integral_fit_result1, x_leading_binned1, x_leading_fit1, popt_fit1, pcov_fit1, x_data_range1, event_range1 = None, None, None, None, None, None, None, None
            if run_double_channel:
                #result1 = fit_event(event1, False, 400, 800, conv)
                result1 = fit_event(event1, debug_fit, 400, 800, conv, 15., (0., 1950., 20., 20000., 250.))
                chi2_fit1, integral_fit_result1, x_leading_binned1, x_leading_fit1, popt_fit1, pcov_fit1, x_data_range1, event_range1 = result1
                print ( "Fit resuts Ch1" )
                print ( popt_fit1 )
                print ( pcov_fit1 ) 
            
            if fit_selected_events:
                axes[i_row, i_col].plot( ( x_data_range0 if show_channel == 0 else x_data_range1 ), 
                                         ( event_range0 if show_channel == 0 else event_range1 ), 'ko' )

            sigma_fit0, tau0_fit0, a1_fit0, tau1_fit0 = None, None, None, None
            integral_fit0 = None
            if chi2_fit0:
                sigma_fit0 = popt_fit0[2]
                tau0_fit0  = popt_fit0[4]
                a1_fit0    = popt_fit0[5]
                tau1_fit0  = popt_fit0[6]

                integral_fit0 = integral_fit_result0[0]

                if fit_selected_events and show_channel == 0:
                    model_exp_RC_res = ModelExpRCRes()

                    func_ = lambda x: model_exp_RC_res(x, *popt_fit0)
    
                    x_min_func = x_data_range0[0]
                    x_max_func = x_data_range0[-1]
                    x_data_lin = np.linspace(x_min_func,x_max_func,1000)
    
                    axes[i_row, i_col].plot( x_data_lin, func_(x_data_lin) )
            else:
                chi2_fit0, integral_fit_result0, x_leading_binned0, x_leading_fit0, popt_fit0, pcov_fit0 = None, None, None, None, None, None
                sigma_fit0, tau0_fit0, a1_fit0, tau1_fit0 = None, None, None, None
                integral_fit0 = None

            sigma_fit1, tau0_fit1, a1_fit1, tau1_fit1 = None, None, None, None
            integral_fit1 = None
            if run_double_channel and chi2_fit1:
                sigma_fit1 = popt_fit1[2]
                tau0_fit1  = popt_fit1[4]
                a1_fit1    = popt_fit1[5]
                tau1_fit1  = popt_fit1[6]

                integral_fit1 = integral_fit_result1[0]

                if fit_selected_events and show_channel != 0:
                    model_exp_RC_res = ModelExpRCRes()

                    func_ = lambda x: model_exp_RC_res(x, *popt_fit1)
    
                    x_min_func = x_data_range1[0]
                    x_max_func = x_data_range1[-1]
                    x_data_lin = np.linspace(x_min_func,x_max_func,1000)
    
                    axes[i_row, i_col].plot( x_data_lin, func_(x_data_lin) )
            else:
                chi2_fit1, integral_fit_result1, x_leading_binned1, x_leading_fit1, popt_fit1, pcov_fit1 = None, None, None, None, None, None
                sigma_fit1, tau0_fit1, a1_fit1, tau1_fit1 = None, None, None, None
                integral_fit1 = None

            #df.loc[i,'Baseline'] = baseline
            #df.loc[i,'Integral'] = integrate(event,xmin,xmax,conv)
            #df.loc[i,'Mean']  = mean(event,xmin,xmax,conv)
            #df.loc[i,'Max']  = GetMax(event)
            if run_double_channel:
                df_out = pd.concat( [ df_out, pd.DataFrame.from_records( [
                    {'Baseline_ch0' : val_baseline0, 
                     'Integral_ch0' : integrate(event0,xmin0,xmax0,conv), 
                     'Mean_ch0' : mean(event0,xmin0,xmax0,conv), 
                     'Max_ch0': GetMax(event0),
                     'FitIntegral_ch0' : integral_fit0, 
                     'BinnedLE_ch0' : x_leading_binned0, 
                     'FitLE_ch0' : x_leading_fit0, 
                     'FitChi2_ch0' : chi2_fit0, 
                     'FitSigma_ch0' : sigma_fit0, 
                     'FitTau0_ch0' : tau0_fit0, 
                     'FitTau1_ch0' : tau1_fit0, 
                     'FitA1_ch0' : a1_fit0, 
                     'Baseline_ch1' : val_baseline1, 
                     'Integral_ch1' : integrate(event1,xmin1,xmax1,conv), 
                     'Mean_ch1' : mean(event1,xmin1,xmax1,conv), 
                     'Max_ch1': GetMax(event1),
                     'FitIntegral_ch1' : integral_fit1, 
                     'BinnedLE_ch1' : x_leading_binned1, 
                     'FitLE_ch1' : x_leading_fit1, 
                     'FitChi2_ch1' : chi2_fit1, 
                     'FitSigma_ch1' : sigma_fit1, 
                     'FitTau0_ch1' : tau0_fit1, 
                     'FitTau1_ch1' : tau1_fit1, 
                     'FitA1_ch1' : a1_fit1 
                    } ] ) ], 
                    ignore_index=True 
                    ) 
            else:
                df_out = pd.concat( [ df_out, pd.DataFrame.from_records( [
                    {'Baseline' : val_baseline0, 
                     'Integral' : integrate(event0,xmin0,xmax0,conv), 
                     'Mean' : mean(event0,xmin0,xmax0,conv), 
                     'Max': GetMax(event0),
                     'FitIntegral' : integral_fit0, 
                     'BinnedLE' : x_leading_binned0, 
                     'FitLE' : x_leading_fit0, 
                     'FitChi2' : chi2_fit0, 
                     'FitSigma' : sigma_fit0, 
                     'FitTau0' : tau0_fit0, 
                     'FitTau1' : tau1_fit0, 
                     'FitA1' : a1_fit0
                    } ] ) ], 
                    ignore_index=True 
                    ) 

            i_col += 1
            if i_col >= cols:
                i_row += 1
                i_col = 0

        if run_double_channel:
            print( df_out['Integral_ch0'] )
            print( df_out['FitIntegral_ch0'] )
            print( df_out['BinnedLE_ch0'] )
            print( df_out['FitLE_ch0'] )
            print( df_out['Integral_ch1'] )
            print( df_out['FitIntegral_ch1'] )
            print( df_out['BinnedLE_ch1'] )
            print( df_out['FitLE_ch1'] )
        else:
            print( df_out['Integral'] )
            print( df_out['FitIntegral'] )
            print( df_out['BinnedLE'] )
            print( df_out['FitLE'] )

        if not fit_selected_events:
            #fig, axes = plt.subplots(2)
            #axes[0].hist( df['Mean'], bins=50 )
            #plt.hist( df_out['Integral'], bins=200 )
    
            df_out_dropna = df_out.dropna()

            if run_double_channel:
                print( df_out_dropna['Integral_ch0'] )
                print( df_out_dropna['FitIntegral_ch0'] )
                print( df_out_dropna['BinnedLE_ch0'] )
                print( df_out_dropna['FitLE_ch0'] )
                print( df_out_dropna['Integral_ch1'] )
                print( df_out_dropna['FitIntegral_ch1'] )
                print( df_out_dropna['BinnedLE_ch1'] )
                print( df_out_dropna['FitLE_ch1'] )
            else:
                print( df_out_dropna['Integral'] )
                print( df_out_dropna['FitIntegral'] )
                print( df_out_dropna['BinnedLE'] )
                print( df_out_dropna['FitLE'] )

            if run_double_channel:
                diff_BinnedLE = df_out_dropna['BinnedLE_ch1'] - df_out_dropna['BinnedLE_ch0']
                diff_FitLE = df_out_dropna['FitLE_ch1'] - df_out_dropna['FitLE_ch0']   
                print ( diff_BinnedLE )
                print ( diff_FitLE )
 
                fig, axes = plt.subplots(2)
                axes[0].hist( diff_BinnedLE, bins=50, range=(-10., 40.) )
                axes[1].hist( diff_FitLE, bins=50, range=(-10., 40.) )
            else: 
                fig, axes = plt.subplots(2, 4, figsize=(20,10))
                axes[0,0].hist( df_out_dropna['FitChi2'], bins=50, range=(0., 100.) )
                axes[0,1].hist( df_out_dropna['FitIntegral'], bins=50, range=(0., 100.e+03) )
                axes[0,2].hist( df_out_dropna['BinnedLE'], bins=50, range=(1750., 1950.) )
                axes[0,3].hist( df_out_dropna['FitLE'], bins=50, range=(1750., 1950.) )
                axes[1,0].hist( df_out_dropna['FitSigma'], bins=50, range=(0.,100.) )
                axes[1,1].hist( df_out_dropna['FitTau0'], bins=50, range=(0., 200.) )
                axes[1,2].hist( df_out_dropna['FitTau1'], bins=50, range=(0., 800.) )        
                axes[1,3].hist( df_out_dropna['FitA1'], bins=50, range=(0., 1.) )

        df_out.to_hdf('df_out.h5', key='df_out', mode='w')

        plt.show( block=True )

if __name__ == '__main__':
    main()
