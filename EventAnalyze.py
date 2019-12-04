import argparse
import logging
logger = logging.getLogger("EventAnalyze")
logger.setLevel(logging.DEBUG)

import numpy as np
import h5py
import matplotlib.pyplot as plt

from scipy.optimize import curve_fit
import scipy.integrate as integrate
from scipy.special import erfc

class Event():
	def __init__(self, _integral, _threshold, _func, _x_data_range, _event_range):
		self.integral = _integral
		self.threshold = _threshold
		self.get = _func
		self.x_data_range = _x_data_range
		self.event_range = _event_range

def baseline(vals, n_low, n_high):
    vals_low = vals[:n_low]
    vals_high = vals[(vals.size - n_high):-1]
    vals_mean = ( np.sum(vals_low) + np.sum(vals_high) ) / ( vals_low.size + vals_high.size )
    return vals_mean

def model_exp_gaus_res(x, *pars):
    C = pars[0]
    t0 = pars[1]
    sigma = pars[2]
    A = pars[3]
    tau = pars[4]

    val = A * ( 1./(2*tau) ) * np.exp( (sigma**2)/(2*(tau**2)) - (x - t0)/tau ) * erfc( sigma/(np.sqrt(2.)*tau) - (x - t0)/(np.sqrt(2.)*sigma) )
    val += C

    return val

def model_exp_RC_res(x, *pars):

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

def fit_function(func_, p0_, bounds_, x_data, y_data):

    popt, pcov = curve_fit( func_, x_data, y_data, p0=p0_, bounds=bounds_ )

    return ( popt, pcov )

def fit_event(event, debug=False, i_xmin=400, i_xmax=800):
    event_cpy = event.copy()

    val_baseline = baseline(event_cpy, 200, 200)
    if debug:
        logger.debug ( "Baseline: {:.2f}".format( val_baseline ) )

    invert = True
    if invert:
        event_cpy -= val_baseline
        event_cpy *= -1

    x_data_range = np.arange(i_xmin,i_xmax)*4.

    x_min = x_data_range[0]
    x_max = x_data_range[-1]

    if debug:
        logger.debug ("x_min e x_max: {:.2f} {:.2f}".format(x_min, x_max))

    event_range = event_cpy[i_xmin:i_xmax]

    event_range[ (event_range < 0) ] = 0

    if debug:
        logger.debug ( "event_range")
        logger.debug (event_range)
        logger.debug ( "square of event_range:")
        logger.debug ( np.sqrt(event_range))

    p0_def_exp_gaus_res = (0., 1950., 20., 10000., 250.)

    bounds_exp_gaus_res = ( (-np.inf,     0.,     0.,     0.,     0.),
                            ( np.inf, np.inf, np.inf, np.inf, np.inf) )

    if debug:
        logger.debug("p0 : {0}\n bounds: {1} ".format( p0_def_exp_gaus_res, bounds_exp_gaus_res) )

    try:
        popt_0, pcov_0 = fit_function(model_exp_gaus_res, p0_def_exp_gaus_res, bounds_exp_gaus_res, x_data_range, event_range)
    except (RuntimeError, ValueError) as err:
        logger.error( err )
        popt_0 = np.zeros(5)
        pcov_0 = np.zeros((5,5))

    if debug:
        logger.debug("pop_t0: {0}\n pcov_0: {1}".format(popt_0, pcov_0) )

    if (popt_0 == np.zeros(5)).all() == True:
        #return (None, None, None, None, None, None, x_data_range, event_range)
        return None

    p0_exp_RC_res = list(popt_0[0:4]) + [20.,0.50] + [ popt_0[4] ]

    bounds_exp_RC_res = ( (-np.inf,     0.,     0.,     0.,     0., 0.,     0.),
                          ( np.inf, np.inf, np.inf, np.inf, np.inf, 1., np.inf) )

    if debug:
        logger.debug ("p0_exp_RC_res: {0}\n bounds_exp_RC_res: {1}".format(p0_exp_RC_res, bounds_exp_RC_res) )

    try:
        popt_1, pcov_1 = fit_function(model_exp_RC_res, p0_exp_RC_res, bounds_exp_RC_res, x_data_range, event_range)
    except (RuntimeError, ValueError) as err:
        logger.error( err )
        popt_1 = np.zeros(7)
        pcov_1 = np.zeros((7,7))

    if debug:
        logger.debug("popt_1: {0}\n pcov_1: {1}".format( popt_1, pcov_1 ))

    if (popt_1 == np.zeros(7)).all() == True:
        #return (None, None, None, None, None, None, x_data_range, event_range)
        return None

    x_data_lin = np.linspace(x_min,x_max,1000)

    pars_no_baseline_0 = [0.] + list( popt_0[1:] )

    func_0 = lambda x: model_exp_gaus_res(x, *pars_no_baseline_0)

    func_0_full = lambda x: model_exp_gaus_res(x, *popt_0)

    pars_no_baseline_1 = [0.] + list( popt_1[1:] )

    func_1 = lambda x: model_exp_RC_res(x, *pars_no_baseline_1)

    func_1_full = lambda x: model_exp_RC_res(x, *popt_1)

    chi2 = np.sum( ( event_range - func_1_full(x_data_range) )**2 ) / len(event_range)
    if debug:
        logger.debug("chi2: {:.2f}".format(chi2))

    integral_result = integrate.quad( func_1, x_min, x_max )
    if debug:
        logger.debug ("integral: {0}".format(integral_result))

    threshold = 30.
    x_threshold_binned = -1.
    x_sel_binned = x_data_range[event_range > threshold]
    if x_sel_binned.size > 0: x_threshold_binned = x_sel_binned[0]
    if debug:
        logger.debug ("treshold {0}".format(x_threshold_binned))

    x_threshold_func = -1.
    x_sel_func = x_data_lin[func_1(x_data_lin) > threshold]
    if x_sel_func.size > 0: x_threshold_func = x_sel_func[0]
    if debug:
        logger.debug ("treshold_func: {0}".format(x_threshold_func))

    return Event(integral_result, x_threshold_func, func_1_full, x_data_range, event_range)

def main():

    parser = argparse.ArgumentParser(description = 'Programa que recebe waveforms e extrai suas informações')
    parser.add_argument('-df', action = 'store', dest = 'waveforms', required = True, help = 'Arquivo waveform do pandas' )
    parser.add_argument('-b', action = 'store_true', dest = 'debug', required = False, help = 'Flag de debug' )
    arguments = parser.parse_args()

    logging.basicConfig(filename='EventAnalyze.log', filemode='w')
    with h5py.File(arguments.waveforms,'r') as f:

        events = f['Vals']
        event = events[0]

        fig, axes = plt.subplots(2)
 
        result = fit_event(event[0], arguments.debug)
 
        axes[0].plot( result.x_data_range, result.event_range, 'ko' )
        axes[0].plot( result.x_data_range, result.get(result.x_data_range) )
 
        result = fit_event(event[1], arguments.debug)
 
        axes[1].plot( result.x_data_range, result.event_range, 'ko' )
        axes[1].plot( result.x_data_range, result.get(result.x_data_range) )
        plt.show()
 
        return

if __name__ == '__main__':
    main()
