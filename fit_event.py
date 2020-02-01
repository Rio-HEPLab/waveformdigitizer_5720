
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import scipy.integrate as integrate
from scipy.special import erfc

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

def fit_function(func_, p0_, bounds_, x_data, y_data):
    
    max_nfev_ = 1000 * x_data.size
    popt, pcov = curve_fit( func_, x_data, y_data, p0=p0_, bounds=bounds_, max_nfev=max_nfev_ )

    return ( popt, pcov )

def fit_event(event_processed, debug=True, i_xmin=400, i_xmax=800, conv=1.0, threshold=30., p0_def_exp_gaus_res=(0., 1950., 20., 10000., 250.)):
    #event_cpy = event.copy()

    #val_baseline = baseline(event_cpy, 200, 200)
    #if debug: 
    #    print ( "Baseline: {:.2f}".format( val_baseline ) )
    #    print ( "\n" )
    
    #invert = True
    #if invert:
    #    event_cpy -= val_baseline
    #    event_cpy *= -1
        
    #i_xmin=400
    #i_xmax=800
    
    #x_data_range = np.arange(i_xmin,i_xmax)*4.
    x_data_range = np.arange(i_xmin,i_xmax) * conv
    
    x_min = x_data_range[0]
    x_max = x_data_range[-1]
    
    if debug:
        print ( "X(min), X(max):" ) 
        print (x_min, x_max)
        print ( "\n" )
    
    #event_range = event_cpy[i_xmin:i_xmax]
    event_range = event_processed[i_xmin:i_xmax]
    
    event_range[ (event_range < 0) ] = 0
    
    if debug: 
        print ( "Waveform:" ) 
        print ( event_range )
        print ( "\n" )
        print ( np.sqrt(event_range) )
        print ( "\n" )
    
    #p0_def_exp_gaus_res = (0., 1950., 20., 10000., 250.)
    
    bounds_exp_gaus_res = ( (-np.inf,     0.,     0.,     0.,     0.),
                            ( np.inf, np.inf, np.inf, np.inf, np.inf) )
    
    if debug: 
        print ( "Exp. + Gaus. fit initial values:" ) 
        print( p0_def_exp_gaus_res, bounds_exp_gaus_res )
        print ( "\n" )
        
    try:
        if debug: print ( "Running Exp. + Gaus. fit." ) 
        popt_0, pcov_0 = fit_function(model_exp_gaus_res, p0_def_exp_gaus_res, bounds_exp_gaus_res, x_data_range, event_range)
    except (RuntimeError, ValueError) as err:
        print( err )
        print ( "\n" )
        popt_0 = np.zeros(5)
        pcov_0 = np.zeros((5,5))
        
    if debug: 
        print ( "Exp. + Gaus. fit results:" ) 
        print( popt_0, pcov_0 )
        print ( "\n" )
    
    if (popt_0 == np.zeros(5)).all() == True:
        return (None, None, None, None, None, None, x_data_range, event_range)
    
    model_exp_RC_res = ModelExpRCRes()
    
    p0_exp_RC_res = list(popt_0[0:4]) + [20.,0.50] + [ popt_0[4] ]
    
    bounds_exp_RC_res = ( (-np.inf,     0.,     0.,     0.,     0., 0.,     0.),
                          ( np.inf, np.inf, np.inf, np.inf, np.inf, 1., np.inf) )
    
    if debug: 
        print ( "Full fit initial values:" ) 
        print ( p0_exp_RC_res, bounds_exp_RC_res )
        print ( "\n" )
    
    try:
        if debug: print ( "Running full fit." ) 
        popt_1, pcov_1 = fit_function(model_exp_RC_res, p0_exp_RC_res, bounds_exp_RC_res, x_data_range, event_range)
    except (RuntimeError, ValueError) as err:
        print( err )
        print ( "\n" )
        popt_1 = np.zeros(7)
        pcov_1 = np.zeros((7,7))
        
    if debug: 
        print ( "Full fit results:" ) 
        print( popt_1, pcov_1 )
        print ( "\n" )     
    
    if (popt_1 == np.zeros(7)).all() == True: 
        return (None, None, None, None, None, None, x_data_range, event_range)
    
    x_data_lin = np.linspace(x_min,x_max,1000)
    
    pars_no_baseline_0 = [0.] + list( popt_0[1:] )
    
    func_0 = lambda x: model_exp_gaus_res(x, *pars_no_baseline_0)
    
    func_0_full = lambda x: model_exp_gaus_res(x, *popt_0)
    
    pars_no_baseline_1 = [0.] + list( popt_1[1:] )
    
    func_1 = lambda x: model_exp_RC_res(x, *pars_no_baseline_1)
    
    func_1_full = lambda x: model_exp_RC_res(x, *popt_1)
    
    if debug:
        plt.figure(figsize=(10,5))
        plt.plot( x_data_range, event_range, 'ko' )
        #plt.plot( x_data_range, model_exp_gaus_res(x_data_range, *popt) )
        plt.plot( x_data_lin, func_0_full(x_data_lin) )
        plt.plot( x_data_lin, func_1_full(x_data_lin) )
        plt.plot( x_data_lin, func_0(x_data_lin) )
        plt.plot( x_data_lin, func_1(x_data_lin) )
        
    chi2 = np.sum( ( event_range - func_1_full(x_data_range) )**2 ) / len(event_range)
    if debug: 
        print ( "Fit chi2:" ) 
        print (chi2)
        print ( "\n" )
        
    integral_result = integrate.quad( func_1, x_min, x_max )
    if debug: 
        print ( "Fit integral result:" ) 
        print (integral_result)
        print ( "\n" )
    
    #threshold = 30.
    x_threshold_binned = -1.
    x_sel_binned = x_data_range[event_range > threshold]
    if x_sel_binned.size > 0: x_threshold_binned = x_sel_binned[0]
    if debug: 
        print ( "Binned LE:" ) 
        print (x_threshold_binned)
        print ( "\n" )
    
    x_threshold_func = -1.
    x_sel_func = x_data_lin[func_1(x_data_lin) > threshold]
    if x_sel_func.size > 0: x_threshold_func = x_sel_func[0]
    if debug: 
        print ( "Fit LE:" ) 
        print (x_threshold_func)
        print ( "\n" )
    
    #return (chi2, integral_result, x_threshold_binned, x_threshold_func, popt_1, pcov_1)
    return (chi2, integral_result, x_threshold_binned, x_threshold_func, popt_1, pcov_1, x_data_range, event_range)


