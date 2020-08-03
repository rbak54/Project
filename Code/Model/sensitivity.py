"""carrying out a sensitivity analysis for R0"""

import sympy
from sympy import *
import matplotlib.pylab as plot
sigma, mu, c_r, epsilon, d, h, gamma, p = var('sigma, mu, c_r, epsilon, d, h, gamma, p ',real=True, positive = true)
q = var('q',real=True, positive = False)
R0 = (sigma / (sigma + mu)) * c_r * epsilon *  (( (-1 / q) * (1 - exp( q * d)))/((- 1 / q) * (1 - exp( q * d)) + h)) / ((mu + gamma)*( 1 / ( 1 - p )))



def diff_func(equation,variable,parms):
    """finds partial derivative over equation for particular variable"""
    dR0dv=diff(equation,variable)
    dR0dv_proportion=dR0dv/R0
    for i in dR0dv_proportion.free_symbols: 
        dR0dv_proportion=dR0dv_proportion.subs(i,parms[str(i)])
    return(dR0dv_proportion)


def loop_diff(parms,equation):
    """loop through parameters in equation(R0), call diff_func and store importance of each parameter"""
    importance_dictionary={}
    for k in equation.free_symbols:
        if abs_req==True:
            importance_dictionary[str(k)]=abs(diff_func(R0,k,parms))
        else:
            importance_dictionary[str(k)]=diff_func(R0,k,parms)
    return(importance_dictionary)


def plot_parm_imp(parms,equation):
    """plot for one set of variables"""
    rank=loop_diff(parms,equation)
    #sorted(rank, key=rank.get,reverse=True)
    f1=plot.figure()
    variable= [i for i in rank]
    value=[rank[str(i)] for i in rank]
    plot.bar(variable,height=value,width=0.8)
    plot.show()

#rank=loop_diff(parm_values,R0)
#sorted(ranabs_req=False
abs_req=False

parm_values ={"sigma": 0.68, "mu": 2.06e-5, "c_r": 18.54 , "epsilon": 0.05, "d": 4/24, "h": 0.25/24, "gamma": 0.25, "p": 0.005, "q":-39.87 }

plot_parm_imp(parms ={"sigma": 0.68, "mu": 2.06e-5, "c_r": 18.54 , "epsilon": 0.05, "d": 4/24, "h": 0.25/24, "gamma": 0.25, "p": 0.005, "q":-39.87 },equation=R0)
plot_parm_imp(parms ={"sigma": 0.68, "mu": 2.06e-5, "c_r": 18.54 , "epsilon": 0.1, "d": 4/24, "h": 0.25/24, "gamma": 0.25, "p": 0.005, "q":-39.87 }
,equation=R0)
plot_parm_imp(parms ={"sigma": 0.68, "mu": 2.06e-5, "c_r": 18.54 , "epsilon": 0.05, "d": 0.25/24, "h": 0.25/24, "gamma": 0.25, "p": 0.005, "q":-39.87 }
,equation=R0)