from scipy.optimize import curve_fit
import numpy as np
import matplotlib.pyplot as plt

import sklearn.metrics

def data_fit(TIC,normalizer):
    normalizedLogParams, normalizedLogParamCov = curve_fit(lognormal, TIC[0], TIC[1], p0=(1.0, 0.0,1.0),bounds=([0.,0., 0.], [np.inf, np.inf, np.inf]),method='trf')#p0=(1.0,3.0,0.5,0.1) ,**kwargs
    # popt = np.around(normalizedLogParams, decimals=1);
    popt = normalizedLogParams

    auc = popt[0]
    mu = popt[1]
    sigma = popt[2]
    mtt = np.exp(mu+(sigma**2/2))
    wholeCurve = lognormal(TIC[0], auc, mu, sigma)
    print("brute force auc:", sklearn.metrics.auc(TIC[0], wholeCurve))
    print("mu:", mu)
    print("sigma:",sigma)
    tp = np.exp(mu - (sigma**2))
    pe = np.max(wholeCurve)
    params = np.array(np.around(np.array([pe, auc, tp, mtt, 0]), decimals=1))
    params[0] = np.around(pe, decimals=3)

    wholeCurve *= normalizer;
    return params, popt, wholeCurve;

def lognormal(x, auc, mu, sigma):      
    curve_fit=(auc/(2.5066*sigma*x))*np.exp((-1/2)*(((np.log(x)-mu)/sigma)**2)) 
    return np.nan_to_num(curve_fit)