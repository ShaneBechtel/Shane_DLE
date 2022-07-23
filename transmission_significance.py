import numpy as np
import matplotlib.pyplot as plt
import numpy.random as rand
from scipy.stats import norm
from IPython import embed


def trans_sig_comp(dist,data,noise_levels,gamma_vals):

    N_pix = len(data)

    weighting = norm(0,2)
    #weights = weighting.pdf(dist)
    weights = 1.0

    data_sig = (data*weights).sum()/N_pix

    sig_vals = (rand.normal(loc=0.0,scale=noise_levels,size=(10000,len(noise_levels)))*weights).sum(axis=1)/N_pix
    sig_vals.sort()
    y = np.arange(len(sig_vals))/len(sig_vals)

    data_percentile = np.where(data_sig>sig_vals)[0][-1]/len(sig_vals)

    embed()

    plt.figure(figsize=(10,7))
    plt.plot(sig_vals,y,label='Simulated Transmission Enhancement')
    plt.hlines(data_percentile,sig_vals.min(),sig_vals.max(),label='Data Significance')
    plt.legend()
    plt.show()
