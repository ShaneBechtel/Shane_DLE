import numpy as np
import matplotlib.pyplot as plt
import numpy.random as rand
from scipy.stats import norm
from IPython import embed


def trans_sig_comp(dist,data,noise_levels,gamma_vals):

    N_pix = len(data)

    # weighting = norm(0,2)
    # weights = weighting.pdf(dist)
    weights = 1.0

    data_sig = (data*weights).sum()/N_pix

    qso_trans = np.exp(-1/gamma_vals)

    sig_vals = (qso_trans + rand.normal(loc=0.0,scale=noise_levels,size=(10000,len(noise_levels)))) * weights
    sig_vals = sig_vals.sum(axis=1)/N_pix
    sig_vals.sort()
    y = np.arange(len(sig_vals))/len(sig_vals)

    try:
        data_percentile = np.where(data_sig>sig_vals)[0][-1]/len(sig_vals)
    except:
        data_percentile = 0

    embed()

    (mu, sigma) = norm.fit(sig_vals)

    plt.figure(figsize=(10,7))
    plt.plot(y,sig_vals,label='Simulated Transmission Enhancement')
    plt.vlines(data_percentile,sig_vals.min(),sig_vals.max(),'r',label='Data Significance')
    plt.vlines(0.84,sig_vals.min(),sig_vals.max(),color='g',label='1 Sigma')
    plt.vlines(0.16, sig_vals.min(), sig_vals.max(), color='g')
    plt.legend()
    plt.show()
    plt.close()


    vals,bins,patches = plt.hist(sig_vals, bins=30, label='Simulated Transmission Enhancement')
    plt.vlines(data_sig, 0, vals.max(), color='r',label='Data Significance')
    plt.vlines(mu+sigma,0,vals.max(),color='g',label='1 Sigma Statistic Value')
    plt.vlines(mu-sigma,0,vals.max(),color='g',label='1 Sigma Statistic Value')
    plt.legend()
    plt.show()
    plt.close()
