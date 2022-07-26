import numpy as np
import matplotlib.pyplot as plt
import numpy.random as rand
from astropy.table import Table
from scipy.stats import norm
from IPython import embed


def trans_sig_comp(dist, data, noise_levels, boost_func, skewer_path):
    N_pix_data = len(data)

    # weighting = norm(0,2)
    # weights = weighting.pdf(dist)

    data_sig = data.sum() / N_pix_data

    # Skewer Realizations

    params = Table.read(skewer_path, hdu=1)
    skewers = Table.read(skewer_path, hdu=2)

    Lbox_cMpc = params['Lbox'] / params['lit_h']
    drpix = Lbox_cMpc / params['Ng']
    rvec_cMpc = np.arange(params['Ng']) * drpix

    qso_idx = int(len(rvec_cMpc) / 2)
    rvec_cMpc -= rvec_cMpc[qso_idx]

    sim_boost = boost_func(rvec_cMpc)

    #tau_skewers = skewers['TAU'] / sim_boost
    tau_skewers = skewers['TAU']
    trans_skewers = np.exp(-tau_skewers)

    summary_range = (rvec_cMpc >= dist[0]) & (rvec_cMpc <= dist[-1])
    trans_stat = trans_skewers[:, summary_range].sum(axis=1) / summary_range.sum()

    sig_vals = rand.normal(loc=0.0, scale=noise_levels, size=(1000, len(noise_levels))).sum(axis=1) / len(noise_levels)
    sig_vals += trans_stat
    sig_vals.sort()

    y = np.arange(len(sig_vals)) / len(sig_vals)

    try:
        data_percentile = np.where(data_sig > sig_vals)[0][-1] / len(sig_vals)
    except:
        data_percentile = 0

    embed()

    (mu, sigma) = norm.fit(sig_vals)

    plt.figure(figsize=(10, 7))
    plt.plot(sig_vals, y, label='Simulated Transmission Enhancement')
    plt.vlines(data_sig, 0, 1, 'r', label='Data Significance')
    plt.vlines(sig_vals[int(len(sig_vals)*0.84)], 0, 1, color='g', label='1 Sigma')
    plt.vlines(sig_vals[int(len(sig_vals)*0.16)], 0, 1, color='g')
    plt.legend()
    plt.show()
    plt.close()

    vals, bins, patches = plt.hist(sig_vals, bins=30, label='Simulated Transmission Enhancement')
    plt.vlines(data_sig, 0, vals.max(), color='r', label='Data Significance')
    plt.vlines(mu + sigma, 0, vals.max(), color='g', label='1 Sigma Statistic Value')
    plt.vlines(mu - sigma, 0, vals.max(), color='g')
    plt.legend()
    plt.show()
    plt.close()

    return rvec_cMpc, trans_skewers
