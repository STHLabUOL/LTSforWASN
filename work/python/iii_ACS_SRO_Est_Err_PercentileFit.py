#%% 
import numpy as np
import os
import scipy.io as sio
from scipy.optimize import curve_fit
from modules.sro_estimation import DXCPPhaT
from modules.acs import ACS
from tqdm import tqdm
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
from mpl_toolkits.mplot3d import Axes3D

from config import *

# %% PARAMETERS: GENERAL

db_path = PATH_DB_TRAIN  # or conf.pathDB_train
data = sio.loadmat(os.path.join(db_path, 'parameters.mat'), mat_dtype=True)
fs_Hz = data['fs_Hz'][0][0]
sigLen_s = data['sigLen_s'][0][0]
n_setups = int(data['n_setups'][0][0])
sigLen_smp = sigLen_s * fs_Hz


data_results = sio.loadmat(os.path.join(db_path, 'results_python/est_results_python.mat'), mat_dtype=True)
chi = data_results['chi']
SRO_Est = data_results['SRO_Est']

# %%
chi_all = []
sro_absErr_all = []
for nn in tqdm(range(n_setups)):

    data_setup = sio.loadmat(os.path.join(db_path, f'mat_files/setup_{nn+1}.mat'), mat_dtype=True)
    SRO_true = data_setup['parAsync'][0][0][3][0][1]

    chi_ = np.squeeze(chi[nn, 100:])

    SRO_Est_ = np.squeeze(SRO_Est[nn,100:])
    SRO_absErr = np.abs(SRO_Est_-SRO_true)

    chi_all.extend(chi_)
    sro_absErr_all.extend(SRO_absErr)

chi_all = np.array(chi_all)
sro_absErr_all = np.array(sro_absErr_all)


# %% Get percentile values required for curve fitting

binwidth = 0.01
nth_prctile = 95
omit_start = 0
omit_end = 0
yNumBins = 150

bin_edges = np.arange(min(chi_all), max(chi_all), binwidth)
bin_edges = bin_edges[omit_start:]
if omit_end > 0:
    bin_edges = bin_edges[:-omit_end]

bin_centers = bin_edges[1:] - binwidth/2

sro_absErr_prctile_values = np.zeros((np.size(bin_centers)))
for ii in range(1,np.size(bin_edges)):
    select = (chi_all > bin_edges[ii-1]) & (chi_all <= bin_edges[ii])
    x = sro_absErr_all[select]
    sro_absErr_prctile_values[ii-1] = np.percentile(x, nth_prctile)


# %% Fit exp and find approx. zero-corssing chi0

def exp_func(x, a, b, c):
    return a - b*np.exp(-c*x)

def getStartPoint():
    # finds suitable starting parameters for exp. curve fitting task
    a1 = np.ones((np.size(bin_centers),1))
    a2 = np.expand_dims(-np.exp(-bin_centers), axis=1)
    A = np.concatenate((a1, a2), axis=1)
    b = np.expand_dims(np.log10(sro_absErr_prctile_values), axis=1)
    params_, _, _, _ = np.linalg.lstsq(A, b, rcond=None)

    params = []
    for p in params_:
        params.append(p[0])
    params.append(1)

    return params

params, _ = curve_fit(exp_func, bin_centers, np.log10(sro_absErr_prctile_values), p0=getStartPoint())
a, b, c = params

fit_x = np.linspace(bin_centers[0], bin_centers[-1], 10000)
fit_y = exp_func(fit_x, a, b, c)#np.power(10, exp_func(fit_x, a, b, c))
fit_z = np.zeros_like(fit_y)

chi0 = fit_x[np.argmin(np.abs(np.log10(fit_y)))]
print(f'Curve zero-crossing at approx. chi0={chi0}')


# %% Plot bivariate histogram

def movmean(data, window_size):
    moving_avg = np.zeros(len(data))
    # Calculate moving average, adjusting window size at the edges
    for i in range(len(data)):
        start = max(0, i - window_size // 2)
        end = min(len(data), i + window_size // 2 + 1)
        moving_avg[i] = np.mean(data[start:end])
    return moving_avg


#%%

yEdges = np.linspace(np.log10(min(sro_absErr_all)), np.log10(max(sro_absErr_all)), yNumBins)
H, xEdges, yEdges = np.histogram2d(chi_all, np.log10(sro_absErr_all), bins=(bin_edges, yEdges))                     
H = H/np.size(H)
yCenters = yEdges[1:] - (yEdges[1]-yEdges[0])/2
xCenters = xEdges[1:] - (xEdges[1]-xEdges[0])/2

fig = plt.figure()
ax1 = fig.add_subplot(111, projection='3d')
#ax1.set_yscale('log')
ylim = (-3, 2)
xlim = (0.04, 0.15)
zlim = (0, 0.2)
for ii in range(1,H.shape[0]):

    x = xCenters[ii] * np.ones(np.size(yCenters))
    y = yCenters#np.power(10, yCenters)
    z = movmean(H[ii,:], 5)#np.convolve(H[ii,:], np.ones(5)/5, mode='valid') #moving average for smoothing
    z[0] = 0
    z[-1] = 0

    # crop polygon to avoid outline clipping
    select = (y > ylim[0]) & (y < ylim[1])
    x = x[select]
    y = y[select]
    z = z[select]

    verts = [list(zip(x, y, z))]
    ax1.add_collection3d(Poly3DCollection(verts, alpha=0.4, edgecolor='k'))

    ax1.set_xlim(xlim)
    ax1.set_ylim(ylim)
    ax1.set_zlim(zlim)
    ax1.view_init(elev=30, azim=20)

plt.plot(fit_x, fit_y, color='tab:orange', linewidth=2, label='UAE$_{\epsilon,0}$')
plt.plot(bin_centers, np.log10(sro_absErr_prctile_values), marker='x', linestyle='None', color='tab:orange', label='{$\chi$;AE$_{\epsilon}^{95}$}')
plt.plot([xlim[0], xlim[1]], [0, 0], linestyle='--', color='k', label='AE$_{\epsilon,0}$=1ppm')
plt.plot([chi0, chi0], [ylim[0], ylim[1]], linestyle='-', color='k', label='$\chi_0$')


plt.xlabel('$\chi$')
plt.ylabel('log|$\epsilon-\hat{\epsilon}$| [ppm]')
plt.legend()
plt.show()



# %%
