#%% 
import numpy as np
import os
import scipy.io as sio
from scipy import signal
from modules.sro_estimation import DXCPPhaT
from modules.acs import ACS
from tqdm import tqdm
import matplotlib.pyplot as plt
from config import *

# %% PARAMETERS: GENERAL

db_path = PATH_DB_TRAIN  # or conf.pathDB_train
data = sio.loadmat(os.path.join(db_path, 'parameters.mat'), mat_dtype=True)
fs_Hz = data['fs_Hz'][0][0]
sigLen_s = data['sigLen_s'][0][0]
n_setups = int(data['n_setups'][0][0])
sigLen_smp = sigLen_s * fs_Hz

# %% PARAMETERS: DXCP-PhaT
parDXCPPhaT = DXCPPhaT.defaultParams
sigLen_frames = int(np.floor((sigLen_smp - parDXCPPhaT["FFTsize_dxcp"]) / parDXCPPhaT["FFTshift_dxcp"]) + 1)


# %% ESTIMATE SRO AND PRE-ACS

# Prepare Output (save results matrix-style)
SRO_Est = np.zeros((n_setups, sigLen_frames))
chi = np.zeros((n_setups, sigLen_frames))

# Loop Iteration parameters and load corresponding signals
for nn in range(n_setups):

    dxcpPhaT = DXCPPhaT()
    acs = ACS()

    # Load signals
    filename = os.path.join(db_path, f'mat_files/setup_{nn+1}.mat')
    data = sio.loadmat(filename)
    sig_z = data['sig_z']
    print(f'Processing for setup {nn+1}')

    sig_z_frames = np.lib.stride_tricks.sliding_window_view(sig_z, window_shape=(parDXCPPhaT["FFTshift_dxcp"],2))[::parDXCPPhaT["FFTshift_dxcp"]]

    for ii in tqdm(range(sig_z_frames.shape[0])):
    #for ii in tqdm(range(100)):

        res = dxcpPhaT.process_data(sig_z_frames[ii,0,:,:])
        if ii <= 2: # first 2^13 FFT analysis frame is available only after three iterations
            continue
        acs.process(res['GCSD2_avg'][:4097])

        SRO_Est[nn, ii-3] = res['SROppm_est_out']
        chi[nn, ii-3] = acs.preACS



# %% SAVE RESULTS

if True:
    results_path = os.path.join(db_path, 'results_python')
    if not os.path.isdir(results_path):
        os.makedirs(results_path)
    results_file = os.path.join(results_path, 'est_results_python.mat')

    if not os.path.isfile(results_file):
        sio.savemat(results_file, {'SRO_Est': SRO_Est, 'chi': chi})
    else:
        raise Exception('Warning: Did not save results. A file with the same name already exists!')

# %%
