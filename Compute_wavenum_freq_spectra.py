import sys, glob
import operator
import numpy as np
import xarray as xr
import numpy as np
import xscale.spectral.fft as xfft

sys.path.insert(0, "/home/ajayi/EnergyCasade/Scripts/Wavnum_Freq_Spec/")
import Wavenum_freq_spec_func as wfs

sys.path.insert(0, "/home/ajayi/EnergyCasade/Scripts/")
from SmallBox import smallbox
for rbox in smallbox:
    box = rbox
    
# - define file path
Daily_data_dir = '/mnt/meom/workdir/ajayi/SSH/' 
Daily_12Months_data_file = Daily_data_dir + 'NATL60-CJM165_y*.1d.SSH.nc'
Daily_JFM_data_file = Daily_data_dir + 'NATL60-CJM165_y2013m0[1-3]*.1d.SSH.nc'
Daily_JAS_data_file = Daily_data_dir + 'NATL60-CJM165_y2013m0[7-9]*.1d.SSH.nc'

Hourly_data_dir = '/mnt/meom/MODEL_SET/NATL60/NATL60-CJM165-S/1h/SSH/'
Hourly_JFM_data_file = Hourly_data_dir + 'NATL60-CJM165_y2013m0[1-3]*.1h_SSH.nc'
Hourly_JAS_data_file = Hourly_data_dir + 'NATL60-CJM165_y2013m0[7-9]*.1h_SSH.nc'

# - Save dataset to this folder
OutputFolder = '/home/ajayi/EnergyCasade/Data/Wavenum_freq_spectrum/'

# - Define filr name and Outputfile
Datafilepath = Hourly_JAS_data_file
OutputFile = 'SSH_JAS_wavenum_freq_spec_from_hourly_outputs.nc'


# - open dataset
data = xr.open_mfdataset(Datafilepath,chunks={'time_counter': 100, 'x':200})
ssh = data['sossheig']
# - extract the needed box
sshbox = ssh[:,box.jmin:box.jmax,box.imin:box.imax]
# - get dx and dy
dx,dy = wfs.get_dx_dy(sshbox[0])
#... Remove NaN ...
ssh_No_NaN = sshbox.interpolate_na(dim='y')


#... Detrend data in all dimension ...
ssh_dtr = wfs.detrendn(ssh_No_NaN,axes=[0,1,2])
#... Apply hanning windowing ...') 
ssh_wdw = wfs.apply_window(ssh_dtr, ssh_dtr.dims, window_type='hanning')
#... Get fourier transform ... 
sshhat = xfft.fft(ssh_wdw, dim=('time_counter', 'x', 'y'), dx={'x': dx, 'y': dx}, sym=True)
#... Get power spectra density ... 
ssh_psd = xfft.psd(sshhat)
#... Get frequency and wavenumber ... 
frequency,kx,ky = wfs.get_f_kx_ky(sshhat)
#... Get istropic wavenumber ... 
wavenumber,kradial = wfs.get_wavnum_kradial(kx,ky)
#... Get numpy array ... 
ssh_psd_np = ssh_psd.values
#... Get 2D frequency-wavenumber field ... 
wavenum_freq_spectrum = wfs.get_f_k_in_2D(kradial,wavenumber,ssh_psd_np)

# Save to Netscdf file
# - build dataarray
spectrum_da = xr.DataArray(wavenum_freq_spectrum,dims=['frequency','wavenumber'],name="wavenum_freq_spec",coords=[frequency,wavenumber])
spectrum_da.attrs['Name'] = OutputFile
# - To dataset then to Netcdf
spectrum_da.to_dataset().to_netcdf(path=OutputFolder+OutputFile,mode='w',engine='scipy')