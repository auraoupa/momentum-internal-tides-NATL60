


## path for mdules

import sys
sys.path.insert(0,"/home/albert/lib/python")

import numpy as np
import xarray as xr

import cartopy.crs as ccrs
import cartopy.feature as cfeature
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER

import numpy.ma as ma

import matplotlib.cm as mplcm

seq_cmap = mplcm.Blues
div_cmap = mplcm.seismic

import matplotlib.pyplot as plt

import seaborn as sns
sns.set(color_codes=True)


## Dataset

dirtide="/mnt/meom/MODEL_SET/eNATL60/eNATL60-BLBT02-S/1h/ACO/"
dirnotide="/mnt/meom/MODEL_SET/eNATL60/eNATL60-BLB002-S/1h/ACO/"



filescurltide="eNATL60ACO-BLBT02_1h_*_socurloverf_*.nc"
filesstraintide="eNATL60ACO-BLBT02_1h_*_sostrainoverf_*.nc"

tfilescurltide=dirtide+filescurltide
tfilesstraintide=dirtide+filesstraintide

filescurlnotide="eNATL60ACO-BLB002_1h_*_socurloverf_*.nc"
filesstrainnotide="eNATL60ACO-BLB002_1h_*_sostrainoverf_*.nc"

tfilescurlnotide=dirnotide+filescurlnotide
tfilesstrainnotide=dirnotide+filesstrainnotide


#!ls $tfilescurlnotide
#!ls $tfilesstraintide




fcurltide = xr.open_mfdataset(tfilescurltide)
curltide=fcurltide['socurloverf'][0:10]
navlat= fcurltide['nav_lat'][0]
navlon= fcurltide['nav_lon'][0]

fcurlnotide = xr.open_mfdataset(tfilescurlnotide)
curlnotide=fcurlnotide['socurloverf'][0]

fstraintide = xr.open_mfdataset(tfilesstraintide)
straintide=fstraintide['sostrainoverf'][0:10]

fstrainnotide = xr.open_mfdataset(tfilesstrainnotide)
strainnotide=fstrainnotide['sostrainoverf'][0:10]




box=(-31,-28,33,36)
domain=(box[0]<navlon)*(navlon<box[1])*(box[2]<navlat)*(navlat<box[3])
where=np.where(domain)

lats=navlat[where]
lons=navlon[where]

ind = np.unravel_index(np.argmin(lats, axis=None), lats.shape)
jmin = where[0][ind[0]]
ind = np.unravel_index(np.argmax(lats, axis=None), lats.shape)
jmax = where[0][ind[0]]
ind = np.unravel_index(np.argmin(lons, axis=None), lons.shape)
imin = where[1][ind[1]]
ind = np.unravel_index(np.argmax(lons, axis=None), lons.shape)
imax = where[1][ind[1]]

curltidebox=curltide[:,jmin:jmax+1,imin:imax+1].stack(z=('x', 'y','time_counter'))
curlnotidebox=curlnotide[:,jmin:jmax+1,imin:imax+1].stack(z=('x', 'y','time_counter'))
straintidebox=straintide[:,jmin:jmax+1,imin:imax+1].stack(z=('x', 'y','time_counter'))
strainnotidebox=strainnotide[:,jmin:jmax+1,imin:imax+1].stack(z=('x', 'y','time_counter'))

weights_curltide = np.ones_like(curltidebox)/float(len(curltidebox))
weights_curlnotide = np.ones_like(curlnotidebox)/float(len(curlnotidebox))
weights_straintide = np.ones_like(straintidebox)/float(len(straintidebox))
weights_strainnotide = np.ones_like(strainnotidebox)/float(len(strainnotidebox))



# In[102]:


fig = plt.figure(figsize=(18.0, 12.0))
axes1 = fig.add_subplot(1, 1, 1)

axes1.hist(curltidebox,100, alpha = 0.5,range=(-1,1),color='r', weights=weights_curltide,label='tide')
axes1.hist(curlnotidebox,100, alpha = 0.5,range=(-1,1),color='b', weights=weights_curlnotide, label='no tide')
plt.ticklabel_format(style='sci', axis='x', scilimits=(0,0))
plt.xlabel('surface curl/f')
plt.legend()


# In[103]:


fig = plt.figure(figsize=(18.0, 12.0))
axes1 = fig.add_subplot(1, 1, 1)
sns.distplot(curltidebox, color="r",label='tide',kde=True,hist=True,norm_hist=True)
sns.distplot(curlnotidebox, color="b",label='no tide',kde=True,hist=True,norm_hist=True)
plt.legend()


# In[108]:


fig = plt.figure(figsize=(18.0, 12.0))
axes1 = fig.add_subplot(1, 1, 1)

axes1.hist(straintidebox,100, alpha = 0.5,range=(0,1),color='r', weights=weights_straintide,label='tide')
axes1.hist(strainnotidebox,100, alpha = 0.5,range=(0,1),color='b', weights=weights_strainnotide, label='no tide')
plt.ticklabel_format(style='sci', axis='x', scilimits=(0,0))
plt.xlabel('surface strain/f')
plt.legend()


# In[109]:


fig = plt.figure(figsize=(18.0, 12.0))
axes1 = fig.add_subplot(1, 1, 1)
sns.distplot(straintidebox, color="r",label='tide',kde=True,hist=True,norm_hist=True)
sns.distplot(strainnotidebox, color="b",label='no tide',kde=True,hist=True,norm_hist=True)
plt.legend()

