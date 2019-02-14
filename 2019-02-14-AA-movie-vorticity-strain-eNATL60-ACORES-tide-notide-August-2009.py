
# coding: utf-8

# # Movies of surface vorticity and strain in eNATL60 simulation with tides and no tides in ACORES region 
# 
# __author__ : Aurélie Albert, Julien Le Sommer (MEOM), Andy Hogg (ANU)
# 
# __context__ : paper "On the momentum flux of internal tides" by Shakespeare & Hogg
# 
# __date__ : February 2019
# 
# __purpose__ : distribution of surface vorticity and strain values in eNATL60 simulation with tides and no tides in the ACORES region 
# 
# __detailed description__ : 
# Vorticity is defined as :
# 
# $$curl = \partial{_x}{v} - \partial{_y}{u} $$
# 
# Strain is defined as :
# 
# $$strain = \sqrt{ (\partial{_x}{v}+\partial{_y}{u})^2 + (\partial{_x}{u}-\partial{_y}{v})^2 } $$
# 
# Both quantities are scale by f.
# 
# __practical steps__ :
# 
#   * Input data are hourly surface currents from eNATL60-BLB002 simulation (no-tide) and eNTAL60-BLBT02 simulation (tide) 
#   
#   * Surface vorticity and strain over f are computed with cdfcurl and cdfstrain cdftool : https://github.com/meom-group/CDFTOOLS) for the first two weeks of August 2009
#   
#   
# __external libraries needed to run this script__ : 
# 
#  
# __licence__ : This work is licensed under a <a rel="license" href="http://creativecommons.org/licenses/by/4.0/">Creative Commons Attribution 4.0 International License</a>.

# In[1]:


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

get_ipython().magic(u'matplotlib inline')


# In[3]:


## Dataset

dirtide="/mnt/albert/equipes/IGE/meom/MODEL_SET/eNATL60/eNATL60-BLBT02-S/1h/ACO/"
dirnotide="/mnt/albert/equipes/IGE/meom/MODEL_SET/eNATL60/eNATL60-BLB002-S/1h/ACO/"



# In[4]:


filescurltide="eNATL60ACO-BLBT02_1h_*_socurloverf_*.nc"
filesstraintide="eNATL60ACO-BLBT02_1h_*_sostrainoverf_*.nc"

tfilescurltide=dirtide+filescurltide
tfilesstraintide=dirtide+filesstraintide

filescurlnotide="eNATL60ACO-BLB002_1h_*_socurloverf_*.nc"
filesstrainnotide="eNATL60ACO-BLB002_1h_*_sostrainoverf_*.nc"

tfilescurlnotide=dirnotide+filescurlnotide
tfilesstrainnotide=dirnotide+filesstrainnotide


# In[5]:


def plot_surf(sub,data,lon,lat,vmin,vmax,cmap,title,date):
    
    ax = plt.subplot(sub,projection=ccrs.PlateCarree(central_longitude=0))
    ax.set_extent((-36, -26, 25, 40))
    land = cfeature.GSHHSFeature(scale='intermediate',
                                 levels=[1],
                                 facecolor=cfeature.COLORS['land'])
    ax.add_feature(land)
    gl = ax.gridlines(draw_labels=True, linestyle=':', color='black',
                      alpha=0.5)
    gl.xlabels_top = False
    gl.ylabels_right = False
    gl.xformatter = LONGITUDE_FORMATTER
    gl.yformatter = LATITUDE_FORMATTER
    ax.tick_params('both',labelsize=22)

    plt.pcolormesh(lon,lat,data,cmap=cmap,vmin=vmin,vmax=vmax )
    plt.colorbar(orientation='horizontal',pad=0.1)
    plt.title(title)
    ax.text(0.57, -0.08, date, transform=ax.transAxes)

    ax.plot([-31, -28], [33, 33],color='black',linewidth=4)
    ax.plot([-31, -28], [36, 36],color='black',linewidth=4)
    ax.plot([-31, -31], [33, 36],color='black',linewidth=4)
    ax.plot([-28, -28], [33, 36],color='black',linewidth=4)


# In[7]:


def plot_comp_vort(k):
    fig=plt.figure(figsize=(20,15))

    fcurltide = xr.open_mfdataset(tfilescurltide)
    curltide=fcurltide['socurloverf'][k]
    navlat= fcurltide['nav_lat'][k]
    navlon= fcurltide['nav_lon'][k]
    plot_surf(121,curltide,navlon,navlat,-1,1,div_cmap,'Surface vorticity over f, with tide',str(curl15.time_counter.values))

    fcurlnotide = xr.open_mfdataset(tfilescurlnotide)
    curlnotide=fcurlnotide['socurloverf'][0]
    plot_surf(122,curlnotide,navlon,navlat,-1,1,div_cmap,'Surface vorticity over f, no tide',str(curl15.time_counter.values))
    plt.title('surfcurloverf_eNATL60_tide-notide_ACO_'+str(k)+'.png')


# In[10]:


def plot_comp_strain(k):
    fig=plt.figure(figsize=(20,15))

    fstraintide = xr.open_mfdataset(tfilesstraintide)
    straintide=fstraintide['sostrainoverf'][k]
    navlat= fstraintide['nav_lat'][k]
    navlon= fstraintide['nav_lon'][k]
    plot_surf(121,straintide,navlon,navlat,-1,1,div_cmap,'Surface strain over f, with tide',str(curl15.time_counter.values))

    fstrainnotide = xr.open_mfdataset(tfilesstrainnotide)
    strainnotide=fstrainnotide['sostrainoverf'][0]
    plot_surf(122,strainnotide,navlon,navlat,-1,1,div_cmap,'Surface strain over f, no tide',str(curl15.time_counter.values))
    plt.title('surfstrainoverf_eNATL60_tide-notide_ACO_'+str(k)+'.png')


# In[8]:


plot_comp_vort(0)

