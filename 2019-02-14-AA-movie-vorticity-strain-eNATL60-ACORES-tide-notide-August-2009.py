


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




def plot_comp_vort(k):
    fig=plt.figure(figsize=(20,15))

    fcurltide = xr.open_mfdataset(tfilescurltide)
    curltide=fcurltide['socurloverf'][k]
    navlat= fcurltide['nav_lat'][k]
    navlon= fcurltide['nav_lon'][k]
    plot_surf(121,curltide,navlon,navlat,-1,1,div_cmap,'Surface vorticity over f, with tide',str(curltide.time_counter.values))

    fcurlnotide = xr.open_mfdataset(tfilescurlnotide)
    curlnotide=fcurlnotide['socurloverf'][0]
    plot_surf(122,curlnotide,navlon,navlat,-1,1,div_cmap,'Surface vorticity over f, no tide',str(curlnotide.time_counter.values))
    plt.savefig('plots/surfcurloverf_eNATL60_tide-notide_ACO_'+str(k)+'.png')


def plot_comp_strain(k):
    fig=plt.figure(figsize=(20,15))

    fstraintide = xr.open_mfdataset(tfilesstraintide)
    straintide=fstraintide['sostrainoverf'][k]
    navlat= fstraintide['nav_lat'][k]
    navlon= fstraintide['nav_lon'][k]
    plot_surf(121,straintide,navlon,navlat,-1,1,div_cmap,'Surface strain over f, with tide',str(straintide.time_counter.values))

    fstrainnotide = xr.open_mfdataset(tfilesstrainnotide)
    strainnotide=fstrainnotide['sostrainoverf'][0]
    plot_surf(122,strainnotide,navlon,navlat,-1,1,div_cmap,'Surface strain over f, no tide',str(strainnotide.time_counter.values))
    plt.savefig('plots/surfstrainoverf_eNATL60_tide-notide_ACO_'+str(k)+'.png')


for k in np.arange(100):
	plot_comp_vort(k)
	plot_comp_strain(k)

!convert -delay 120 -loop 0  surfcurloverf_eNATL60_tide-notide_ACO_?.png surfcurloverf_eNATL60_tide-notide_ACO_??.png surfcurloverf_eNATL60_tide-notide_ACO.gif
!convert -delay 120 -loop 0  surfstrainoverf_eNATL60_tide-notide_ACO_?.png surfstrainoverf_eNATL60_tide-notide_ACO_??.png surfstrainoverf_eNATL60_tide-notide_ACO.gif
