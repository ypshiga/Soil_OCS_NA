#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon May  4 21:02:29 2020

@author: yshiga
"""
from mpl_toolkits.basemap import Basemap, cm
import numpy as np
import matplotlib.pyplot as plt
from netCDF4 import Dataset


data_dir = "/Users/yshiga/Documents/Data/NARR/"
cmap_vals = plt.cm.get_cmap('RdYlGn_r')
# Plotting code to create animation showing soil ocs fluxes at NARR 
# resolution (32km and 3 hourly)
i=2006
for j in range(7,8):

    # for i in range(2006,2007):  # loop through years 2006 - 2017
    #     for j in range(1,13):  # loop through months 1-12
    filename_ocs = data_dir + "soil_ocs/soil_ocs_combo.{0}{1:02d}.nc".format(i, j)
    filename_t = data_dir + "tsoil/tsoil.{0}{1:02d}.nc".format(i, j)
    
    nc_t = Dataset(filename_t, 'r')  # open soil temp for given month
    nc_ocs = Dataset(filename_ocs, 'r')  # open soil temp for given month
    COS_native_combo = nc_ocs.variables['COS_native_combo'][:, :, :]  # time step, level, lat, lon
    
    lat_0=nc_t.variables['lat'][:]
    lon_0=nc_t.variables['lon'][:]


    for ii in range(0,COS_native_combo.shape[0]):  
        fig=plt.figure(figsize=(19.2, 10.8))
        ax = fig.add_axes([0.1,0.1,0.8,0.8])
        #width = 7500000; lon_0 = -105; lat_0 = 45
        
        m = Basemap(projection='laea',lon_0=-105,lat_0=45.,\
            llcrnrlat=9,urcrnrlat=58,\
            llcrnrlon=-135,urcrnrlon=-28,\
            rsphere=6371200.,resolution='l',area_thresh=10000)
        
        # m = Basemap(width=width,height=width,projection='aeqd',
        #     lat_0=lat_0,lon_0=lon_0)
        # draw coastlines, state and country boundaries, edge of map.
        m.drawcoastlines()
        m.drawstates()
        m.drawcountries()
        # draw parallels and meridians.
        # label parallels on right and top
        # meridians on bottom and left
        parallels = np.arange(0.,81,10.)
        # labels = [left,right,top,bottom]
        m.drawparallels(parallels,labels=[False,True,True,False])
        meridians = np.arange(10.,351.,20.)
        m.drawmeridians(meridians,labels=[True,False,False,True])
        ny = COS_native_combo.shape[1]; nx = COS_native_combo.shape[2]
        
        #ny = data.shape[0]; nx = data.shape[1]
        lons, lats = m.makegrid(nx, ny) # get lat/lons of ny by nx evenly space grid.
        x, y = m(lon_0, lat_0) # compute map proj coordinates.
        # draw filled contours.
        levels = np.linspace(-15, 15, 10)
        cs = m.contourf(x,y,COS_native_combo[ii,:,:],levels=levels,cmap='RdYlGn_r',extend="both")
        cs.cmap.set_under(cmap_vals(0))
        cs.cmap.set_over(cmap_vals(.99))
        cs.set_clim(-15, 15)
        m.drawmapboundary(fill_color='white')
        
        fontname = 'Open Sans'
        fontsize = 14
        
        # Positions for the label and time counter
        date_x = -134
        date_y = 13
        date_spacing = 65
        # Date text
        xpt,ypt = m(date_x,date_y)
        fig.tight_layout(pad=0.5)
    
        plt.text(xpt+100000,ypt+100000,"Year: {0}\nMonth: {1}\nDay Counter : {2:.3f}".format(i,j,np.divide((ii+1),8.0)),color='black',
                fontname=fontname, fontsize=fontsize*1.3,)
        cbar = m.colorbar(cs,location='bottom',pad="5%",)
        fig.savefig("/Users/yshiga/Documents/Data/soil_ocs/Jul/frames_{0}.png".format(ii),dpi=100,frameon=False,facecolor='white')
        ax.clear
        plt.close()