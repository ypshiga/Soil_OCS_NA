# -*- coding: utf-8 -*-
"""
Created on Fri Aug 10 14:52:19 2018
Updated May 4 2020
@author: Yoichi Shiga
"""
# Soil OCS fluxes over North America
# using equations from Whelan et al 2016
# and NARR reanalysis soil temperature and soil moisture (top subsurface level)
# Creates soil flux maps for Forests, Grasslands, Croplands, and Wetlands
# Saves NARR resolution (.3 degree ~32km) soil OCS fluxes combining land types
# based on NARR native land cover type map which is based on SiB

import numpy as np
import os.path
from netCDF4 import Dataset
import pygrib


###############################################
##### Soil COS Emission to Atmosphere###########
##############################################
# The exponential shape of this curve
# was first noted by Liu et al 2010
# and observed to death by Maseyk et al 2014

# These fits are from lab-based soil incubation experiments
# For forest and for US agricultural soil
# Formula cos_flux = Ae^(B*Temperature_in_C)

# np.exp(regout.params[0]+regout.params[1]*graph_temp)

# r2 = 0.99 for WRC, Willow Creek Forest Tall Tower Site
wrc_A = 0.073723303
wrc_B = np.log(1.133994086)


# r2 = 0.86 for BOND, Bondville Fluxnet site, a soy/corn field
bond_A = 0.436558797
bond_B = np.log(1.103409334)

# r2 = 0.87 for STUNT, Stunt Ranch UC Rerserve, CA, savannah
stunt_0 = -4.26448842
stunt_1 = 0.1084323


def COS_ag_production(soil_temp):
    return bond_A*np.exp(bond_B*soil_temp)


def COS_forest_production(soil_temp):
    return wrc_A*np.exp(wrc_B*soil_temp)  # in pmol m^-2 sec^-1


def COS_grass_production(soil_temp):
    return np.exp(stunt_0+stunt_1*soil_temp)  # in pmol m^-2 sec^-1

###############################################
# Soil COS Uptake from Atmosphere###########
##############################################
# This is an NO production model borrowed from Behrendt et al 2014
# Repurposed for COS consumption
# Van Diest and Kesselmeier 2008 re-purposed a different NO production model,
# but the spirit of the pursuit is the same

# Calculate max COS uptake, soilw at max uptake, and uptake at 35% VWC
# with soil temperature (in C) as input
# 35 was chosen arbitrarily to reduce the number of variables when curve
# fitting
# and because the soil moisture at max uptake is generally <40

# Things blow up at 0 C, so this won't support freezing temperatures
# Roisin Commane suggests we can approximate fluxes at 0 C with
# 0 COS pmol m^-2 sec^-1

# For each temperature (bin or individual value), this calculates the
# soil moisture response curve
# The first equation is a curve "shape" value, a


def curve_shape_a(opt_flux, flux_at_theta_g, theta_g, theta_opt):
    Rj = opt_flux/flux_at_theta_g
    outs = np.log(Rj)
    inversed = (np.log(theta_opt/theta_g)) + (theta_g/theta_opt - 1)
    ins = 1/inversed
    return outs*ins


def flux_at_all_theta(all_theta, opt_flux, flux_at_theta_g, theta_g, theta_opt):
    a = curve_shape_a(opt_flux, flux_at_theta_g, theta_g, theta_opt)
    power_fct_increasing = (all_theta/theta_opt) ** a
    exp_fct_decreasing = np.exp(-1 * a * ((all_theta/theta_opt) - 1))
    return opt_flux * power_fct_increasing * exp_fct_decreasing


def flux_theta_g_constant(all_theta, opt_flux, flux_at_theta_g, theta_g, theta_opt):
    a = curve_shape_a(opt_flux, flux_at_theta_g, 35, theta_opt)
    power_fct_increasing = (all_theta/theta_opt) ** a
    exp_fct_decreasing = np.exp(-1 * a * ((all_theta/theta_opt) - 1))
    return opt_flux * power_fct_increasing * exp_fct_decreasing

################################################
###### empirically derived coefficients##########
####### for uptake over varying temperatures####
##############################################
# These equations were derived for COS fluxes from
# multiple subsamples of a soil
# From the Bondville Fluxnet site
# with the production component (COS_ag_production function above)
# subtracted out

# for Optimum (maximum) uptake at some temperature
# polynomial fit
# r^2 = 0.98


opt_flux_a = -0.009855866
opt_flux_b = 0.197187006
opt_flux_c = -9.323708203
# SW at which Optimum COS uptake occurs
# r^2 = 0.68

opt_sw_m = 0.287408209
opt_sw_b = 14.50082209
# COS uptake at 35 VWC%
# r^2 = 0.87


other_flux_a = -0.011921136
other_flux_b = 0.110444854
other_flux_c = -1.178897129

# arbitrarily chosen, a soil moisture > optimum soil moisture
other_sw = 35


def opti_flux(temperature):
    return temperature**2 * opt_flux_a + temperature*opt_flux_b + opt_flux_c


def other_flux(temperature):
    return temperature**2 * other_flux_a + temperature * other_flux_b + other_flux_c


def sw_at_opt_flux(temperature):
    return temperature * opt_sw_m + opt_sw_b


########################################################
##### Soil COS Uptake and Emission Together#############
#######################################################
# Now, add these two together
# For each temperature bin/individual temperature
    
# ** Modified to exclude temperature filter (I now do this after caclulating
# all fluxes) ***
def COS_over_T_and_SW_ag(temperature, soilw):
#    if temperature <=0:
#        return 0, 0
#    else:
        opt_flux = opti_flux(temperature)
        flux_at_theta_g = other_flux(temperature)
        theta_g = other_sw
        theta_opt = sw_at_opt_flux(temperature)
        biotic_flux = flux_at_all_theta(soilw, opt_flux, flux_at_theta_g, theta_g, theta_opt)
        abiotic_flux = COS_ag_production(temperature)
        # return biotic_flux + abiotic_flux
        return biotic_flux, biotic_flux + abiotic_flux


def COS_over_T_and_SW_forest(temperature, soilw):
#    if temperature <=0:
#        return 0, 0
#    else:
        opt_flux = opti_flux(temperature)
        flux_at_theta_g = other_flux(temperature)
        theta_g = other_sw
        theta_opt = sw_at_opt_flux(temperature)
        biotic_flux = flux_at_all_theta(soilw, opt_flux, flux_at_theta_g, theta_g, theta_opt)
        abiotic_flux = COS_forest_production(temperature)
        #return biotic_flux + abiotic_flux
        return biotic_flux, biotic_flux + abiotic_flux

    
def COS_over_T_and_SW_grass(temperature, soilw):
#    if temperature <=0:
#        return 0, 0
#    else:
        opt_flux = opti_flux(temperature)
        flux_at_theta_g = other_flux(temperature)
        theta_g = other_sw
        theta_opt = sw_at_opt_flux(temperature)
        biotic_flux = flux_at_all_theta(soilw, opt_flux, flux_at_theta_g, theta_g, theta_opt)
        abiotic_flux = COS_grass_production(temperature)
        # return biotic_flux + abiotic_flux
        return biotic_flux, biotic_flux + abiotic_flux
    
def wetland(t):
    #at 0, should be 0
    #at 35 should be 120
    # I now change to t=0 to zero after 
#    if t <= 0:
#        return 0
#    else:
    return np.exp(-1.46700192)*np.exp(0.18836697*t)

# Directory for data
data_dir = "/Users/yshiga/Documents/Data/NARR/"
# get vegetation type mask
# NARR saves this out as a grb file
veg_type_name = data_dir + "rr-fixed.grb"
grb_file = pygrib.open(veg_type_name)
veg_type_ind = []
for grb in grb_file:
    if grb.parameterName == "225" : 
        veg_type_ind.append(grb.values)
veg_type_ind = np.array(veg_type_ind)

# Veg type in NARR from SiB       
# 0	Water Bodies
# 1	Evergreen Broadleaf Trees 
# 2	Broadleaf Deciduous Trees
# 3	Deciduous and Evergreen Trees
# 4	Evergreen Needleleaf Trees
# 5	Deciduous Needleleaf Trees
# 6	Ground Cover with Trees and Shrubs
# 7	Groundcover Only
# 8	Broadleaf Shrubs with Perennial Ground Cover
# 9	Broadleaf Shrubs with Bare Soil
# 10	Groundcover with Dwarf Trees and Shrubs
# 11	Bare Soil
# 12	Agriculture or C3 Grassland
# 13	Persistent Wetland
# 14	Ice Cap and Glacier
# 15	Missing Data
     
# create veg type masks   
# Forest = 1,2,3,4,5,6,10
forest_ind=np.array((veg_type_ind[0,:,:]==1, veg_type_ind[0,:,:]==2, veg_type_ind[0,:,:]==3, veg_type_ind[0,:,:]==4, veg_type_ind[0,:,:]==5, veg_type_ind[0,:,:]==6, veg_type_ind[0,:,:]==10))
forest_ind=np.logical_or.reduce(forest_ind)
forest_ind= forest_ind[np.newaxis,np.newaxis,:,:]
# Grasslands = 7,8, 9,11
grass_ind=np.array((veg_type_ind[0,:,:]==7, veg_type_ind[0,:,:]==8, veg_type_ind[0,:,:]==9, veg_type_ind[0,:,:]==11))
grass_ind=np.logical_or.reduce(grass_ind)
grass_ind= grass_ind[np.newaxis,np.newaxis,:,:]
# Croplands = 12
crop_ind=np.array((veg_type_ind[0,:,:]==12))
crop_ind= crop_ind[np.newaxis,np.newaxis,:,:]
# Wetlands = 13
wetland_ind=np.array((veg_type_ind[0,:,:]==13))
wetland_ind= wetland_ind[np.newaxis,np.newaxis,:,:]
# Zero = 14, 0, 15
zero_ind=np.array((veg_type_ind[0,:,:]==0, veg_type_ind[0,:,:]==14, veg_type_ind[0,:,:]==15))
zero_ind=np.logical_or.reduce(zero_ind)
zero_ind= zero_ind[np.newaxis,np.newaxis,:,:]

for i in range(2006,2007):  # loop through years 2006 - 2017
    for j in range(1,13):  # loop through months 1-12
        # Create netcdf file name
        filename_ocs = data_dir + "soil_ocs/soil_ocs_combo.{0}{1:02d}.nc".format(i, j)
        if os.path.isfile(filename_ocs):
             print "File for {0}{1:02d} already exists".format(i,j)
        else:
            # Load temperature and soil moisture data
            filename_t = data_dir + "tsoil/tsoil.{0}{1:02d}.nc".format(i, j)
            filename_s = data_dir + "soill/soill.{0}{1:02d}.nc".format(i, j)
            print(i)
            print(j)
            nc_t = Dataset(filename_t, 'r')  # open soil temp for given month
            soil_temp = nc_t.variables['tsoil'][:, 0, :, :]  # time step, level, lat, lon
            soil_temp = soil_temp - 273.15  # convert to Celsius
            nc_s = Dataset(filename_s, 'r')  # open soil moisture for given month
            soilw = np.multiply(nc_s.variables['soill'][:, 0, :, :],100)  # time step, level, lat, lon
            ind_0 = soil_temp <= 0  # find where soil temp is less than or equal to zero
            ind_hot = soil_temp > 45  # find where soil temp above 45
            N = ind_0[np.newaxis,:,:,:] # create index for zero temp
            N = np.tile(N,(1,1,1)) # tile to match dimensiosn
            N_hot = ind_hot[np.newaxis,:,:,:] # create index for temp above 45
            N_hot = np.tile(N_hot,(1,1,1)) # tile to match dimension
            soil_temp=np.where(N_hot,45,soil_temp) # set max soil temp to 45
            # Calculate soil OCS fluxes for Ag, forest, grasslands, and wetlands
            COS_ag = COS_over_T_and_SW_ag(soil_temp, soilw)
            COS_forest = COS_over_T_and_SW_forest(soil_temp, soilw)
            COS_grass = COS_over_T_and_SW_grass(soil_temp, soilw)
            COS_wet =  wetland(soil_temp)
            # If temp is less than zero set soil OCS fluxes to zero
            COS_ag_no_zero = np.where(N, 0, COS_ag[1])
            COS_forest_no_zero = np.where(N, 0, COS_forest[1])
            COS_grass_no_zero = np.where(N, 0, COS_grass[1])
            COS_uptake_no_zero = np.where(N, 0, COS_grass[0]) # all "uptake" use same equation 
            COS_wet_no_zero = np.where(N, 0, COS_wet[0]) # wetlands
            # If it is not land set soil OCS fluxes to nan
            land_var_name = data_dir + "land.nc"
            nc_land = Dataset(land_var_name, 'r')
            land_mask = nc_land.variables['land'][0, :, :]  # land mask
            land_ind = land_mask == 1 #create index for land
            N1 = land_ind[np.newaxis,np.newaxis,:,:] #create index for land adding dimensions
            dim_one=len(COS_ag_no_zero[0]) # number of time steps in month
            N1 = np.tile(N1,(len(COS_ag_no_zero),dim_one,1,1)) # tile to match dimensions
            # remove land - set to NaN
            COS_ag_land = np.where(N1 == False, np.nan, COS_ag_no_zero)
            COS_forest_land = np.where(N1 == False, np.nan, COS_forest_no_zero)
            COS_grass_land = np.where(N1 == False, np.nan, COS_grass_no_zero)
            COS_uptake_land = np.where(N1 == False, np.nan, COS_uptake_no_zero)
            COS_wet_land = np.where(N1 == False, np.nan, COS_wet_no_zero)
            #subset based on vegetation type
            # index for veg types
            forest_ind_long = np.tile(forest_ind,(dim_one,1,1))
            zero_ind_long = np.tile(zero_ind,(dim_one,1,1))
            wetland_ind_long = np.tile(wetland_ind,(dim_one,1,1))
            crop_ind_long = np.tile(crop_ind,(dim_one,1,1))
            grass_ind_long = np.tile(grass_ind,(dim_one,1,1))
            # start with forest
            COS_native_combo = np.ma.copy(COS_forest_land)
            # replace grasslands
            COS_native_combo[grass_ind_long] = COS_grass_land[grass_ind_long]
            # replace croplands
            COS_native_combo[crop_ind_long] = COS_ag_land[crop_ind_long]
            # replace wetlands
            COS_native_combo[wetland_ind_long] = COS_wet_land[wetland_ind_long]
            # replace other "zero"
            COS_native_combo[zero_ind_long] = 0
            # remove non land
            COS_native_combo = np.where(N1 == False, np.nan, COS_native_combo)
            # variables to exclude when saving Netcdf
            toexclude_dim = ["level"]
            toexclude_var = ['level', 'tsoil','Lambert_Conformal']
            #save netcdf
            with Dataset.Dataset(filename_t) as src, Dataset.Dataset(filename_ocs, "w") as dst:
                # copy attributes
                for name in src.ncattrs():
                    dst.setncattr(name, src.getncattr(name))  # copy global attributes all at once via dictionary
                # copy dimensions
                for name, dimension in src.dimensions.items():
                    # if name not in toexclude_dim:
                        dst.createDimension(
                                name, (len(dimension) if not dimension.isunlimited() else None))
                # copy all file data except for the excluded
                for name, variable in src.variables.iteritems():
                    if name not in toexclude_var:
                       #  print name
                        x = dst.createVariable(name, variable.datatype, variable.dimensions)
                        dst.variables[name][:] = src.variables[name][:]
                COS_native_combo = COS_native_combo[0,:,:]
                dst.createVariable("COS_native_combo", "f4", (u'time', u'y', u'x'), zlib=True)
                dst.variables['COS_native_combo'][:] = COS_native_combo

