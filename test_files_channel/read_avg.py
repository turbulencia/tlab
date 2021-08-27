import numpy as np
import sys
import matplotlib.pyplot as plt
import os
import my_pylib as mp
import netCDF4  as nc 
import warnings; warnings.filterwarnings("ignore", category=DeprecationWarning) 

#-----------------------------------------------------------------------------#
# path to avg-fields
# path    = str(os.path.dirname(__file__) + '/../test_little_channel/' )
path    = str(os.path.dirname(__file__) + '/../test_yamo_180/' )
name    = 'avg'

# index
index_nc   = 100000
index_flow = 100000
file_nc    = str(path+name+str(index_nc)+'.nc')

# plot settings 
plt.rcParams['figure.dpi'] = 250 
size    = (8,6)
shading = 'gouraud'#'gouraud'
figs    = 'figs' 
plt.close('all')

#-----------------------------------------------------------------------------#
# read grid
grid = mp.DnsGrid(path+'grid')

# read flow fields
flow = mp.Field(path,var='flow',index=index_flow)
flow.read_3d_field()

# read avg.nc netcdf field
avg = nc.Dataset(file_nc, 'r')
# avg.variables.keys() # see all variables stored
#-----------------------------------------------------------------------------#
# simulation propteries

# forcing type
forcing_cpg = True  # constant pressure gradient
forcing_cfr = False # constant flow rate

re = 180            # reynoldsnumber in dns.ini file
if forcing_cpg:
    re_tau = re
    re_cl  = (re_tau/0.116)**(1/0.88)
    fcpg   = ( re_tau / re_cl )**2    # constant streamwise pressure gradient
if forcing_cfr:
    re_cl  = re 
    re_tau = re_cl**0.88 * 0.116  
    
nu    = 1 / re_cl  # kinematic viscosity, always!
# ro    = 1          # Rossby number
# pr    = 1          # Prandtl number
rho   = 1          # density
delta = 1          # channel half height
#-----------------------------------------------------------------------------#
# load variables

# time and vertical distance
t=avg.variables['t'][:]
y=avg.variables['y'][:]

# velocities
u=avg.variables['rU'][:,:] # the first index is time, the second is vertical node
v=avg.variables['rV'][:,:]
w=avg.variables['rW'][:,:]

# Re-stresses
# uv=avg.variables['Rxy'][:,:]

# pressure
p=avg.variables['rP'][:,:]

#-----------------------------------------------------------------------------#
# balance
# dudy2 = mp.derivative2(grid.y, u.mean(axis=0), grid.ny)
# duvdy = mp.derivative1(grid.y, uv.mean(axis=0), grid.ny)
# dpdx  = nu * dudy2 - duvdy[:,0]

#-----------------------------------------------------------------------------#
# plots - averages

# prabolic u_mean ini profile
ucl   = 1                                            # centerline velocity 
ycl   = grid.y.max() / 2                             # centerline position
u_par = - (ucl / ycl**2 ) * (grid.y - ycl)**2 + ucl  # parabolic ini velocity profile

# velocities
plt.figure(figsize=size)
plt.xlabel("y")
plt.ylabel("u_i")
plt.plot(y, u_par, '.', label='u_mean_ini')
plt.plot(y, u.mean(axis=0), label='u_mean')
plt.plot(y, v.mean(axis=0), label='v_mean')
# plt.plot(y, w.mean(axis=0), label='w_mean')
# plt.plot(y, p.mean(axis=0), label='p_mean')
plt.legend()
# plt.xscale('log')
# plt.title('horizontal - velocities')
plt.grid(True)
# plt.xlim(1,10000)
# plt.ylim(0,25)
plt.show()
# plt.savefig("plots_neutral/velocities.svg")

#-----------------------------------------------------------------------------#
# plots - 2d flow fields

plt.figure(figsize=size)
plt.title('2d-plot -- xy-plane -- velocity u')
plt.xlabel("x")
plt.ylabel("y")
plt.pcolormesh(grid.x,grid.y,flow.u[:,:,0].T, shading=shading ,cmap='RdBu_r')#, norm=midnorm)
plt.colorbar()
plt.show()
#
plt.figure(figsize=size)
plt.title('2d-plot -- xy-plane -- velocity v')
plt.xlabel("x")
plt.ylabel("y")
plt.pcolormesh(grid.x,grid.y,flow.v[:,:,0].T, shading=shading ,cmap='RdBu_r')#, norm=midnorm)
plt.colorbar()
plt.show()
#
# plt.figure(figsize=size)
# plt.title('2d-plot -- xy-plane -- velocity w')
# plt.xlabel("x")
# plt.ylabel("y")
# plt.pcolormesh(grid.x,grid.y,flow.w[:,:,0].T, shading=shading ,cmap='RdBu_r')#, norm=midnorm)
# plt.colorbar()
# plt.show()

#-----------------------------------------------------------------------------#
sys.exit()
#-----------------------------------------------------------------------------#

# velocities
plt.figure(figsize=size)
plt.xlabel("y")
plt.ylabel("terms")
plt.plot(y, nu*dudy2, label='nu*dudy2')
plt.plot(y, -duvdy,   label='-duvdy')
plt.plot(y, -dpdx,    label='-dpdx')
plt.legend()
plt.grid(True)
# plt.xlim(1,10000)
# plt.ylim(0,25)
plt.show()
# plt.savefig("plots_neutral/velocities.svg")


'''
# # Re-stresses / vel. fluc.
# uu=datafile_nc.variables['Rxx'][:,:] # the first index is time, the second is vertical node
# uv=datafile_nc.variables['Rxy'][:,:] # the first index is time, the second is vertical node
# uw=datafile_nc.variables['Rxz'][:,:] # the first index is time, the second is vertical node
# vv=datafile_nc.variables['Ryy'][:,:] # the first index is time, the second is vertical node
# ww=datafile_nc.variables['Rzz'][:,:] # the first index is time, the second is vertical node
# vw=datafile_nc.variables['Ryz'][:,:] # the first index is time, the second is vertical node
# k = 0.5 * (uu + vv + ww)             # turbulent kinetic energy

# # frction velocity
# u_tau=datafile_nc.variables['FrictionVelocity'][:] # the first index is time
# u_tau_mean    = np.zeros(len(t))
# u_tau_mean[:] = u_tau.mean()

# # friction angle 
# alpha=datafile_nc.variables['FrictionAngle'][:] # the first index is time
# alpha_mean    = np.zeros(len(t))
# alpha_mean[:] = alpha.mean()

#-----------------------------------------------------------------------------#
# compute
#-----------------------------------------------------------------------------#

# stepsizes
delta_y = y[1:] - y[:-1]
delta_t = t[1:] - t[:-1]

# wallunits with Ro=f=G=1 (default values)
delta_v = nu / u_tau_mean[0]     # viscous length scale
re_tau  = u_tau_mean[0]**2 / nu
y_plus  = y / delta_v

# viscous sublayer
y_vsl = np.linspace(0,10,1001)
u_vsl = y_vsl

# log-law
kappa = 0.405
b = 5
y_log = np.linspace(10,1000,1001)
u_log = kappa**-1 * np.log(y_log) + b 

## compute u_tau with central difference scheme

# first velocity derivations in y
dudy    = np.zeros(y.size) # initialise
dwdy    = np.zeros(y.size)                                
dudy[0] = (u.mean(axis=0)[1] - u.mean(axis=0)[0]) / (y[1] - y[0]) # lower wall
dwdy[0] = (w.mean(axis=0)[1] - w.mean(axis=0)[0]) / (y[1] - y[0])          
dudy[y.size-1] = (u.mean(axis=0)[y.size-1] - u.mean(axis=0)[y.size-2]) / (y[y.size-1]  - y[y.size-2]) # upper wall
dwdy[y.size-1] = (w.mean(axis=0)[y.size-1] - w.mean(axis=0)[y.size-2]) / (y[y.size-1]  - y[y.size-2]) 
for i in range(1, y.size-1):
    dudy[i] = (u.mean(axis=0)[i+1] - u.mean(axis=0)[i-1]) / (y[i+1] - y[i-1]) # in between 
    dwdy[i] = (w.mean(axis=0)[i+1] - w.mean(axis=0)[i-1]) / (y[i+1] - y[i-1]) 

# tau
tau_xy_vis = nu*dudy
tau_yz_vis = nu*dwdy
tau_vis    = np.sqrt(tau_xy_vis**2 + tau_yz_vis**2) 

tau_xy_turb = - uv.mean(axis=0)
tau_yz_turb = - vw.mean(axis=0)
tau_turb    = np.sqrt(tau_xy_turb**2 + tau_yz_turb**2)

tau_xy  = tau_xy_vis + tau_xy_turb
tau_yz  = tau_yz_vis + tau_yz_turb
tau_tot = np.sqrt(tau_xy**2 + tau_yz**2) 

tau_w              = tau_tot[0]
u_tau_comp_mean    = np.zeros(len(t))
u_tau_comp_mean[:] = np.sqrt(tau_w / rho)

# u_tau_mean tlab
u_tau_mean_tlab    = np.zeros(len(y))
u_tau_mean_tlab[:] = u_tau.mean()

# friction angle
alpha_com = np.degrees(np.arctan(tau_yz[0] / tau_xy[0]))
alpha_com_mean    = np.zeros(len(t))
alpha_com_mean[:] = alpha_com.mean()

#sys.exit()
#-----------------------------------------------------------------------------#
# plot variables
#-----------------------------------------------------------------------------#
plt.close("all")
# plot settings 
plt.rcParams['figure.dpi'] = 250 
size = (8,6)
#sys.exit()
#-----------------------------------------------------------------------------#

# grid and time stepping
plt.figure(figsize=size)
plt.xlabel("steps")
plt.ylabel("delta_y+")
plt.plot(delta_y / delta_v, label='delta_y+')
#plt.plot(w.mean(axis=0), y, label='w_mean-velocity')
plt.legend()
plt.title('grid - vertical step size (delta_y+)')
plt.grid(True)
plt.xlim(0,1200)
plt.ylim(0,25)
plt.show()
plt.savefig("plots_neutral/delta_y+.svg")

plt.figure(figsize=size)
plt.xlabel("steps")
plt.ylabel("delta_t")
plt.plot(delta_t, label='delta_t')
#plt.plot(w.mean(axis=0), y, label='w_mean-velocity')
plt.legend()
plt.title('time steps (delta_t)')
plt.grid(True)
plt.xlim(0,60)
plt.ylim(0.003,0.0034)
plt.show()
plt.savefig("plots_neutral/delta_t.svg")

#-----------------------------------------------------------------------------#

# velocities
plt.figure(figsize=size)
plt.xlabel("y+")
plt.ylabel("u_i+")
plt.plot(y_plus, u.mean(axis=0)/u_tau[0], label='u_mean')
plt.plot(y_plus, w.mean(axis=0)/u_tau[0], label='w_mean')
plt.plot(y_plus, np.sqrt(u.mean(axis=0)**2 + w.mean(axis=0)**2)/u_tau[0], label='sqrt(u**2 + w**2)')
plt.plot(y_vsl, u_vsl, linestyle = 'dotted', label='visc. sublayer')
plt.plot(y_log, u_log, linestyle = 'dotted', label='log layer')
plt.legend()
plt.xscale('log')
plt.title('horizontal - velocities')
plt.grid(True)
plt.xlim(1,10000)
plt.ylim(0,25)
plt.show()
plt.savefig("plots_neutral/velocities.svg")

#-----------------------------------------------------------------------------#

# Reymolds stresses
plt.figure(figsize=size)
plt.xlabel("y+")
plt.ylabel("u'u'+")
plt.plot(y_plus, uu.mean(axis=0)/u_tau[0]**2, label='uu_mean')
plt.plot(y_plus, vv.mean(axis=0)/u_tau[0]**2, label='vv_mean')
plt.plot(y_plus, ww.mean(axis=0)/u_tau[0]**2, label='ww_mean')
plt.plot(y_plus, uv.mean(axis=0)/u_tau[0]**2, label='uv_mean')
plt.plot(y_plus, uw.mean(axis=0)/u_tau[0]**2, label='uw_mean')
plt.plot(y_plus, vw.mean(axis=0)/u_tau[0]**2, label='vw_mean')
plt.plot(y_plus, k.mean(axis=0) /u_tau[0]**2, label='tke_mean')
plt.legend() 
plt.xscale('log')
plt.title('Reynolds stresses')
plt.grid(True)
plt.xlim(1,10000)
plt.ylim(-1,8)
plt.show()
plt.savefig("plots_neutral/re_log.svg")

plt.figure(figsize=size)
plt.xlabel("y+")
plt.ylabel("u'u'+")
plt.plot(y_plus, uu.mean(axis=0)/u_tau[0]**2, label='uu_mean')
plt.plot(y_plus, vv.mean(axis=0)/u_tau[0]**2, label='vv_mean')
plt.plot(y_plus, ww.mean(axis=0)/u_tau[0]**2, label='ww_mean')
plt.plot(y_plus, uv.mean(axis=0)/u_tau[0]**2, label='uv_mean')
plt.plot(y_plus, uw.mean(axis=0)/u_tau[0]**2, label='uw_mean')
plt.plot(y_plus, vw.mean(axis=0)/u_tau[0]**2, label='vw_mean')
plt.plot(y_plus, k.mean(axis=0) /u_tau[0]**2, label='tke_mean')
plt.legend() 
plt.xscale('linear')
plt.title('Reynolds stresses')
plt.grid(True)
plt.xlim(0,10000)
plt.ylim(-1,8)
plt.show()
plt.savefig("plots_neutral/re_lin.svg")

#-----------------------------------------------------------------------------#

# plot tau
plt.figure(figsize=size)
plt.title('tau_ij+')
plt.xlabel("y+")
plt.ylabel("tau_ij / tau_w")
plt.plot(y_plus, tau_turb/u_tau[0]**2, linestyle='solid', color = 'orange', label='tau_tub')
plt.plot(y_plus, tau_xy_turb/u_tau[0]**2, linestyle='dotted', color = 'orange', label='tau_tub_xy')
plt.plot(y_plus, tau_yz_turb/u_tau[0]**2, linestyle='dashed', color = 'orange', label='tau_tub_yz')
plt.plot(y_plus, tau_vis/u_tau[0]**2, linestyle='solid', color = 'green', label='tau_vis')
plt.plot(y_plus, tau_xy_vis/u_tau[0]**2, linestyle='dotted', color = 'green', label='tau_vis_xy')
plt.plot(y_plus, tau_yz_vis/u_tau[0]**2, linestyle='dashed', color = 'green', label='tau_vis_yz')
plt.plot(y_plus, tau_tot/u_tau[0]**2, linestyle='solid', color = 'blue', label='tau_tot')
plt.plot(y_plus, tau_xy/u_tau[0]**2, linestyle='dotted', color = 'blue', label='tau_tot_xy')
plt.plot(y_plus, tau_yz/u_tau[0]**2, linestyle='dashed', color = 'blue', label='tau_tot_yz')
plt.legend() 
plt.grid(True)
plt.xscale('log')
plt.xlim(1,10000)
plt.ylim(0,1)
plt.show()
plt.savefig("plots_neutral/tau.svg")

# plot tau
plt.figure(figsize=size)
plt.title('tau_ij+')
plt.xlabel("y+")
plt.ylabel("tau_ij / tau_w")
plt.plot(y_plus, tau_turb/u_tau[0]**2, linestyle='solid', color = 'orange', label='tau_tub')
plt.plot(y_plus, tau_xy_turb/u_tau[0]**2, linestyle='dotted', color = 'orange', label='tau_tub_xy')
plt.plot(y_plus, tau_yz_turb/u_tau[0]**2, linestyle='dashed', color = 'orange', label='tau_tub_yz')
plt.plot(y_plus, tau_vis/u_tau[0]**2, linestyle='solid', color = 'green', label='tau_vis')
plt.plot(y_plus, tau_xy_vis/u_tau[0]**2, linestyle='dotted', color = 'green', label='tau_vis_xy')
plt.plot(y_plus, tau_yz_vis/u_tau[0]**2, linestyle='dashed', color = 'green', label='tau_vis_yz')
plt.plot(y_plus, tau_tot/u_tau[0]**2, linestyle='solid', color = 'blue', label='tau_tot')
plt.plot(y_plus, tau_xy/u_tau[0]**2, linestyle='dotted', color = 'blue', label='tau_tot_xy')
plt.plot(y_plus, tau_yz/u_tau[0]**2, linestyle='dashed', color = 'blue', label='tau_tot_yz')
plt.legend() 
plt.grid(True)
plt.xscale('linear')
plt.xlim(1,10000)
plt.ylim(0,1)
plt.show()
plt.savefig("plots_neutral/tau_linear.svg")

plt.figure(figsize=size)
plt.title('tau_ij+')
plt.xlabel("y/delta")
plt.ylabel("tau_ij / tau_w")
plt.plot(y/u_tau_mean[0], tau_turb/u_tau[0]**2, linestyle='solid', color = 'orange', label='tau_tub')
#plt.plot(y/u_tau_mean[0], tau_xy_turb/u_tau[0]**2, linestyle='dotted', color = 'orange', label='tau_tub_xy')
#plt.plot(y/u_tau_mean[0], tau_yz_turb/u_tau[0]**2, linestyle='dashed', color = 'orange', label='tau_tub_yz')
plt.plot(y/u_tau_mean[0], tau_vis/u_tau[0]**2, linestyle='solid', color = 'green', label='tau_vis')
#plt.plot(y/u_tau_mean[0], tau_xy_vis/u_tau[0]**2, linestyle='dotted', color = 'green', label='tau_vis_xy')
#plt.plot(y/u_tau_mean[0], tau_yz_vis/u_tau[0]**2, linestyle='dashed', color = 'green', label='tau_vis_yz')
plt.plot(y/u_tau_mean[0], tau_tot/u_tau[0]**2, linestyle='solid', color = 'blue', label='tau_tot')
#plt.plot(y/u_tau_mean[0], tau_xy/u_tau[0]**2, linestyle='dotted', color = 'blue', label='tau_tot_xy')
#plt.plot(y/u_tau_mean[0], tau_yz/u_tau[0]**2, linestyle='dashed', color = 'blue', label='tau_tot_yz')
plt.legend() 
plt.grid(True)
plt.xscale('linear')
plt.xlim(0,0.4)
plt.ylim(0,1)
plt.show()
plt.savefig("plots_neutral/tau_delta.svg")

#-----------------------------------------------------------------------------#

# plot hodograph
plt.figure(figsize=size)
plt.xlabel("u/G")
plt.ylabel("w/G")
plt.plot(u.mean(axis=0), - w.mean(axis=0))#, label='')
#plt.legend()
plt.title('hodograph')
plt.grid(True)
#plt.xlim(0,1)
#plt.ylim(0,0.35)
plt.show()
plt.savefig("plots_neutral/hodograph.svg")

# friction angle
plt.figure(figsize=size)
plt.xlabel("time")
plt.ylabel("alpha [deg]")
plt.plot(t[:],alpha[:], label='alpha')
plt.plot(t[:],alpha_mean[:], label='alpha_mean')
plt.plot(t[:],alpha_com_mean[:], label='alpha_computed')
plt.legend()
plt.grid(True)
plt.title('friction angle')
plt.show()
plt.savefig("plots_neutral/alpha.svg")

# friction velocity
plt.figure(figsize=size)
plt.xlabel("time")
plt.ylabel("u_tau")
plt.plot(t[:],u_tau[:], label='u_tau')
plt.plot(t[:],u_tau_mean, label='u_tau_mean')
plt.plot(t[:],u_tau_comp_mean[:], label='u_tau_computed')
plt.legend()
plt.grid(True)
plt.title('u_tau')
plt.show()
plt.savefig("plots_neutral/u_tau.svg")

#-----------------------------------------------------------------------------#
plt.close("all")

'''