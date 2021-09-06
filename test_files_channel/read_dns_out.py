import numpy as np
# import sys
import matplotlib.pyplot as plt
import os
import my_pylib as mp
# from scipy import integrate
# import matplotlib.colors as mcolors
#---------------------------------------------------------------------------#
# path to dns.out file
path = str(os.path.dirname(__file__) + '/../test_yamo_180/' )

# plot settings
plt.rcParams['figure.dpi'] = 250 
size    = (8,6)
shading = 'nearest'#'gouraud'
figs    = 'figs' 
plt.close('all')
#---------------------------------------------------------------------------#
# read grid 
grid = mp.DnsGrid(path+'grid')
#---------------------------------------------------------------------------#
# read dns.out
with open(path+'dns.out', 'r') as f:
    dns_out = f.readlines()
    dns_out = np.asanyarray(dns_out)
    #---------------------------------------------#
    # get part of header
    dns_header = dns_out[1]
    # clean up
    dns_out_clean = dns_out
    index_clean   = []
    for i in range(dns_out.size):
        if dns_out[i]==dns_header:
            index_clean.append(i)
    j = 0
    for i in index_clean:
        dns_out_clean = np.delete(dns_out_clean,range(i-1-j,i+3-j), axis=0)
        j = j + 4
    #---------------------------------------------#
    itera      = np.zeros(dns_out_clean.size)
    time_tot   = np.zeros(dns_out_clean.size)
    time_dt    = np.zeros(dns_out_clean.size)
    cfl_number = np.zeros(dns_out_clean.size)
    dif_number = np.zeros(dns_out_clean.size)
    visc       = np.zeros(dns_out_clean.size)
    dil_min    = np.zeros(dns_out_clean.size)
    dil_max    = np.zeros(dns_out_clean.size)
    ubulk      = np.zeros(dns_out_clean.size)
    force      = np.zeros(dns_out_clean.size)
    #---------------------------------------------#
    i = 0
    for line in dns_out_clean:
        data     = line.split()
        #-----------------------------------------#
        itera[i]      = int(data[1])
        time_tot[i]   = float(data[2])
        time_dt[i]    = float(data[3])
        cfl_number[i] = float(data[4])
        dif_number[i] = float(data[5])
        visc[i]       = float(data[6])
        dil_min[i]    = float(data[7])
        dil_max[i]    = float(data[8])
        ubulk[i]      = float(data[9])
        force[i]      = float(data[10])
        i += 1
del dns_out, dns_out_clean, data, f, line
#---------------------------------------------------------------------------#
# plots

# plot ubulk - iteration steps
plt.figure(figsize=size)
plt.title('streamwise bulk velocity')
plt.xlabel("iterations")
plt.ylabel("ubulk")
plt.ylim(0.6,0.7)
plt.plot(itera, ubulk,label="ubulk")
plt.legend(loc=1)
plt.grid(True)
plt.show()

# # plot ubulk - time
# plt.figure(figsize=size)
# plt.title('streamwise bulk velocity')
# plt.xlabel("time_tot")
# plt.ylabel("ubulk")
# plt.ylim(0.6,0.7)
# plt.plot(time_tot, ubulk,label="ubulk")
# plt.legend(loc=1)
# plt.grid(True)
# plt.show()

# # plot time_dt
# plt.figure(figsize=size)
# plt.title(' ')
# plt.xlabel("iterations")
# plt.ylabel("dt")
# plt.ylim(0.01,0.025)
# plt.plot(itera, time_dt,label="dt")
# plt.legend(loc=1)
# plt.grid(True)
# plt.show()