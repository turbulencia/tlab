#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import netCDF4 as nc
import sys
import numpy as np

if __name__ == "__main__":
    B0 = 0.005       # Surface buoyancy flux in simulation units.
    N0 = 3. **0.5    # Buoyancy frequency in the free atmosphere, in simulation units.
    L0 = (B0/N0**3.)**0.5
    #
    sdata = nc.Dataset(sys.argv[1], 'r')
    sdata.set_auto_mask(False)
    it = sdata.variables['it'][:]
    t  = sdata.variables['t'][:]
    y  = sdata.variables['y'][:]
    b  = sdata.variables['rS'][:,:]
    b_bg  = N0**2.0 *y
    z_enc = np.sqrt( np.trapz(b[:,:]-b_bg[None,:], y, axis=1) *2.0 /N0**2.0 )
    #
    for i in range(np.size(it)):
        print(it[i],t[i],z_enc[i]/L0)
