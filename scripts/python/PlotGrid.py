import numpy   as np   # For array operations.
import pylab   as pl   # For plotting.
from matplotlib import rc

#pl.ion()

filenames       = ["../../5120x1024x5120-1600-20160411/ics/y.dat",
                   "../../5120x1536x5120-1600-20160411/ics/y.dat",
                   "../../5120x2048x5120-1600-20160411/ics/y.dat",
                   "../../5120x2560x5120-1600-20160411/ics/y.dat"]
reference_sizes = [10.8, 16.2, 21.6, 27.0]

pl.figure()
ax = pl.subplot(111)

i = 0
for file in filenames:
    f = np.loadtxt(file)
    n = len(f)

    reference_step = reference_sizes[i] /n

    x  = np.linspace(1, n, num=n)

    f1 = np.gradient(f)
    f2 = np.gradient(f1)

    stretching       = f2 /f1 *100         # in percentage
    normalized_delta = f1 /reference_step

    pl.plot(f,normalized_delta,'-')
    pl.plot(f,stretching,'--')

    i = i +1
    
rc('xtick', direction='in')
rc('ytick', direction='in')

pl.xlabel(r'node position')
pl.ylabel(r'spacing, stretching')
ax.grid(True)

pl.show()
