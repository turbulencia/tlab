import numpy as np
import scipy.integrate as integrate
import matplotlib.pyplot as plt
from template import *

h = 6.626e-34           # Planck constant, in J s
k = 1.380e-23           # Boltzmann constant, in J /K
c = 2.998e8             # speed of light, in m /s

c2 = h*c/k

# Check: Should recover Stefan-Boltzmann
result = integrate.quad(lambda x: x**3./(np.exp(x)-1.0), 1.0e-2, 100.)
# print(result)
B = (15./np.pi**4.)*result[0]
print(B)

Ts = np.linspace(250., 320., num=32)   # temperatures

# band 1
B1 = np.empty_like(Ts)
nu = (15000.-5600., 15000.+5600.)       # frequency interval, in 1 /m
for idx, T in enumerate(Ts):
    x1 = c2*nu[0]/T
    x2 = c2*nu[1]/T
    result = integrate.quad(lambda x: x**3./(np.exp(x)-1.0), x1, x2)
    B1[idx] = (15./np.pi**4.)*result[0]

# band 2
B2 = np.empty_like(Ts)
nu = (150000.-4000., 150000.+4000.)     # frequency interval, in 1 /m
for idx, T in enumerate(Ts):
    x1 = c2*nu[0]/T
    x2 = c2*nu[1]/T
    result = integrate.quad(lambda x: x**3./(np.exp(x)-1.0), x1, x2)
    B2[idx] = (15./np.pi**4.)*result[0]

coef1=np.polyfit(Ts, B1, 2)
print(coef1)
B1_fit = coef1[0]*Ts**2 +coef1[1]*Ts +coef1[2]
print('Inf-error {}'.format(np.amax(np.abs(B1-B1_fit)/B1)))

coef2=np.polyfit(Ts, B2, 2)
B2_fit = coef2[0]*Ts**2+coef2[1]*Ts+coef2[2]
print(coef2)
print('Inf-error {}'.format(np.amax(np.abs(B2-B2_fit)/B2)))

fig, axs = plt.subplots( nrows=1, ncols=2, figsize = (8, 3) )
axs[0].plot(Ts, B1, 'o')
axs[0].plot(Ts, B1_fit, label=r'fit')
axs[0].set_xlabel(r'$T$ ($K$)')
axs[0].set_ylabel(r'$\beta_\mathrm{b,1}(T)$')
axs[0].legend(loc='best')

axs[1].plot(Ts, B2, 'o')
axs[1].plot(Ts, B2_fit, label=r'fit')
axs[1].set_xlabel(r'$T$ ($K$)')
axs[1].set_ylabel(r'$\beta_\mathrm{b,2}(T)$')

plt.tight_layout(pad=0.1)
plt.savefig("radiation1.pdf",bbox_inches='tight')

plt.show()
