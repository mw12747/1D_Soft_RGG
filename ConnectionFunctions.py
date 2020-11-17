import matplotlib.pyplot as plt
import numpy as np

def ray(r_c, r, eta):
    return np.exp(-(r / r_c) ** (eta))

def rgg(r):
    return r <= 1

mu = 1
x_list = np.arange(0, 6, 0.001)

waxman = [ray(mu, x, 1) for x in x_list]
rayleigh = [ray(mu, x, 2) for x in x_list]
example = [ray(mu, x, 10) for x in x_list]
unit = [rgg(x) for x in x_list]

plt.figure()
plt.rcParams.update({'font.size': 12})

# plt.title('Different Connection Functions')

plt.plot(x_list, waxman, ':', color='darkgrey', label='$\eta$ = 1 (Waxman)')
plt.plot(x_list, rayleigh, '-.', color='grey', label='$\eta$ = 2 (Rayleigh)')
plt.plot(x_list, example, '--', color='dimgrey', label='$\eta$ = 10')
plt.plot(x_list, unit, '-', color='k', label='$\eta$ = $\infty$ (Hard RGG)')

plt.xlabel("$r / r_c$")
plt.ylabel("$H(r)$")

plt.xlim(-0.05)
plt.ylim(-0.02)

plt.legend()

plt.grid(color='whitesmoke', which='major', linestyle='-', linewidth='1')
plt.minorticks_on()
plt.grid(color='whitesmoke', which='minor', linestyle='--', linewidth='1')


plt.show()
