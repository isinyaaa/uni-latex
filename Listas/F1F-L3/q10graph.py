import matplotlib.pyplot as plt
import numpy as np

np.seterr(over='ignore')

with open('tabela.txt') as f:
    lines = f.readlines()[6::]
    R_Rsun_ratio = \
        np.array([line.split()[0] for line in lines])\
        .astype(np.longdouble)
    norm_electron_density = \
        np.array([line.split()[1] for line in lines])\
        .astype(np.longdouble)

A_N = 6.0221409e+23
S_R = 6.95700e+8

electron_density = 10 ** norm_electron_density * A_N * 10 ** 6
r = R_Rsun_ratio * S_R
dr = np.array([r[i + 1] - r[i] for i in range(0, len(r) - 1)])


val_index = lambda array1, val: (np.abs(array1 - val)).argmin()


def electron_func(r_sun_frac):
    return np.sum(
        [electron_density[i] * 4 * np.pi * (r[i] ** 2) * dr[i]
         for i in range(0, val_index(r, r_sun_frac * S_R) - 1)])


total = electron_func(1)
print(total)

delta = 1/len(r)
electron_count = \
    np.array([[i, electron_func(i)] for i in np.arange(0, 1, delta)])

electron_count_half_rad = \
    electron_count[val_index(electron_count, 0.5 * total), 0]

print(electron_count_half_rad)

M_E = 5.97e+24
percentage_mass_earth = np.array([
    [32.1, 26, 55.845],  # iron
    [30.1, 8, 15.999],  # oxygen
    [15.1, 14, 28.0855],  # silicon
    [13.9, 12, 24.305],  # magnesium
    [2.9, 16, 32.065],  # sulphur
    [1.8, 28, 58.6934],  # nickel
    [1.5, 20, 40.078],  # calcium
    [1.4, 13, 26.981539]])  # aluminum

electron_count_earth = np.sum([
    (M_E * percentage_mass_earth[i, 0]/100) /
    (percentage_mass_earth[i, 2]/1000) *
    percentage_mass_earth[i, 1] * A_N
    for i in range(1, len(percentage_mass_earth) - 1)])

equivalent_sun_rad = electron_count[val_index(electron_count, electron_count_earth), 0]

print(equivalent_sun_rad)

my_x_ticks = (R_Rsun_ratio.max() - R_Rsun_ratio.min()) / 15
my_y_ticks = (norm_electron_density.max() - norm_electron_density.min()) / 20

plt.xticks(np.arange(R_Rsun_ratio.min(), R_Rsun_ratio.max(), my_x_ticks))
plt.yticks(np.arange(norm_electron_density.min(), norm_electron_density.max(), my_y_ticks))
plt.plot(R_Rsun_ratio, norm_electron_density)

plt.xlabel('log(n_e/N_A)')
plt.ylabel('R/R_sun')
plt.title('BS05(OP) Electron Density\n2005, ApJ, 621, L85 (astro-ph/0412440)')
plt.show()
