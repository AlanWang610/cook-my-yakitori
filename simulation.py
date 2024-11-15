import numpy as np
import scipy as sp
import matplotlib.pyplot as plt

def flux(temperature_C):
    # convert to kelvin
    temperature_K = temperature_C + 273.15
    emissivity = 0.77
    # stefan-boltzmann law
    M = 5.67e-8 * emissivity * temperature_K ** 4
    return M

R = 0.01
L=0.1
# build cylindrical mesh
# set the number of divisions in the r direction
n_r = 30
d_r = R / n_r
# set the number of time steps
n_t = 500
# set the thermal conductivity of chicken meat (Siripon et al. 2007)
alpha = 0.4930
# set the time step
d_t = d_r ** 2 / (2 * alpha)
print('dt: ', d_t)
# set the thermal diffusivity of chicken meat (Tavman et al. 2007)
k = alpha / (0.96 * 3.447)
# set the thermal
# calculate flux at the location of cooking
heat_source = 593
# divide by 2 to spread out the heat across the whole cylinder
flux_at_cooking_location = flux(heat_source) / 2
print('flux at cooking location in W/m^2: ', flux_at_cooking_location)

def cooking_simulation(T, n_t, n_r, alpha, d_t, d_r, r, flux_at_cooking_location):
    for t in range(1, n_t):
        T_temp = np.zeros(n_r)
        # update radial grid points using finite differences
        for i in range(1, n_r - 1):
            T_temp[i] = T[t-1, i] + alpha * d_t * ((T[t-1, i+1] - 2*T[t-1, i] + T[t-1, i-1]) / d_r ** 2)
        # update center assuming insulated boundary
        T_temp[0] = T_temp[1]
        # update outer boundary with heat flux
        T_temp[-1] = T_temp[-2] + flux_at_cooking_location * d_r / alpha
        T[t, :] = T_temp
    return T

r = np.linspace(d_r / 2, R - d_r / 2, n_r)
T_initial = 20 * np.ones(len(r))
T = np.zeros((n_t, n_r))
T[0, :] = T_initial

T = cooking_simulation(T, n_t, n_r, alpha, d_t, d_r, r, flux_at_cooking_location)
# plot T over time as a set of colored bars
plt.figure()
plt.imshow(T, aspect='auto')
plt.colorbar().set_label('Temperature (C)')
plt.xlabel('radius (dr)')
plt.ylabel('time (dt)')
plt.title('Yakitori Cooking Simulation')
plt.show()
print('Total time: ', n_t * d_t, 'seconds')