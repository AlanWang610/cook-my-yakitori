import numpy as np
import scipy as sp
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.animation as animation
plt.rcParams['animation.ffmpeg_path'] = r'C:\\ffmpeg\\bin\\ffmpeg.exe'

def flux(temperature_C = 593, outer_temperature_C = 20):
    # convert to kelvin
    temperature_K = temperature_C + 273.15
    outer_temperature_K = outer_temperature_C + 273.15
    emissivity = 0.77
    # stefan-boltzmann law
    M = 5.67e-8 * emissivity * (temperature_K ** 4 - outer_temperature_K ** 4)
    return M

R = 0.01
L= 0.1
# build cylindrical mesh
# set the number of divisions in the r direction
n_r = 30
d_r = R / n_r
# set the number of time steps
n_t = 5000
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
flux_at_cooking_location = flux(heat_source, 20) / 2
print('flux at cooking location in W/m^2: ', flux_at_cooking_location)
# divide by 2 to spread out the heat across the whole cylinder

def cooking_simulation(T, n_t, n_r, alpha, d_t, d_r, r, heat_source, flux):
    for t in range(1, n_t):
        T_temp = np.zeros(n_r)
        # update radial grid points using finite differences
        for i in range(1, n_r - 1):
            T_temp[i] = T[t-1, i] + alpha * d_t * ((T[t-1, i+1] - 2*T[t-1, i] + T[t-1, i-1]) / d_r ** 2)
        # update center assuming insulated boundary
        T_temp[0] = T_temp[1]
        # update outer boundary with heat flux
        T_temp[-1] = T_temp[-2] + (flux(heat_source, T_temp[-1]) / 2) * d_r / alpha
        T[t, :] = T_temp
    return T

r = np.linspace(d_r / 2, R - d_r / 2, n_r)
T_initial = 20 * np.ones(len(r))
T = np.zeros((n_t, n_r))
T[0, :] = T_initial

T = cooking_simulation(T, n_t, n_r, alpha, d_t, d_r, r, heat_source, flux)
# # plot T over time as a set of colored bars
# plt.figure()
# plt.imshow(T, aspect='auto')
# plt.colorbar().set_label('Temperature (C)')
# plt.xlabel('radius (dr)')
# plt.ylabel('time (dt)')
# plt.title('Yakitori Cooking Simulation')
# plt.show()
# print('Total time: ', n_t * d_t, 'seconds')

def plot_cross_section(T, r, n_t, ax):
    T_cross_section = T[n_t, :]
    theta = np.linspace(0, 2 * np.pi, 100)
    r_grid, theta_grid = np.meshgrid(r, theta)
    X = r_grid * np.cos(theta_grid)
    Y = r_grid * np.sin(theta_grid)
    Z = np.tile(T_cross_section, (100, 1))

    ax.clear()  # Clear previous frame
    mesh = ax.pcolormesh(X, Y, Z, shading='auto')
    if not hasattr(plot_cross_section, "colorbar"):  # Add colorbar only once
        plot_cross_section.colorbar = plt.colorbar(mesh, ax=ax, label='Temperature (C)')
    ax.set_xlabel('x (m)')
    ax.set_ylabel('y (m)')
    ax.set_title('Temperature Distribution in Yakitori Cross Section')
    ax.axis('equal')

fig, ax = plt.subplots()

def update_plot(i):
    print(i)
    plot_cross_section(T, r, i, ax)

ani = animation.FuncAnimation(fig, update_plot, frames=range(0, n_t, n_t // 20))
ani.save('yakitori_cooking_simulation.mp4', writer='ffmpeg')
plt.close(fig)