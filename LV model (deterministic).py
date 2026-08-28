# Importing modules
import numpy as np
import matplotlib.pyplot as plt

# Setting up initial conditions
alpha = 0.1                 # prey growth rate
beta = 0.0004               # predation rate on prey
gamma = 0.15                # predator death rate
delta = 0.00035             # predator growth per prey eaten
K = 5000                    # carrying capacity
h = 0.45                    # handling time
x0 = 1000                   # initial prey population
y0 = 500                    # initial predator population

# Setting up time domain--
t_start, t_end = 0, 250                 # start time, end time
dt = 0.01                               # resolution of time
n_steps = int((t_end - t_start)/dt)     # number of steps
# time domain
t = np.linspace(t_start, t_end, n_steps + 1)

# Pre-allocating arrays to store solution
x = np.zeros(n_steps + 1)
y = np.zeros(n_steps + 1)
# setting up the first item in the array as initial conditions
x[0] = x0
y[0] = y0

# Lotka-Volterra Model
def lotka_volterra(x,y):

    # The prey equation
    dxdt = alpha*x*(1 - x/K) - (beta*x*y)/(1 + beta*h*x)
        # carrying capacity and functional response

    # The predator equation
    dydt = (delta*x*y)/(1 + beta*h*x) - gamma*y
    return dxdt, dydt

# Integrating with Runge-Kutta 4th Order Method
# using a for loop to loop through all the time steps
for i in range(n_steps):
    # calculating k1, k2, k3, k4
    k1x, k1y = lotka_volterra(x[i],y[i])
    k2x, k2y = lotka_volterra(x[i] + (k1x*dt)/2, (y[i] + (k1y*dt)/2))
    k3x, k3y = lotka_volterra(x[i] + (k2x*dt)/2, (y[i] + (k2y*dt)/2))
    k4x, k4y = lotka_volterra(x[i] + (k3x*dt), (y[i] + (k3y*dt)))
    # adding results for x and y into the arrays
    x[i+1] = x[i] + (dt/6)*(k1x + 2*k2x + 2*k3x + k4x)
    y[i+1] = y[i] + (dt/6)*(k1y + 2*k2y + 2*k3y + k4y)

# Plotting time series
plt.figure(figsize=(12,8))
plt.plot(t, x, label="Prey Population")
plt.plot(t, y, label="Predator Population")
plt.xlabel('Time')
plt.ylabel('Population Size')
plt.title("Lotka-Volterra Predator-Prey Model")
plt.tight_layout()
plt.savefig("LV_deterministic_timeseries.png", dpi=150)
plt.show()

# Plotting phase portrait
plt.figure(figsize=(10, 10))
plt.plot(x, y)
plt.xlabel("Prey population")
plt.ylabel("Predator population")
plt.title("Phase Portrait")
plt.tight_layout()
plt.savefig("LV_deterministic_phase.png", dpi=150)
plt.show()
