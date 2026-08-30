# Importing modules
import numpy as np
import matplotlib.pyplot as plt

# ---------------------------------------------------------------------------
# Setting up initial conditions
alpha = 0.1                 # prey growth rate
beta = 0.0004               # predation rate on prey
gamma = 0.15                # predator death rate
delta = 0.00035             # predator growth per prey eaten
K = 5000                    # carrying capacity
h = 0.45                    # handling time
x0 = 1000                   # initial prey population
y0 = 500                    # initial predator population
#noise variables
sigma_x = 0.05              # prey noise
sigma_y = 0.035             # predator noise

# ---------------------------------------------------------------------------
# Setting up time domain
t_start, t_end = 0, 250                 # start time, end time
dt = 0.01                               # resolution of time
n_steps = int((t_end - t_start)/dt)     # number of steps
# time domain
t = np.linspace(t_start, t_end, n_steps + 1)

# ---------------------------------------------------------------------------
# Pre-allocating arrays to store solutions
x = np.zeros(n_steps + 1)
y = np.zeros(n_steps + 1)
# setting up the first item in the array as initial conditions
x[0] = x0
y[0] = y0

# ---------------------------------------------------------------------------
# Lotka-Volterra Model
def lotka_volterra(x,y):

    # The prey equation
    dxdt = alpha*x*(1 - x/K) - (beta*x*y)/(1 + beta*h*x)
    # The predator equation
    dydt = (delta*x*y)/(1 + beta*h*x) - gamma*y
    return dxdt, dydt

# ---------------------------------------------------------------------------
# Solving: Euler-Maruyama Method
def euler_maruyama():
    for i in range(n_steps):
        # Setting up the deterministic drift terms
        dxdt, dydt = lotka_volterra(x[i],y[i])

        # Mathematically modelling Brownian Motion
        dW_x = np.random.normal(0, np.sqrt(dt))
        dW_y = np.random.normal(0, np.sqrt(dt))

        # Setting up the stochastic diffusion terms
        noise_x = sigma_x * x[i] * dW_x
        noise_y = sigma_y * y[i] * dW_y

        # Iterate through:
        x[i+1] = x[i] + dxdt * dt + noise_x
        y[i+1] = y[i] + dydt * dt + noise_y

        # Guard against negative populations
        x[i+1] = max(x[i + 1], 0)
        y[i+1] = max(y[i + 1], 0)

# ---------------------------------------------------------------------------
# Monte Carlo Simulation with an arbitrary number of runs
simulations = 200
plt.figure(figsize=(12,8))

for i in range(simulations):
    # run the model
    euler_maruyama()
    # plot prey and predator curves and give them different colours
    plt.plot(t, x, color='steelblue', alpha=0.4, linewidth=0.8)
    plt.plot(t, y, color='darkorange', alpha=0.4, linewidth=0.8)

    # alternatively - code for phase plot:
    # plt.plot(x, y, color='steelblue', alpha=0.5, linewidth=0.8)

# ---------------------------------------------------------------------------
# Generating legend (real data already plotted above)
plt.xlabel('Time')
plt.ylabel('Population Size')
plt.title(f"Lotka-Volterra Monte Carlo ({simulations} simulations)")
plt.tight_layout()
plt.savefig("LV_Monte_Carlo.png", dpi=150)
plt.show()
