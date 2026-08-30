# Stochastic Predator-Prey Modelling

This is a numerical simulation of the Lotka-Volterra predator-prey model in Python with NumPy and Matplotlib. It works through solving a deterministic form through the 4th order Runge-Kutta method, and a stochastic form with the Euler-Maruyama method. Finally, the stochastic Lotka-Volterra model is simulated an arbitrary number of times during a Monte Carlo simulation to analyse qualitatively how variation in noise can effect extinction risk.

<img width="1904" height="1133" alt="image" src="https://github.com/user-attachments/assets/e15fa9e2-4bf7-4663-9f7b-b8cfd24b9f62" />

# Installation
You will need to install the Python modules NumPy and Matplotlib 


# Usage
'lv_deterministic.py' provides the deterministic model solved with 4th-order Runge-Kutta
'lv_stochastic.py' provides the stochastic model with multiplicative noise, solved with Euler-Maruyama
'lv_monte_carlo.py' includes a loop that iterates an arbitrary amount of times to ran a Monte Carlo simulation

All Python files output a png image of the timeseries/phase plot after the model has been simulated


# Overview
The Lotka-Volterra is among the oldest system of equations that exist to model the population dynamics between a predator and a prey. It is a pair of coupled, non-linear ordinary differential equations that must be solved numerically in lieu of a general analytical solution. This project aims to:
- Formulate the two differential equations in Python and solve them numerically
- Improve the biological realism of the model
- Reconstruct the model as a pair of stochastic differential equations to capture environmental randomness and solve them numerically
- Run Monte Carlo simulations across an arbitrary number of realisations

# Methodology
In the deterministic model, the base system, assumed to be:

$$ \frac{dx}{dt} = \alpha x - \beta xy \qquad \frac{dy}{dt} = \delta xy - \gamma y $$

is solved numerically with the 4th-order Runge-Kutta method built from scratch. This choice over Euler's method for solving differential equations is motivated by the Runge-Kutta's lower global error and its superior mathematical accuracy. 

The biological realism is improved with the incorporation of carrying capacity to prevent unbounded growth in the prey population, and a Holling Type-II functional response which caps predation rates at high prey density by accounting for predator handling time. On a phase graph, this is reflected by a phase profile that spirals inward, eventually converging towards a single equilibrium point at which the population of both species remain stable.

In the stochastic model, environmental randomness is introduced through the Wiener system, a mathematical representation of Brownian Motion. The pair of coupled stochastic differential equations (SDEs):

$$ dX = (\alpha x - \beta xy)dt + \sigma_x x dW_x $$
$$ dY = (\delta xy - \gamma y)dt + \sigma_y y dW_y $$

is made up of a deterministic drift term and a stochastic diffusion term. The latter achieves randomness in dW, a Wiener process increment. The value of dW is determined randomly from a normal distribution with mean 0 and variance equal to the time increment (dt). Stochasticity is achieved during every single iteration of the Euler-Maruyama method because a new value of dW is computed each time. This gives every single simulation of the SDEs a unique array of solutions which translates to unpredictability and mimics actual environmental randomness.

Since Brownian motion is continuous everywhere and differentiable nowhere, Euler-Maruyama method  with a random term correctly scaled to the square root of dt is used to solve this pair of SDEs

# Monte Carlo Analysis
Since every simulation of the stochastic model produces a unique sample from a distribution of outcomes, the model is run an arbitrary number of times via iterating through a for loop. Each iteration produces a unique plot, which is overlaid on top of every other plot on the same set of axes. The results reveal:
- Population trajectories of both species (predator and prey) desynchronises over time, described as a "phase drift"
- A widening uncertainty envelope as time progresses as a result of the phase drift
- Extinction might occur when the noise multipliers (sigma) are given a large enough value - an uncertainty the deterministic model cannot produce
However, it is important to observe that the trajectories of both populations are initially tightly synchronised, but then visibly drift out of phase over time causing the spread between runs to widen. With a deterministic model, the initial starting population already determines the future with absolute certainty

# Example results
Deterministic:
<img width="1879" height="1125" alt="image" src="https://github.com/user-attachments/assets/5003c96a-95f4-4d51-befe-3087c9129c16" />

Stochastic (Monte Carlo Simulation):
<img width="1870" height="1098" alt="image" src="https://github.com/user-attachments/assets/b7bc068d-b68d-42f3-bcc5-9f4e154f298f" />

# Reflection
A glaring limitation of this project is that the parameters that govern the growth and encounter rates of both populations are arbitrarily defined, not fitted to real ecological data (such as those between foxes and rabbits, or between lynx and hares). Furthermore, SDEs can be solved with higher accuracy through more sophisticated methods such the Milstein Scheme. 

It is also noteworthy that the noise terms of each population are independent of each other. This means that an environmental event affecting both populations (such as drought) are not accounted for. A natural extension of this project is to research into modelling the growth and encounter rate as periodic to take seasons or climate change into greater consideration - non-autonomous (time-dependent) parameters. A more ambitious perspective is to see if it is possible to extend the model to represent a multi-species systems to simulate food chain interactions.






