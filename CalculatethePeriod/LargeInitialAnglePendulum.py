import numpy as np
import matplotlib.pyplot as plt

# constants
g = 9.81    # gravitational field strength

# function of derivation (doing theta and omega at once)
def derivation(theta, omega, t, l):
    dtheta_dh = omega
    domega_dh = -g * np.sin(theta) / l
    return dtheta_dh, domega_dh

#function of RK4
def rk4(theta, omega, t, h, l):
    print(theta, omega)
    # calculation of k1 k2 k3 k4
    k1_theta, k1_omega = derivation(theta, omega, t, l)
    k2_theta, k2_omega = derivation(theta + 0.5 * k1_theta * h, omega + 0.5 * k1_omega * h, t + 0.5 * h, l)
    k3_theta, k3_omega = derivation(theta + 0.5 * k2_theta * h, omega + 0.5 * k2_omega * h, t + 0.5 * h, l)
    k4_theta, k4_omega = derivation(theta + k3_theta * h, omega + k3_omega * h, t + h, l)

    # calculation of result
    next_theta = theta + (h / 6) * (k1_theta + 2 * k2_theta + 2 * k3_theta + k4_theta)
    next_omega = omega + (h / 6) * (k1_omega + 2 * k2_omega + 2 * k3_omega + k4_omega)
    print(next_theta, next_omega)
    return next_theta, next_omega

def calculatePeriod(theta0, h, l, g):
    k = np.sin(theta0/2)
    sum = 0
    an = 1 * 2 * 3.1415926 * np.sqrt(l/g) #the first element in the series
    n = 0
    while an > 0.1 * h:
        #end when an is smaller than 0.1 * h since h is the minimum time span of our simulation
        sum += an
        n += 1
        an = an * ((2 * n - 1) / (2 * n) * k) ** 2 #calculate a(n+1) with an
    return sum

#control variable: initial theta (from 15 degree to 165 degree, step = 15 degree)
#independent variables: initial omega = 0, l = 1
initial_theta = np.radians(15)
initial_omega = 0
theta_step = np.radians(15)
l = 1
angle_period = []
while (initial_theta <= np.radians(166)):
    ts = [0]
    h = 0.0001
    thetas = [initial_theta]
    omegas = [initial_omega]
    
    #set ending condition: when omega changes between positive and negative for the second time
    count = 0
    while (count < 2):
        t_next = ts[-1] + h
        theta_next, omega_next = rk4(thetas[-1], omegas[-1], ts[-1], h, l)
        if (omegas[-1] * omega_next < 0):
            count += 1
        ts.append(t_next)
        thetas.append(theta_next)
        omegas.append(omega_next)
    angle_period.append(
        (round(np.rad2deg(initial_theta)),
         "%.4f"%calculatePeriod(initial_theta, h, l, g),
         "%.4f"%ts[-1]))
    initial_theta += theta_step

for theta_init, calc_period, simu_period in angle_period:
    print(f"initial theta = {theta_init} degree:\ncalculated period is {calc_period}, simulated period is {simu_period}")