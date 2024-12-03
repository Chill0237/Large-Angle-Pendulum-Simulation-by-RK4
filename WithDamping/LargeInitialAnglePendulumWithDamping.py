import numpy as np
import matplotlib.pyplot as plt

# constants
g = 9.81    # gravitational field strength

# function of derivation (doing theta and omega at once)
def derivation(theta, omega, t, l, c):
    dtheta_dh = omega
    domega_dh = -g * np.sin(theta) / l - c * omega
    return dtheta_dh, domega_dh

#function of RK4
def rk4(theta, omega, t, h, l, c):
    print(theta, omega)
    # calculation of k1 k2 k3 k4
    k1_theta, k1_omega = derivation(theta, omega, t, l, c)
    k2_theta, k2_omega = derivation(theta + 0.5 * k1_theta * h, omega + 0.5 * k1_omega * h, t + 0.5 * h, l, c)
    k3_theta, k3_omega = derivation(theta + 0.5 * k2_theta * h, omega + 0.5 * k2_omega * h, t + 0.5 * h, l, c)
    k4_theta, k4_omega = derivation(theta + k3_theta * h, omega + k3_omega * h, t + h, l, c)

    # calculation of result
    next_theta = theta + (h / 6) * (k1_theta + 2 * k2_theta + 2 * k3_theta + k4_theta)
    next_omega = omega + (h / 6) * (k1_omega + 2 * k2_omega + 2 * k3_omega + k4_omega)
    print(next_theta, next_omega)
    return next_theta, next_omega

#build figure and assignments
fig1 = plt.figure(figsize = (16, 16))
fig1.subplots_adjust(wspace=0.25, hspace=0.25)
omega_theta_c1 = fig1.add_subplot(4, 3, 1)
plt.xlabel("theta(rad)", fontsize = 10)
plt.ylabel("omega(rad)", fontsize = 10)
plt.title("omega-theta, with c=2", fontsize = 10)
theta_t_c1 = fig1.add_subplot(4, 3, 2)
plt.xlabel("t(s)", fontsize = 10)
plt.ylabel("theta(rad)", fontsize = 10)
plt.title("theta-t, with c=2", fontsize = 10)
omega_t_c1 = fig1.add_subplot(4, 3, 3)
plt.xlabel("t(s)", fontsize = 10)
plt.ylabel("omega(rad)", fontsize = 10)
plt.title("omega-t, with c=2", fontsize = 10)
omega_theta_c2 = fig1.add_subplot(4, 3, 4)
plt.xlabel("theta(rad)", fontsize = 10)
plt.ylabel("omega(rad)", fontsize = 10)
plt.title("omega-theta, with c=4", fontsize = 10)
theta_t_c2 = fig1.add_subplot(4, 3, 5)
plt.xlabel("t(s)", fontsize = 10)
plt.ylabel("theta(rad)", fontsize = 10)
plt.title("theta-t, with c=4", fontsize = 10)
omega_t_c2 = fig1.add_subplot(4, 3, 6)
plt.xlabel("t(s)", fontsize = 10)
plt.ylabel("omega(rad)", fontsize = 10)
plt.title("omega-t, with c=4", fontsize = 10)
omega_theta_c3 = fig1.add_subplot(4, 3, 7)
plt.xlabel("theta(rad)", fontsize = 10)
plt.ylabel("omega(rad)", fontsize = 10)
plt.title("omega-theta, with c=6", fontsize = 10)
theta_t_c3 = fig1.add_subplot(4, 3, 8)
plt.xlabel("t(s)", fontsize = 10)
plt.ylabel("theta(rad)", fontsize = 10)
plt.title("theta-t, with c=6", fontsize = 10)
omega_t_c3 = fig1.add_subplot(4, 3, 9)
plt.xlabel("t(s)", fontsize = 10)
plt.ylabel("omega(rad)", fontsize = 10)
plt.title("omega-t, with c=6", fontsize = 10)
omega_theta_c4 = fig1.add_subplot(4, 3, 10)
plt.xlabel("theta(rad)", fontsize = 10)
plt.ylabel("omega(rad)", fontsize = 10)
plt.title("omega-theta, with c=8", fontsize = 10)
theta_t_c4 = fig1.add_subplot(4, 3, 11)
plt.xlabel("t(s)", fontsize = 10)
plt.ylabel("theta(rad)", fontsize = 10)
plt.title("theta-t, with c=8", fontsize = 10)
omega_t_c4 = fig1.add_subplot(4, 3, 12)
plt.xlabel("t(s)", fontsize = 10)
plt.ylabel("omega(rad)", fontsize = 10)
plt.title("omega-t, with c=8", fontsize = 10)

#control variable: initial theta (from 15 degree to 165 degree, step = 15 degree)
#independent variables: initial omega = 0, l = 1, c = 1
initial_theta = np.radians(15)
initial_omega = 0
theta_step = np.radians(15)
l = 1
c = 2
while (initial_theta <= np.radians(166)):
    ts = [0]
    h = 0.0001
    thetas = [initial_theta]
    omegas = [initial_omega]
    
    #set ending condition: when t = 6
    while (ts[-1] <= 6):
        t_next = ts[-1] + h
        theta_next, omega_next = rk4(thetas[-1], omegas[-1], ts[-1], h, l, c)
        ts.append(t_next)
        thetas.append(theta_next)
        omegas.append(omega_next)
    omega_theta_c1.plot(thetas, omegas, linewidth = 1, linestyle = "None", marker = ".", markersize = 0.2)
    theta_t_c1.plot(ts, thetas, linewidth = 1, linestyle = "None", marker = ".", markersize = 0.2)
    omega_t_c1.plot(ts, omegas, linewidth = 1, linestyle = "None", marker = ".", markersize = 0.2)
    initial_theta += theta_step

#control variable: initial theta (from 15 degree to 165 degree, step = 15 degree)
#independent variables: initial omega = 0, l = 1, c = 2
initial_theta = np.radians(15)
initial_omega = 0
theta_step = np.radians(15)
l = 1
c = 4
while (initial_theta <= np.radians(166)):
    ts = [0]
    h = 0.0001
    thetas = [initial_theta]
    omegas = [initial_omega]
    
    #set ending condition: when t = 6
    while (ts[-1] <= 6):
        t_next = ts[-1] + h
        theta_next, omega_next = rk4(thetas[-1], omegas[-1], ts[-1], h, l, c)
        ts.append(t_next)
        thetas.append(theta_next)
        omegas.append(omega_next)
    omega_theta_c2.plot(thetas, omegas, linewidth = 1, linestyle = "None", marker = ".", markersize = 0.2)
    theta_t_c2.plot(ts, thetas, linewidth = 1, linestyle = "None", marker = ".", markersize = 0.2)
    omega_t_c2.plot(ts, omegas, linewidth = 1, linestyle = "None", marker = ".", markersize = 0.2)
    initial_theta += theta_step

#control variable: initial theta (from 15 degree to 165 degree, step = 15 degree)
#independent variables: initial omega = 0, l = 1, c = 3
initial_theta = np.radians(15)
initial_omega = 0
theta_step = np.radians(15)
l = 1
c = 6
while (initial_theta <= np.radians(166)):
    ts = [0]
    h = 0.0001
    thetas = [initial_theta]
    omegas = [initial_omega]
    
    #set ending condition: when t = 6
    while (ts[-1] <= 6):
        t_next = ts[-1] + h
        theta_next, omega_next = rk4(thetas[-1], omegas[-1], ts[-1], h, l, c)
        ts.append(t_next)
        thetas.append(theta_next)
        omegas.append(omega_next)
    omega_theta_c3.plot(thetas, omegas, linewidth = 1, linestyle = "None", marker = ".", markersize = 0.2)
    theta_t_c3.plot(ts, thetas, linewidth = 1, linestyle = "None", marker = ".", markersize = 0.2)
    omega_t_c3.plot(ts, omegas, linewidth = 1, linestyle = "None", marker = ".", markersize = 0.2)
    initial_theta += theta_step

#control variable: initial theta (from 15 degree to 165 degree, step = 15 degree)
#independent variables: initial omega = 0, l = 1, c = 4
initial_theta = np.radians(15)
initial_omega = 0
theta_step = np.radians(15)
l = 1
c = 8
while (initial_theta <= np.radians(166)):
    ts = [0]
    h = 0.0001
    thetas = [initial_theta]
    omegas = [initial_omega]
    
    #set ending condition: when t = 6
    while (ts[-1] <= 6):
        t_next = ts[-1] + h
        theta_next, omega_next = rk4(thetas[-1], omegas[-1], ts[-1], h, l, c)
        ts.append(t_next)
        thetas.append(theta_next)
        omegas.append(omega_next)
    omega_theta_c4.plot(thetas, omegas, linewidth = 1, linestyle = "None", marker = ".", markersize = 0.2)
    theta_t_c4.plot(ts, thetas, linewidth = 1, linestyle = "None", marker = ".", markersize = 0.2)
    omega_t_c4.plot(ts, omegas, linewidth = 1, linestyle = "None", marker = ".", markersize = 0.2)
    initial_theta += theta_step

fig1.savefig("With Damping/LargeInitialAnglePendulumWithDamping.png")