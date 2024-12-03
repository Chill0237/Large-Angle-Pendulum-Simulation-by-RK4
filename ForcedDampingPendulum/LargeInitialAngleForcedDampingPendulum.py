import numpy as np
import matplotlib.pyplot as plt

# constants
g = 9.81    # gravitational field strength

# function of derivation (doing theta and omega at once)
def derivation(theta, omega, t, l, c, alpha_max, omega_alpha):
    dtheta_dh = omega
    domega_dh = -g * np.sin(theta) / l - c * omega + alpha_max * np.sin(omega_alpha * t)
    return dtheta_dh, domega_dh

#function of RK4
def rk4(theta, omega, t, h, l, c, alpha_max, omega_alpha):
    print(theta, omega)
    # calculation of k1 k2 k3 k4
    k1_theta, k1_omega = derivation(theta, omega, t, l, c, alpha_max, omega_alpha)
    k2_theta, k2_omega = derivation(theta + 0.5 * k1_theta * h, omega + 0.5 * k1_omega * h, t + 0.5 * h, l, c, alpha_max, omega_alpha)
    k3_theta, k3_omega = derivation(theta + 0.5 * k2_theta * h, omega + 0.5 * k2_omega * h, t + 0.5 * h, l, c, alpha_max, omega_alpha)
    k4_theta, k4_omega = derivation(theta + k3_theta * h, omega + k3_omega * h, t + h, l, c, alpha_max, omega_alpha)

    # calculation of result
    next_theta = theta + (h / 6) * (k1_theta + 2 * k2_theta + 2 * k3_theta + k4_theta)
    next_omega = omega + (h / 6) * (k1_omega + 2 * k2_omega + 2 * k3_omega + k4_omega)
    print(next_theta, next_omega)
    return next_theta, next_omega

#build figure and assignments
fig = plt.figure(figsize = (16, 16))
fig.subplots_adjust(wspace=0.25, hspace=0.25)
omega_theta_oa1 = fig.add_subplot(4, 3, 1)
plt.xlabel("theta(rad)", fontsize = 10)
plt.ylabel("omega(rad)", fontsize = 10)
plt.title("omega-theta, with omega_alpha = 0.5pi", fontsize = 10)
theta_t_oa1 = fig.add_subplot(4, 3, 2)
plt.xlabel("t(s)", fontsize = 10)
plt.ylabel("theta(rad)", fontsize = 10)
plt.title("theta-t, with omega_alpha = 0.5pi", fontsize = 10)
omega_t_oa1 = fig.add_subplot(4, 3, 3)
plt.xlabel("t(s)", fontsize = 10)
plt.ylabel("omega(rad)", fontsize = 10)
plt.title("omega-t, with omega_alpha = 0.5pi", fontsize = 10)
omega_theta_oa2 = fig.add_subplot(4, 3, 4)
plt.xlabel("theta(rad)", fontsize = 10)
plt.ylabel("omega(rad)", fontsize = 10)
plt.title("omega-theta, with omega_alpha = pi", fontsize = 10)
theta_t_oa2 = fig.add_subplot(4, 3, 5)
plt.xlabel("t(s)", fontsize = 10)
plt.ylabel("theta(rad)", fontsize = 10)
plt.title("theta-t, with omega_alpha = pi", fontsize = 10)
omega_t_oa2 = fig.add_subplot(4, 3, 6)
plt.xlabel("t(s)", fontsize = 10)
plt.ylabel("omega(rad)", fontsize = 10)
plt.title("omega-t, with omega_alpha = pi", fontsize = 10)
omega_theta_oa3 = fig.add_subplot(4, 3, 7)
plt.xlabel("theta(rad)", fontsize = 10)
plt.ylabel("omega(rad)", fontsize = 10)
plt.title("omega-theta, with omega_alpha = 1.5pi", fontsize = 10)
theta_t_oa3 = fig.add_subplot(4, 3, 8)
plt.xlabel("t(s)", fontsize = 10)
plt.ylabel("theta(rad)", fontsize = 10)
plt.title("theta-t, with omega_alpha = 1.5pi", fontsize = 10)
omega_t_oa3 = fig.add_subplot(4, 3, 9)
plt.xlabel("t(s)", fontsize = 10)
plt.ylabel("omega(rad)", fontsize = 10)
plt.title("omega-t, with omega_alpha = 1.5pi", fontsize = 10)
omega_theta_oa4 = fig.add_subplot(4, 3, 10)
plt.xlabel("theta(rad)", fontsize = 10)
plt.ylabel("omega(rad)", fontsize = 10)
plt.title("omega-theta, with omega_alpha = 2pi", fontsize = 10)
theta_t_oa4 = fig.add_subplot(4, 3, 11)
plt.xlabel("t(s)", fontsize = 10)
plt.ylabel("theta(rad)", fontsize = 10)
plt.title("theta-t, with omega_alpha = 2pi", fontsize = 10)
omega_t_oa4 = fig.add_subplot(4, 3, 12)
plt.xlabel("t(s)", fontsize = 10)
plt.ylabel("omega(rad)", fontsize = 10)
plt.title("omega-t, with omega_alpha = 2pi", fontsize = 10)

#control variable: c (from 2 to 8, step = 2)
#independent variables: initial theta = 0, initial omega = 0, l = 1, alpha_max = 1, omega_alpha = 0.5pi
initial_theta = 0
initial_omega = 0
l = 1
c = 2
c_step = 2
alpha_max = 1
omega_alpha = 0.5 * 3.1415926
while (c <= 8):
    ts = [0]
    h = 0.0001
    thetas = [initial_theta]
    omegas = [initial_omega]
    
    #set ending condition: when t = 12
    while (ts[-1] <= 12):
        t_next = ts[-1] + h
        theta_next, omega_next = rk4(thetas[-1], omegas[-1], ts[-1], h, l, c, alpha_max, omega_alpha)
        ts.append(t_next)
        thetas.append(theta_next)
        omegas.append(omega_next)
    omega_theta_oa1.plot(thetas, omegas, linewidth = 1, linestyle = "None", marker = ".", markersize = 0.2)
    theta_t_oa1.plot(ts, thetas, linewidth = 1, linestyle = "None", marker = ".", markersize = 0.2)
    omega_t_oa1.plot(ts, omegas, linewidth = 1, linestyle = "None", marker = ".", markersize = 0.2)
    c += c_step

#control variable: c (from 2 to 8, step = 2)
#independent variables: initial theta = 0, initial omega = 0, l = 1, alpha_max = 1, omega_alpha = pi
initial_theta = 0
initial_omega = 0
l = 1
c = 2
c_step = 2
alpha_max = 1
omega_alpha = 3.1415926
while (c <= 8):
    ts = [0]
    h = 0.0001
    thetas = [initial_theta]
    omegas = [initial_omega]
    
    #set ending condition: when t = 12
    while (ts[-1] <= 12):
        t_next = ts[-1] + h
        theta_next, omega_next = rk4(thetas[-1], omegas[-1], ts[-1], h, l, c, alpha_max, omega_alpha)
        ts.append(t_next)
        thetas.append(theta_next)
        omegas.append(omega_next)
    omega_theta_oa2.plot(thetas, omegas, linewidth = 1, linestyle = "None", marker = ".", markersize = 0.2)
    theta_t_oa2.plot(ts, thetas, linewidth = 1, linestyle = "None", marker = ".", markersize = 0.2)
    omega_t_oa2.plot(ts, omegas, linewidth = 1, linestyle = "None", marker = ".", markersize = 0.2)
    c += c_step

#control variable: c (from 2 to 8, step = 2)
#independent variables: initial theta = 0, initial omega = 0, l = 1, alpha_max = 1, omega_alpha = 1.5 * pi
initial_theta = 0
initial_omega = 0
l = 1
c = 2
c_step = 2
alpha_max = 1
omega_alpha = 1.5 * 3.1415926
while (c <= 8):
    ts = [0]
    h = 0.0001
    thetas = [initial_theta]
    omegas = [initial_omega]
    
    #set ending condition: when t = 12
    while (ts[-1] <= 12):
        t_next = ts[-1] + h
        theta_next, omega_next = rk4(thetas[-1], omegas[-1], ts[-1], h, l, c, alpha_max, omega_alpha)
        ts.append(t_next)
        thetas.append(theta_next)
        omegas.append(omega_next)
    omega_theta_oa3.plot(thetas, omegas, linewidth = 1, linestyle = "None", marker = ".", markersize = 0.2)
    theta_t_oa3.plot(ts, thetas, linewidth = 1, linestyle = "None", marker = ".", markersize = 0.2)
    omega_t_oa3.plot(ts, omegas, linewidth = 1, linestyle = "None", marker = ".", markersize = 0.2)
    c += c_step

#control variable: c (from 2 to 8, step = 2)
#independent variables: initial theta = 0, initial omega = 0, l = 1, alpha_max = 1, omega_alpha = 1.5 * pi
initial_theta = 0
initial_omega = 0
l = 1
c = 2
c_step = 2
alpha_max = 1
omega_alpha = 2 * 3.1415926
while (c <= 8):
    ts = [0]
    h = 0.0001
    thetas = [initial_theta]
    omegas = [initial_omega]
    
    #set ending condition: when t = 12
    while (ts[-1] <= 12):
        t_next = ts[-1] + h
        theta_next, omega_next = rk4(thetas[-1], omegas[-1], ts[-1], h, l, c, alpha_max, omega_alpha)
        ts.append(t_next)
        thetas.append(theta_next)
        omegas.append(omega_next)
    omega_theta_oa4.plot(thetas, omegas, linewidth = 1, linestyle = "None", marker = ".", markersize = 0.2)
    theta_t_oa4.plot(ts, thetas, linewidth = 1, linestyle = "None", marker = ".", markersize = 0.2)
    omega_t_oa4.plot(ts, omegas, linewidth = 1, linestyle = "None", marker = ".", markersize = 0.2)
    c += c_step

fig.savefig("Forced Damping Pendulum/LargeInitialAngleForcedDampingPendulum.png")