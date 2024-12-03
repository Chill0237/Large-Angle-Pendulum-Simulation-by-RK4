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

#build figure and assignments
fig1 = plt.figure(figsize = (16, 16))
fig1.subplots_adjust(wspace=0.25, hspace=0.25)
omega_theta_initial_theta = fig1.add_subplot(4, 3, 1)
plt.xlabel("theta(rad)", fontsize = 10)
plt.ylabel("omega(rad)", fontsize = 10)
plt.title("omega-theta, when different initial thetas", fontsize = 10)
theta_t_initial_theta = fig1.add_subplot(4, 3, 2)
plt.xlabel("t(s)", fontsize = 10)
plt.ylabel("theta(rad)", fontsize = 10)
plt.title("theta-t, when different initial thetas", fontsize = 10)
omega_t_initial_theta = fig1.add_subplot(4, 3, 3)
plt.xlabel("t(s)", fontsize = 10)
plt.ylabel("omega(rad)", fontsize = 10)
plt.title("omega-t, when different initial thetas", fontsize = 10)
omega_theta_initial_omega = fig1.add_subplot(4, 3, 4)
plt.xlabel("theta(rad)", fontsize = 10)
plt.ylabel("omega(rad)", fontsize = 10)
plt.title("omega-theta, when different initial omegas", fontsize = 10)
theta_t_initial_omega = fig1.add_subplot(4, 3, 5)
plt.xlabel("t(s)", fontsize = 10)
plt.ylabel("theta(rad)", fontsize = 10)
plt.title("theta-t, when different initial omegas", fontsize = 10)
omega_t_initial_omega = fig1.add_subplot(4, 3, 6)
plt.xlabel("t(s)", fontsize = 10)
plt.ylabel("omega(rad)", fontsize = 10)
plt.title("omega-t, when different initial omegas", fontsize = 10)
omega_theta_length_theta = fig1.add_subplot(4, 3, 7)
plt.xlabel("theta(rad)", fontsize = 10)
plt.ylabel("omega(rad)", fontsize = 10)
plt.title("omega-theta, when different lengths and set theta", fontsize = 10)
theta_t_length_theta = fig1.add_subplot(4, 3, 8)
plt.xlabel("t(s)", fontsize = 10)
plt.ylabel("theta(rad)", fontsize = 10)
plt.title("theta-t, when different lengths and set theta", fontsize = 10)
omega_t_length_theta = fig1.add_subplot(4, 3, 9)
plt.xlabel("t(s)", fontsize = 10)
plt.ylabel("omega(rad)", fontsize = 10)
plt.title("omega-t, when different lengths and set theta", fontsize = 10)
omega_theta_length_omega = fig1.add_subplot(4, 3, 10)
plt.xlabel("theta(rad)", fontsize = 10)
plt.ylabel("omega(rad)", fontsize = 10)
plt.title("omega-theta, when different lengths and set omega", fontsize = 10)
theta_t_length_omega = fig1.add_subplot(4, 3, 11)
plt.xlabel("t(s)", fontsize = 10)
plt.ylabel("theta(rad)", fontsize = 10)
plt.title("theta-t, when different lengths and set omega", fontsize = 10)
omega_t_length_omega = fig1.add_subplot(4, 3, 12)
plt.xlabel("t(s)", fontsize = 10)
plt.ylabel("omega(rad)", fontsize = 10)
plt.title("omega-t, when different lengths and set omega", fontsize = 10)

#build figure of percent of period and percent of amplitude
fig2 = plt.figure(figsize = (16, 16))

fig2.subplots_adjust(wspace=0.25, hspace=0.25)
theta_t_initial_theta_2 = fig2.add_subplot(4, 2, 1)
plt.xlabel("t(percent of period)", fontsize = 10)
plt.ylabel("theta(percent of amplitude)", fontsize = 10)
plt.title("theta-t, when different initial thetas", fontsize = 10)
omega_t_initial_theta_2 = fig2.add_subplot(4, 2, 2)
plt.xlabel("t(percent of period)", fontsize = 10)
plt.ylabel("omega(percent of amplitude)", fontsize = 10)
plt.title("omega-t, when different initial thetas", fontsize = 10)
theta_t_initial_omega_2 = fig2.add_subplot(4, 2, 3)
plt.xlabel("t(percent of period)", fontsize = 10)
plt.ylabel("theta(percent of amplitude)", fontsize = 10)
plt.title("theta-t, when different initial omegas", fontsize = 10)
omega_t_initial_omega_2 = fig2.add_subplot(4, 2, 4)
plt.xlabel("t(percent of period)", fontsize = 10)
plt.ylabel("omega(percent of amplitude)", fontsize = 10)
plt.title("omega-t, when different initial omegas", fontsize = 10)
theta_t_length_theta_2 = fig2.add_subplot(4, 2, 5)
plt.xlabel("t(percent of period)", fontsize = 10)
plt.ylabel("theta(percent of amplitude)", fontsize = 10)
plt.title("theta-t, when different lengths and set theta", fontsize = 10)
omega_t_length_theta_2 = fig2.add_subplot(4, 2, 6)
plt.xlabel("t(percent of period)", fontsize = 10)
plt.ylabel("omega(percent of amplitude)", fontsize = 10)
plt.title("omega-t, when different lengths and set theta", fontsize = 10)
theta_t_length_omega_2 = fig2.add_subplot(4, 2, 7)
plt.xlabel("t(percent of period)", fontsize = 10)
plt.ylabel("theta(percent of amplitude)", fontsize = 10)
plt.title("theta-t, when different lengths and set omega", fontsize = 10)
omega_t_length_omega_2 = fig2.add_subplot(4, 2, 8)
plt.xlabel("t(percent of period)", fontsize = 10)
plt.ylabel("omega(percent of amplitude)", fontsize = 10)
plt.title("omega-t, when different lengths and set omega", fontsize = 10)

#control variable: initial theta (from 15 degree to 165 degree, step = 15 degree)
#independent variables: initial omega = 0, l = 1
initial_theta = np.radians(15)
initial_omega = 0
theta_step = np.radians(15)
l = 1
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
    omega_theta_initial_theta.plot(thetas, omegas, linewidth = 1, linestyle = "None", marker = ".", markersize = 0.2)
    theta_t_initial_theta.plot(ts, thetas, linewidth = 1, linestyle = "None", marker = ".", markersize = 0.2)
    omega_t_initial_theta.plot(ts, omegas, linewidth = 1, linestyle = "None", marker = ".", markersize = 0.2)

    #resize the graph for better comparison
    ts_2 = [i / ts[-1] * 100 for i in ts]
    max_thetas = max(thetas)
    max_omegas = max(omegas)
    thetas_2 = [i / max_thetas for i in thetas]
    omegas_2 = [i / max_omegas for i in omegas]
    theta_t_initial_theta_2.plot(ts_2, thetas_2, linewidth = 1, linestyle = "None", marker = ".", markersize = 0.2)
    omega_t_initial_theta_2.plot(ts_2, omegas_2, linewidth = 1, linestyle = "None", marker = ".", markersize = 0.2)

    initial_theta += theta_step

#control variable: initial omega (from 30 degree/s to 390 degree/s, step = 30 degree/s)
#independent variables: initial theta = 0, l = 1
initial_theta = 0
initial_omega = np.radians(30)
omega_step = np.radians(30)
l = 1
while (initial_omega <= np.radians(391)):
    ts = [0]
    h = 0.0001
    thetas = [initial_theta]
    omegas = [initial_omega]
    
    #set ending condition: when theta changes between positive and negative for the second time
    count = 0
    while (count < 2):
        t_next = ts[-1] + h
        theta_next, omega_next = rk4(thetas[-1], omegas[-1], ts[-1], h, l)
        if (thetas[-1] * theta_next < 0):
            count += 1
        if (theta_next < np.radians(-180)):
            theta_next += np.radians(360)
        if (theta_next > np.radians(180)):
            theta_next -= np.radians(360)
        ts.append(t_next)
        thetas.append(theta_next)
        omegas.append(omega_next)
    omega_theta_initial_omega.plot(thetas, omegas, linewidth = 1, linestyle = "None", marker = ".", markersize = 0.2)
    theta_t_initial_omega.plot(ts, thetas, linewidth = 1, linestyle = "None", marker = ".", markersize = 0.2)
    omega_t_initial_omega.plot(ts, omegas, linewidth = 1, linestyle = "None", marker = ".", markersize = 0.2)

    #resize the graph for better comparison, ignore the rotating cases
    if (initial_omega <= np.radians(331)):
        ts_2 = [i / ts[-1] * 100 for i in ts]
        max_thetas = max(thetas)
        max_omegas = max(omegas)
        thetas_2 = [i / max_thetas for i in thetas]
        omegas_2 = [i / max_omegas for i in omegas]
        theta_t_initial_omega_2.plot(ts_2, thetas_2, linewidth = 1, linestyle = "None", marker = ".", markersize = 0.2)
        omega_t_initial_omega_2.plot(ts_2, omegas_2, linewidth = 1, linestyle = "None", marker = ".", markersize = 0.2)

    initial_omega += omega_step

#control variable: l (from 0.2 meter to 2 meter, step = 0.2 meter)
#independent variables: initial theta = 165 degree, initial omega = 0 degree/s
initial_theta = np.radians(165)
initial_omega = 0
l_step = 0.2
l = 0.2
while (l <= 2):
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
    omega_theta_length_theta.plot(thetas, omegas, linewidth = 1, linestyle = "None", marker = ".", markersize = 0.2)
    theta_t_length_theta.plot(ts, thetas, linewidth = 1, linestyle = "None", marker = ".", markersize = 0.2)
    omega_t_length_theta.plot(ts, omegas, linewidth = 1, linestyle = "None", marker = ".", markersize = 0.2)

    #resize the graph for better comparison
    ts_2 = [i / ts[-1] * 100 for i in ts]
    max_thetas = max(thetas)
    max_omegas = max(omegas)
    thetas_2 = [i / max_thetas for i in thetas]
    omegas_2 = [i / max_omegas for i in omegas]
    theta_t_length_theta_2.plot(ts_2, thetas_2, linewidth = 1, linestyle = "None", marker = ".", markersize = 0.2)
    omega_t_length_theta_2.plot(ts_2, omegas_2, linewidth = 1, linestyle = "None", marker = ".", markersize = 0.2)

    l += l_step

#control variable: l (from 0.2 meter to 2 meter, step = 0.2 meter)
#independent variables: initial theta = 0 degree, initial omega = 245 degree/s
initial_theta = 0
initial_omega = np.radians(245)
l_step = 0.2
l = 0.2
while (l <= 2):
    ts = [0]
    h = 0.0001
    thetas = [initial_theta]
    omegas = [initial_omega]
    
    #set ending condition: when omega changes between positive and negative for the second time
    count = 0
    while (count < 2):
        t_next = ts[-1] + h
        theta_next, omega_next = rk4(thetas[-1], omegas[-1], ts[-1], h, l)
        if (thetas[-1] * theta_next < 0):
            count += 1
        ts.append(t_next)
        thetas.append(theta_next)
        omegas.append(omega_next)
    omega_theta_length_omega.plot(thetas, omegas, linewidth = 1, linestyle = "None", marker = ".", markersize = 0.2)
    theta_t_length_omega.plot(ts, thetas, linewidth = 1, linestyle = "None", marker = ".", markersize = 0.2)
    omega_t_length_omega.plot(ts, omegas, linewidth = 1, linestyle = "None", marker = ".", markersize = 0.2)

    #resize the graph for better comparison
    ts_2 = [i / ts[-1] * 100 for i in ts]
    max_thetas = max(thetas)
    max_omegas = max(omegas)
    thetas_2 = [i / max_thetas for i in thetas]
    omegas_2 = [i / max_omegas for i in omegas]
    theta_t_length_omega_2.plot(ts_2, thetas_2, linewidth = 1, linestyle = "None", marker = ".", markersize = 0.2)
    omega_t_length_omega_2.plot(ts_2, omegas_2, linewidth = 1, linestyle = "None", marker = ".", markersize = 0.2)

    l += l_step

fig1.savefig("Large Initial Angle Pendulum/LargeInitialAnglePendulum.png")
fig2.savefig("Large Initial Angle Pendulum/LargeInitialAnglePendulum2.png")