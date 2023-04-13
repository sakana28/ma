import numpy as np
from scipy.interpolate import interp1d

def bearing_signal_model_dist(d, D, contact_angle, n, fault_type, fc, fd, fm, fr, N, fs, SNR_dB, q_fault, q_stiffness, q_rotation):
    """
    Generation of a simulated signal for distributed fault in rolling element bearing.

    Inputs:
        d (float): bearing roller diameter [mm]
        D (float): pitch circle diameter [mm]
        contact_angle (float): contact angle [rad]
        n (int): number of rolling elements
        fault_type (str): fault type selection: 'inner', 'outer', 'ball'
        fc (float): row vector containing the carrier component of the speed
        fd (float): row vector containing the frequency deviation
        fm (float): row vector containing the modulation frequency
        fr (list): rotation frequency profile
        N (int): number of points per revolution
        fs (float): sampling frequency [Hz]
        SNR_dB (float): signal to noise ratio [dB]
        q_fault (float): amplitude modulation at the fault frequency
        q_stiffness (float): amplitude value of the deterministic component related to the stiffness variation
        q_rotation (float): amplitude value of the deterministic component related to the bearing rotation

    Returns:
        tuple: (t, x, x_noise, fr_time)
            t (numpy.ndarray): time signal [s]
            x (numpy.ndarray): simulated bearing signal without noise
            x_noise (numpy.ndarray): simulated bearing signal with noise
            fr_time (numpy.ndarray): speed profile in the time domain [Hz]
    """

    if fault_type == 'inner':
        geometry_parameter = 1 / 2 * (1 + d / D * np.cos(contact_angle))
    elif fault_type == 'outer':
        geometry_parameter = 1 / 2 * (1 - d / D * np.cos(contact_angle))
    elif fault_type == 'ball':
        geometry_parameter = 1 / (2 * n) * (1 - (d / D * np.cos(contact_angle)) ** 2) / (d / D)
    else:
        raise ValueError("Invalid fault type")

    L_theta = len(fr)
    theta = np.arange(0, L_theta) * 2 * np.pi / N
    theta_time = np.zeros(len(fr))

    for index in range(1, len(fr)):
        theta_time[index] = theta_time[index - 1] + (2 * np.pi / N) / (2 * np.pi * fr[index])

    L = int(np.floor(theta_time[-1] * fs))
    t = np.arange(0, L) / fs

    fr_time = interp1d(theta_time, fr, kind='spline')(t)

    x_rotation = q_rotation * np.cos(fc / fc * theta + fd / fc * (np.cumsum(np.cos(fm / fc * theta) / N)))
    x_rotation_time = interp1d(theta_time, x_rotation, kind='spline')(t)

    tau_stiffness = n / 2 * (1 - d / D * np.cos(contact_angle))
    x_stiffness = q_stiffness * np.cos(fc / fc * tau_stiffness * theta + fd / fc * tau_stiffness * (np.cumsum(np.cos(fm / fc * tau_stiffness * theta) / N)))
    x_stiffness_time = interp1d(theta_time, x_stiffness, kind='spline')(t)

    tau_fault = n * geometry_parameter
    q = 1 + q_fault * np.sin(fc / fc * tau_fault * theta + fd / fc * geometry_parameter * (np.cumsum(np.cos(fm / fc * geometry_parameter * theta) / N)))
    q_time = interp1d(theta_time, q, kind='spline')(t)
    x_fault_time = np.random.randn(1, L)
    x_fault_time = x_fault_time * q_time

    x = x_fault_time + x_stiffness_time + x_rotation_time

    SNR = 10 ** (SNR_dB / 10)
    Esym = np.sum(np.abs(x) ** 2) / L
    N0 = Esym / SNR
    noise_sigma = np.sqrt(N0)
    nt = noise_sigma * np.random.randn(1, L)
    x_noise = x + nt

    return t, x, x_noise, fr_time

# Example usage: TBD
d = 10
D = 50
contact_angle = np.deg2rad(15)
n = 8
fault_type = 'inner'
fc = 1000
fd = 50
fm = 25
fr = [100, 200, 300]
N = 1024
fs = 10000
SNR_dB = 30
q_fault = 0.1
q_stiffness = 0.2
q_rotation = 0.3

t, x, x_noise, fr_time = bearing_signal_model_dist(d, D, contact_angle, n, fault_type, fc, fd, fm, fr, N, fs, SNR_dB, q_fault, q_stiffness, q_rotation)