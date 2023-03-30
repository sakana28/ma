import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d
from sdof_response import sdof_response

def bearing_signal_model_local(d, D, contact_angle, n, fault_type, fr, fc, fd, fm, N, variance_factor, fs, k, zita, fn, Lsdof, SNR_dB, q_amp_mod=1):
    """
    Generation of a simulated signal for localized fault in rolling element bearing
    
    Parameters:
    d (float): Bearing roller diameter [mm]
    D (float): Pitch circle diameter [mm]
    contact_angle (float): Contact angle [rad]
    n (int): Number of rolling elements
    fault_type (str): Fault type selection: inner, outer, ball [string]
    fr (np.array): Row vector containing the rotation frequency profile
    fc (np.array): Row vector containing the carrier component of the speed
    fd (np.array): Row vector containing the frequency deviation
    fm (np.array): Row vector containing the modulation frequency
    N (int): Number of points per revolution
    variance_factor (float): Variance for the generation of the random contribution (ex. 0.04)
    fs (float): Sample frequency of the time vector
    k (float): SDOF spring stiffness [N/m]
    zita (float): SDOF damping coefficient
    fn (float): SDOF natural frequency [Hz]
    Lsdof (int): Length of the in number of points of the SDOF response
    SNR_dB (float): Signal to noise ratio [dB]
    q_amp_mod (float): Amplitude modulation due to the load (ex. 0.3, default: 1)

    Returns:
     t = time signal [s]
     x = simulated bearing signal without noise
     xNoise = simulated bearing signal with noise
     frTime = speed profile in the time domain [Hz]
     meanDeltaT = theoretical mean of the inter-arrival times
     varDeltaT = theoretical variance of the inter-arrival times
     menDeltaTimpOver = real mean of the inter-arrival times
     varDeltaTimpOver = real variance of the inter-arrival times
     errorDeltaTimp = generated error in the inter-arrival times

    """

    if fault_type == 'inner':
        geometry_parameter = 1 / 2 * (1 + d/D * np.cos(contact_angle))  # inner race fault
    elif fault_type == 'outer':
        geometry_parameter = 1 / 2 * (1 - d/D * np.cos(contact_angle))  # outer race fault
    elif fault_type == 'ball':
        geometry_parameter = 1 / (2 * n) * (1 - (d/D * np.cos(contact_angle))**2) / (d/D)  # outer race fault

    Ltheta = len(fr)
    theta = np.arange(0, Ltheta) * 2 * np.pi / N

    delta_theta_fault = 2 * np.pi / (n * geometry_parameter)
    number_of_impulses = int(np.floor(theta[-1] / delta_theta_fault))
    mean_delta_theta = delta_theta_fault
    var_delta_theta = (variance_factor * mean_delta_theta) ** 2
    delta_theta_fault = np.sqrt(var_delta_theta) * np.random.randn(1, number_of_impulses - 1) + mean_delta_theta
    theta_fault = np.concatenate(([0], np.cumsum(delta_theta_fault)))
    fr_theta_fault = np.interp(theta_fault, theta, fr) #change all interp to interp1d

    delta_t_imp = delta_theta_fault / (2 * np.pi * fr_theta_fault[1:])
    t_t_imp = np.concatenate(([0], np.cumsum(delta_t_imp)))

    L = int(np.floor(t_t_imp[-1] * fs))  # signal length
    t = np.arange(0, L) / fs
    fr_time = np.interp(t, t_t_imp, fr_theta_fault)

    delta_t_imp_index = np.round(delta_t_imp * fs).astype(int)
    error_delta_t_imp = delta_t_imp_index / fs - delta_t_imp

    index_impulses = np.concatenate(([1], np.cumsum(delta_t_imp_index)))
    index = len(index_impulses)
    while index_impulses[index - 1] / fs > t[-1]:
        index -= 1
    index_impulses = index_impulses[:index]

    mean_delta_t = np.mean(delta_t_imp)
    var_delta_t = np.var(delta_t_imp)
    mean_delta_t_imp_over = np.mean(delta_t_imp_index / fs)
    var_delta_t_imp_over = np.var(delta_t_imp_index / fs)

    x = np.zeros(L)
    x[index_impulses] = 1

    # amplitude modulation
    if fault_type == 'inner':

        if len(fc) > 1:
            theta_time = np.zeros(len(fr))
            for index in range(1, len(fr)):
                theta_time[index] = theta_time[index - 1] + (2 * np.pi / N) / (2 * np.pi * fr[index])
            fc_time = np.interp(t, theta_time, fc)
            fd_time = np.interp(t, theta_time, fd)
            fm_time = np.interp(t, theta_time, fm)

            q = 1 + q_amp_mod * np.cos(2 * np.pi * fc_time * t + 2 * np.pi * fd_time * np.cumsum(np.cos(2 * np.pi * fm_time * t) / fs))
        else:
            q = 1 + q_amp_mod * np.cos(2 * np.pi * fc * t + 2 * np.pi * fd * np.cumsum(np.cos(2 * np.pi * fm * t) / fs))
        x = q * x

    sdof_resp_time = sdof_response(fs, k, zita, fn, Lsdof)
    x = np.convolve(sdof_resp_time, x, mode='same')

    L = len(x)
    np.random.seed(0)  # set the random generator seed to default (for comparison only)
    SNR = 10 ** (SNR_dB / 10)  # SNR to linear scale
    E_sym = np.sum(np.abs(x) ** 2) / L  # Calculate actual symbol energy
    N0 = E_sym / SNR  # Find the noise spectral density
    noise_sigma = np.sqrt(N0)  # Standard deviation for AWGN Noise when x is real
    nt = noise_sigma * np.random.randn(1, L)  # computed noise
    x_noise = x + nt.flatten()  # received signal

    return t, x, x_noise, fr_time, mean_delta_t, var_delta_t, mean_delta_t_imp_over, var_delta_t_imp_over, error_delta_t_imp

# Example usage

d = 21.4
D = 203
contact_angle = 9*np.pi/180
n = 23
fault_type = 'inner'

N = 2048  
Ltheta = 10000 * N  
theta = np.arange(Ltheta) * 2 * np.pi / N
fc = np.array([10])
fd = 0.08 * fc
fm = 0.1 * fc
fr = fc + 2 * np.pi * fd * (np.cumsum(np.cos(fm * theta))/N)
variance_factor = 0.04
fs = 20000
k =  2e13
zita = 0.05
fn = 6e3
Lsdof = 2**8
SNR_dB = 0
q_amp_mod = 0.3

t, x, x_noise, fr_time, mean_delta_t, var_delta_t, mean_delta_t_imp_over, var_delta_t_imp_over, error_delta_t_imp = bearing_signal_model_local(
    d, D, contact_angle, n, fault_type, fr, fc, fd, fm, N, variance_factor, fs, k, zita, fn, Lsdof, SNR_dB, q_amp_mod)
print(x)
plt.figure().set_figwidth(15)
plt.plot(t, x)
plt.xlabel('time')
plt.ylabel('amplitude')
plt.title('without noise')
plt.show()




