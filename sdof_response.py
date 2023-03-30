import numpy as np
import matplotlib.pyplot as     plt
def sdof_response(fs, k, zita, fn, Lsdof):
    """
    Acceleration of a SDOF system

    Input:
    fs : Sample frequency [Hz]
    k : Spring stiffness [N/m]
    zita : Damping coefficient
    fn : Natural frequency [Hz]
    Lsdof (int): Desired signal length [points]

    Output:
    np.array: Acceleration (row vector)
    """
    m = k / (2 * np.pi * fn) ** 2
    F = 1
    A = F / m
    omegan = 2 * np.pi * fn
    omegad = omegan * np.sqrt(1 - zita ** 2)

    t = np.arange(0, Lsdof) / fs
    # system response
    xt = A / omegad * np.exp(-zita * omegan * t) * np.sin(omegad * t)  # displacement
    xd = np.concatenate(([0], np.diff(xt) * fs))  # velocity
    sdof_resp_time = np.concatenate(([0], np.diff(xd) * fs))  # acceleration

    #Test Plot
    """
    plt.plot(t, xt)
    plt.xlabel('X Value')
    plt.ylabel('Y Value')
    plt.title('SDOF')
    plt.show()
    """

    return sdof_resp_time
