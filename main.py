import matplotlib.pyplot as plt
import numpy as np
from typing import NamedTuple

# from task_13.Fermi import calculate_fermi_level
from task_13.BandBend import calculate_band_bend


if __name__ == '__main__':
    class Semiconductor(NamedTuple):
        me: float
        mh: float
        Jd: float
        E_g: float
        Nd: float
        Nas: float
        Eas: float
        epsilon: float
        E_f: float

    Si = Semiconductor(me=0.36, mh=0.49, Jd=0.045, E_g=1.12, Nd=1e17, Nas=1e8, Eas=0.1, epsilon=11.7, E_f=1.)

    phi = calculate_band_bend(epsilon=Si.epsilon, Nd=Si.Nd, t=300, Nas=Si.Nas, Eas=Si.Eas, Eout=0,
                                       Ef=Si.E_f, phi0=0, phi1=0.5, method='dichotomy', tolerance=1e-7)

    phi = calculate_band_bend(epsilon=Si.epsilon, Nd=Si.Nd, t=300, Nas=Si.Nas, Eas=Si.Eas, Eout=0,
                              Ef=Si.E_f, phi0=0.05, method='fixed-point', tolerance=1e-7)

    print(11)
