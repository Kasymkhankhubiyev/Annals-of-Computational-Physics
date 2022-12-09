import matplotlib.pyplot as plt
import numpy as np
from typing import NamedTuple

from task_13.Fermi import calculate_fermi_level
from task_13.BandBend import calculate_band_bend


class Semiconductor:
    def __init__(self, me, mh, Jd, E_g, Nd, Nas, Eas, epsilon, E_f):
        self.me, self.mh, self.Jd, self.E_g = me, mh, Jd, E_g
        self.Nd, self.Nas, self.Eas = Nd, Nas, Eas,
        self.epsilon, self.E_f = epsilon, E_f


if __name__ == '__main__':

    Si = Semiconductor(me=0.36, mh=0.49, Jd=0.045, E_g=1.12, Nd=1e17, Nas=1e8, Eas=0.1, epsilon=11.7, E_f=1.)

    fermi_dichotomy = calculate_fermi_level(me=Si.me, mh=Si.mh, t=300, Jd=Si.Jd, Ef0=0.57, Ef1=1.12, Ec=Si.E_g,
                                            Nd=Si.Nd, method='dichotomy', tolerance=1e-7)

    fermi_fixed_point = calculate_fermi_level(me=Si.me, mh=Si.mh, t=300, Jd=Si.Jd, Ef0=1., Ec=Si.E_g,
                                              Nd=Si.Nd, method='fixed-point', tolerance=1e-7)

    fermi_newtown = calculate_fermi_level(me=Si.me, mh=Si.mh, t=300, Jd=Si.Jd, Ef0=1., Ec=Si.E_g,
                                          Nd=Si.Nd, method='newtown', tolerance=1e-7)

    fermi_secant = calculate_fermi_level(me=Si.me, mh=Si.mh, t=300, Jd=Si.Jd, Ef0=1., Ec=Si.E_g,
                                         Nd=Si.Nd, method='secant', delta=0.1, tolerance=1e-7)

    print('\nИзгиб зоны с уровнем Ферми полученным методом Дихотомии:\n')
    Si.E_f = fermi_dichotomy[0]

    phi_dichotomy_d = calculate_band_bend(epsilon=Si.epsilon, Nd=Si.Nd, t=300, Nas=Si.Nas, Eas=Si.Eas, Eout=0,
                                          Ef=Si.E_f, phi0=0, phi1=0.5, method='dichotomy', tolerance=1e-7)

    phi_fixed_point_d = calculate_band_bend(epsilon=Si.epsilon, Nd=Si.Nd, t=300, Nas=Si.Nas, Eas=Si.Eas, Eout=0,
                                            Ef=Si.E_f, phi0=0.05, method='fixed-point', tolerance=1e-7)

    phi_newtown_d = calculate_band_bend(epsilon=Si.epsilon, Nd=Si.Nd, t=300, Nas=Si.Nas, Eas=Si.Eas, Eout=0,
                                        Ef=Si.E_f, phi0=0.05, method='newtown', tolerance=1e-7)

    phi_secant_d = calculate_band_bend(epsilon=Si.epsilon, Nd=Si.Nd, t=300, Nas=Si.Nas, Eas=Si.Eas, Eout=0,
                                       Ef=Si.E_f, phi0=0.05, method='secant', tolerance=1e-7, delta=1e-2)

    print('\nИзгиб зоны с уровнем Ферми полученным методом Простых итераций:\n')
    Si.E_f = fermi_fixed_point[0]

    phi_dichotomy_fp = calculate_band_bend(epsilon=Si.epsilon, Nd=Si.Nd, t=300, Nas=Si.Nas, Eas=Si.Eas, Eout=0,
                                           Ef=Si.E_f, phi0=0, phi1=0.5, method='dichotomy', tolerance=1e-7)

    phi_fixed_point_fp = calculate_band_bend(epsilon=Si.epsilon, Nd=Si.Nd, t=300, Nas=Si.Nas, Eas=Si.Eas, Eout=0,
                                             Ef=Si.E_f, phi0=0.05, method='fixed-point', tolerance=1e-7)

    phi_newtown_fp = calculate_band_bend(epsilon=Si.epsilon, Nd=Si.Nd, t=300, Nas=Si.Nas, Eas=Si.Eas, Eout=0,
                                         Ef=Si.E_f, phi0=0.05, method='newtown', tolerance=1e-7)

    phi_secant_fp = calculate_band_bend(epsilon=Si.epsilon, Nd=Si.Nd, t=300, Nas=Si.Nas, Eas=Si.Eas, Eout=0,
                                        Ef=Si.E_f, phi0=0.05, method='secant', tolerance=1e-7, delta=1e-2)

    print('\nИзгиб зоны с уровнем Ферми полученным методом Ньютона:\n')
    Si.E_f = fermi_newtown[0]

    phi_dichotomy_n = calculate_band_bend(epsilon=Si.epsilon, Nd=Si.Nd, t=300, Nas=Si.Nas, Eas=Si.Eas, Eout=0,
                                          Ef=Si.E_f, phi0=0, phi1=0.5, method='dichotomy', tolerance=1e-7)

    phi_fixed_point_n = calculate_band_bend(epsilon=Si.epsilon, Nd=Si.Nd, t=300, Nas=Si.Nas, Eas=Si.Eas, Eout=0,
                                            Ef=Si.E_f, phi0=0.05, method='fixed-point', tolerance=1e-7)

    phi_newtown_n = calculate_band_bend(epsilon=Si.epsilon, Nd=Si.Nd, t=300, Nas=Si.Nas, Eas=Si.Eas, Eout=0,
                                        Ef=Si.E_f, phi0=0.05, method='newtown', tolerance=1e-7)

    phi_secant_n = calculate_band_bend(epsilon=Si.epsilon, Nd=Si.Nd, t=300, Nas=Si.Nas, Eas=Si.Eas, Eout=0,
                                       Ef=Si.E_f, phi0=0.05, method='secant', tolerance=1e-7, delta=1e-2)

    print('\nИзгиб зоны с уровнем Ферми полученным методом Секущих:\n')
    Si.E_f = fermi_secant[0]

    phi_dichotomy_s = calculate_band_bend(epsilon=Si.epsilon, Nd=Si.Nd, t=300, Nas=Si.Nas, Eas=Si.Eas, Eout=0,
                                          Ef=Si.E_f, phi0=0, phi1=0.5, method='dichotomy', tolerance=1e-7)

    phi_fixed_point_s = calculate_band_bend(epsilon=Si.epsilon, Nd=Si.Nd, t=300, Nas=Si.Nas, Eas=Si.Eas, Eout=0,
                                            Ef=Si.E_f, phi0=0.05, method='fixed-point', tolerance=1e-7)

    phi_newtown_s = calculate_band_bend(epsilon=Si.epsilon, Nd=Si.Nd, t=300, Nas=Si.Nas, Eas=Si.Eas, Eout=0,
                                        Ef=Si.E_f, phi0=0.05, method='newtown', tolerance=1e-7)

    phi_secant_s = calculate_band_bend(epsilon=Si.epsilon, Nd=Si.Nd, t=300, Nas=Si.Nas, Eas=Si.Eas, Eout=0,
                                       Ef=Si.E_f, phi0=0.05, method='secant', tolerance=1e-7, delta=1e-2)
