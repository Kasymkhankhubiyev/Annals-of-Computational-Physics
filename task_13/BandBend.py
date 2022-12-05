import numpy as np
from task_13.exceptions import CantMatchMethod, CantRunDichotomyMethod


def bend_function(epsilon: float, phi: float, Nd: float, Nas: float, Eas: float, Ef: float, t: float, Eout: float):
    """
    math: $\sqrt{\frac{\epsilon \phi_s N_d}{2 \pi e^2}} = N_as \frac{1}{1 + exp(\frac{Eas + \phi_s + E_f}{kT})} +
            \frac{E_{out}}{4\pi e}$
    """
    k = 1  # Boltzmann constant
    e = 1  # electron charge
    n_as = Nas / (1. + np.exp((Eas + phi + Ef) / (k * t))) + Eout / (4 * np.pi * e)
    w = (epsilon * phi * Nd / (2 * np.pi * e**2)) ** 0.5
    return w - n_as


def _dichotomy_method(epsilon: float, phi0: float, phi1: float, nd: float, n_as: float, e_as: float, e_f: float,
                      t: float, e_out: float, counter: int, tolerance=1e-7):
    return 0, 0


def _fixed_point_method(epsilon: float, phi: float, nd: float, n_as: float, e_as: float, e_f: float, t: float,
                        e_out: float, count: int,  tolerance=1e-7):
    pass


def _newtown_method(epsilon: float, phi: float, nd: float, n_as: float, e_as: float, e_f: float, t: float,
                    e_out: float, count: int, tolerance=1e-7):
    pass


def _secant_method(epsilon: float, phi: float, nd: float, n_as: float, e_as: float, e_f: float, t: float,
                   e_out: float, counter: int, tolerance=1e-7):
    pass


def calculate_band_bend(epsilon: float, Nd: float, t: float, Nas: float, Eas: float, Eout: float, method: str,
                        Ef: float, phi0: float, tolerance=1e-7, phi1=None) -> tuple:

    methods = ['dichotomy', 'newtown', 'fixed-point']
    bend, counter = None, 0

    try:
        if method == 'dichotomy':
            if phi1 is not None and phi1 > phi0:
                bend, counter = _dichotomy_method(epsilon=epsilon, phi0=phi0, phi1=phi1, nd=Nd, t=t, n_as=Nas, e_as=Eas,
                                                  e_f=Ef, e_out=Eout, tolerance=tolerance, counter=counter)
            else:
                raise CantRunDichotomyMethod(phi0=phi0, phi1=phi1)
        elif method == 'newtown':
            pass
        elif method == 'fixed-point':
            pass
        else:
            raise CantMatchMethod(message=method, methods=methods)

        print(f'Fermi level = {bend}, \t needed {counter} iterations')
        return bend, counter

    except CantMatchMethod as e:
        print(e.args)
    except CantRunDichotomyMethod as e:
        print(e.args)
