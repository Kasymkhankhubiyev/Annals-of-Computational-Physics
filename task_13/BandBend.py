import numpy as np
from task_13.exceptions import CantMatchMethod

"""
math: $\sqrt{\frac{\epsilon \phi_s N_d}{2 \pi e^2}} = N_as \frac{1}{1 + exp(\frac{Eas + \phi_s + E_f}{kT})} + 
        \frac{E_{out}}{4\pi e}$
"""


def _dichotomy_method():
    return 0, 0


def _fixed_point_method():
    pass


def _newtown_method():
    pass


def _secant_method():
    pass


def calculate_band_bend(Eg: float, epsilon: float, mc: float, mv: float, Ed: float, Nd: float, t: float,
                        Nas: float, Eas: float, Eout: float, method: str):

    methods = ['dichotomy', 'newtown', 'fixed-point']
    bend, counter = None, 0

    try:
        if method == 'dichotomy':
            bend, counter = _dichotomy_method()
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
