import numpy as np
from task13.exceptions import CantMatchMethod


def _dichotomy_method():
    return 0, 0


def _fixed_point_method():
    pass


def _newtown_method():
    pass


def _secant_method():
    pass


def calculate_band_bend(method: str):
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
