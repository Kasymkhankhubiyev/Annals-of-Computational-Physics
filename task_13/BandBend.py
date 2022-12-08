import numpy as np
from task_13.exceptions import CantMatchMethod, CantRunDichotomyMethod
# from fompy.constants import e


def _bend_function(epsilon: float, phi: float, Nd: float, Nas: float, Eas: float,
                   Ef: float, t: float, Eout: float) -> float:
    """
    :math: $\sqrt{\frac{\epsilon \phi_s N_d}{2 \pi e^2}} =
            N_as \frac{1}{1 + exp(\frac{Eas + \phi_s + E_f}{kT})} + \frac{E_{out}}{4\pi e}$
    """
    # eV = 1.60218e-12  # 1 eV = 1.60218e-12 arg
    k = 1.381e-16  # Boltzmann constant arg/K
    # e = 4.803e-10  # electron charge
    e, eV = 1, 1
    n_as = Nas * (1. / (1. + np.exp((Eas + phi - Ef) * eV / (k * 6.24e11 * t)))) + Eout * 3.3 * 1e-5 / (4 * np.pi * e)
    w = (epsilon * phi * eV * Nd / (2 * np.pi * e**2)) ** 0.5
    # print(w, n_as)
    return w - n_as


def _diff_funcrtion(epsilon: float, phi: float, Nd: float, Nas: float, Eas: float,
                    Ef: float, t: float) -> float:
    """
    returns a value of the derivative function in the given point.

    :math: $\frac{\epsilon N_d}{4 \pi e^2} \factor \frac{1}{\frac{\epsilon \phi_s N_d}{2 \pi e^2}} +
        N_as \frac{1}{(1 + exp(\frac{Eas + \phi_s + E_f}{kT}))^2} * exp(\frac{Eas + \phi_s + E_f}{kT})* \frac{1}{kT}$
    :param epsilon:
    :param phi:
    :param Nd:
    :param Nas:
    :param Eas:
    :param Ef:
    :param t:
    :param Eout:
    :return: dw - dN_as_plus
    """
    # eV = 1.60218e-12  # 1 eV = 1.60218e-12 arg
    k = 1.381e-16  # Boltzmann constant  arg/K
    e, eV = 1, 1  # electron charge
    diff_n_as = -1 * (Nas / (1. + np.exp((Eas + phi - Ef) * eV / (k * 6.24e11 * t)))**2) * \
                np.exp((Eas + phi * eV - Ef) / (k * 6.24e11 * t)) * (1 / (k * 6.24e11 * t))
    diff_w = 0.5 * epsilon * Nd / (2 * np.pi * e**2) / (epsilon * phi * eV * Nd / (2 * np.pi * e**2)) ** 0.5

    return diff_w - diff_n_as


def _dichotomy_method(epsilon: float, phi0: float, phi1: float, nd: float, n_as: float, e_as: float, e_f: float,
                      t: float, e_out: float, counter: int, tolerance=1e-7) -> tuple:
    counter += 1
    f_a = _bend_function(epsilon=epsilon, phi=phi0, Nd=nd, Nas=n_as, Eas=e_as, Ef=e_f, t=t, Eout=e_out)
    # f_b = f(x=right, m=m, a=a, u0=u0)
    phi_i = (phi0 + phi1) / 2
    f_i = _bend_function(epsilon=epsilon, phi=phi_i, Nd=nd, Nas=n_as, Eas=e_as, Ef=e_f, t=t, Eout=e_out)

    # print(f_a, f_i)
    if np.abs(f_i) < tolerance:  # достигли точность
        return phi_i, counter
    elif f_a * f_i > 0:  # нет нулей
        return _dichotomy_method(epsilon=epsilon, phi0=phi_i, phi1=phi1, nd=nd, n_as=n_as, e_as=e_as, e_f=e_f, t=t,
                                 e_out=e_out, counter=counter, tolerance=tolerance)
        # return _dichotomy_method(epsilon=epsilon, phi0=phi1, phi1=phi_i, nd=nd, n_as=n_as, e_as=e_as, e_f=e_f, t=t,
        #                          e_out=e_out, counter=counter, tolerance=tolerance)
    else:
        return _dichotomy_method(epsilon=epsilon, phi0=phi0, phi1=phi_i, nd=nd, n_as=n_as, e_as=e_as, e_f=e_f, t=t,
                                 e_out=e_out, counter=counter, tolerance=tolerance)
        # return _dichotomy_method(epsilon=epsilon, phi0=phi_i, phi1=phi1, nd=nd, n_as=n_as, e_as=e_as, e_f=e_f, t=t,
        #                          e_out=e_out, counter=counter, tolerance=tolerance)


def _fixed_point_method(epsilon: float, phi: float, phi_fixed: float, nd: float, n_as: float, e_as: float, e_f: float, t: float,
                        e_out: float, count: int,  tolerance=1e-7):
    diff = _diff_funcrtion(epsilon=epsilon, phi=phi_fixed, Nd=nd, Nas=n_as, Eas=e_as, Ef=e_f, t=t)
    # print(f'differencial: {diff}')
    _lambda = 1/diff
    # print((f'lambda: {_lambda}'))

    phi_i = -_lambda * _bend_function(epsilon=epsilon, phi=phi, Nd=nd, Nas=n_as, Eas=e_as, Ef=e_f, t=t, Eout=e_out)
    # print(f'phi_i: {phi_i}')

    if np.abs(_bend_function(epsilon=epsilon, phi=phi, Nd=nd, Nas=n_as, Eas=e_as, Ef=e_f, t=t, Eout=e_out)) <= tolerance:
        return phi, count
    else:
        return _fixed_point_method(epsilon=epsilon, phi=phi+phi_i, phi_fixed=phi_fixed, nd=nd, n_as=n_as, e_as=e_as,
                                   e_f=e_f, t=t, e_out=e_out, count=count+1, tolerance=tolerance)


def _newtown_method(epsilon: float, phi: float, nd: float, n_as: float, e_as: float, e_f: float, t: float,
                    e_out: float, count: int, tolerance=1e-7):
    _lambda = -_bend_function(epsilon=epsilon, phi=phi, Nd=nd, Nas=n_as, Eas=e_as, Ef=e_f, t=t, Eout=e_out) / \
              _diff_funcrtion(epsilon=epsilon, phi=phi, Nd=nd, Nas=n_as, Eas=e_as, Ef=e_f, t=t)

    if np.abs(_lambda) <= tolerance:
        return phi, count
    else:
        return _newtown_method(epsilon=epsilon, phi=phi+_lambda, nd=nd, n_as=n_as, e_as=e_as, e_f=e_f, t=t, e_out=e_out,
                               count=count + 1, tolerance=tolerance)


def _secant_method(epsilon: float, phi: float, nd: float, n_as: float, e_as: float, e_f: float, t: float,
                   e_out: float, counter: int, tolerance=1e-7):
    pass


def bend_methods() -> list:
    return ['dichotomy', 'newtown', 'fixed-point', 'secant']


def calculate_band_bend(epsilon: float, Nd: float, t: float, Nas: float, Eas: float, Eout: float, method: str,
                        Ef: float, phi0: float, tolerance=1e-7, phi1=None) -> tuple:

    methods = ['dichotomy', 'newtown', 'fixed-point', 'secant']
    bend, counter = None, 0

    try:
        if method == 'dichotomy':
            if phi1 is not None and phi1 > phi0:
                bend, counter = _dichotomy_method(epsilon=epsilon, phi0=phi0, phi1=phi1, nd=Nd, t=t, n_as=Nas, e_as=Eas,
                                                  e_f=Ef, e_out=Eout, tolerance=tolerance, counter=counter)
            else:
                raise CantRunDichotomyMethod(phi0=phi0, phi1=phi1)
        elif method == 'fixed-point':
            bend, counter = _fixed_point_method(epsilon=epsilon, phi=phi0, phi_fixed=phi0, nd=Nd, n_as=Nas,
                                                e_as=Eas, e_f=Ef, e_out=Eout, count=counter, tolerance=tolerance, t=t)
        elif method == 'newtown':
            bend, counter = _newtown_method(epsilon=epsilon, phi=phi0, nd=Nd, n_as=Nas, e_as=Eas, e_f=Ef, t=t,
                                            e_out=Eout, count=counter, tolerance=tolerance)
        elif method == 'secant':
            pass
        else:
            raise CantMatchMethod(message=method, methods=methods)

        print(f'bend width = {bend}, \t needed {counter} iterations')
        return bend, counter

    except CantMatchMethod as e:
        print(e.args)
    except CantRunDichotomyMethod as e:
        print(e.args)
