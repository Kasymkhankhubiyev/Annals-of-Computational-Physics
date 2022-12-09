import numpy as np
from task_13.exceptions import CantMatchMethod, CantRunDichotomyMethod, CatchZeroDelta


def _bend_function(epsilon: float, phi: float, Nd: float, Nas: float, Eas: float,
                   Ef: float, t: float, Eout: float) -> float:
    """
    :math: $\sqrt{\frac{\epsilon \phi_s N_d}{2 \pi e^2}} =
            N_as \frac{1}{1 + exp(\frac{Eas + \phi_s + E_f}{kT})} + \frac{E_{out}}{4\pi e}$

    :param epsilon: dielectric constant
    :param phi: band bend value in eV
    :param Nd: donors' concentration
    :param Nas: surface acceptors' concentration
    :param Eas: field created by surface acceptors
    :param Ef: Fermi level in eV
    :param t: temperature in Kelvin
    :param Eout: External field value in V/m
    :return: a value of the difference of different formulas for surface charge
    """
    k = 1.381e-16  # Boltzmann constant arg/K
    e, eV = 1, 1
    n_as = Nas * (1. / (1. + np.exp((Eas + phi - Ef) * eV / (k * 6.24e11 * t)))) + Eout * 3.3 * 1e-5 / (4 * np.pi * e)
    w = (epsilon * phi * eV * Nd / (2 * np.pi * e**2)) ** 0.5
    return w - n_as


def _diff_funcrtion(epsilon: float, phi: float, Nd: float, Nas: float, Eas: float,
                    Ef: float, t: float) -> float:
    """
    returns a value of the derivative function in the given point.

    :math: $\frac{\epsilon N_d}{4 \pi e^2} \factor \frac{1}{\frac{\epsilon \phi_s N_d}{2 \pi e^2}} +
        N_as \frac{1}{(1 + exp(\frac{Eas + \phi_s + E_f}{kT}))^2} * exp(\frac{Eas + \phi_s + E_f}{kT})* \frac{1}{kT}$
    :param epsilon: dielectric constant
    :param phi: band bend value in eV
    :param Nd: donors' concentration
    :param Nas: surface acceptors' concentration
    :param Eas: field created by surface acceptors
    :param Ef: Fermi level in eV
    :param t: temperature in Kelvin
    :param Eout: External field value in V/m
    :return: a value of the derivative function in the given point
    """
    k = 1.381e-16  # Boltzmann constant  arg/K
    e, eV = 1, 1  # electron charge
    diff_n_as = -1 * (Nas / (1. + np.exp((Eas + phi - Ef) * eV / (k * 6.24e11 * t)))**2) * \
                np.exp((Eas + phi * eV - Ef) / (k * 6.24e11 * t)) * (1 / (k * 6.24e11 * t))
    diff_w = 0.5 * epsilon * Nd / (2 * np.pi * e**2) / (epsilon * phi * eV * Nd / (2 * np.pi * e**2)) ** 0.5

    return diff_w - diff_n_as


def _dichotomy_method(epsilon: float, phi0: float, phi1: float, nd: float, n_as: float, e_as: float, e_f: float,
                      t: float, e_out: float, counter: int, tolerance=1e-7) -> tuple:
    """

    :param epsilon: dielectric constant
    :param phi0: lowest bend value in eV
    :param phi1: highest bend value in eV
    :param nd: donors' concentration
    :param n_as: surface acceptors' concentration
    :param e_as: field created by surface acceptors
    :param e_f: Fermi level in eV
    :param t: temperature in Kelvin
    :param e_out: External field value in V/m
    :param counter: counts iteration amount needed to calculate a solution
    :param tolerance: tolerant error
    :return: band bend value in eV and iterations steps amount
    """
    counter += 1
    f_a = _bend_function(epsilon=epsilon, phi=phi0, Nd=nd, Nas=n_as, Eas=e_as, Ef=e_f, t=t, Eout=e_out)
    phi_i = (phi0 + phi1) / 2
    f_i = _bend_function(epsilon=epsilon, phi=phi_i, Nd=nd, Nas=n_as, Eas=e_as, Ef=e_f, t=t, Eout=e_out)

    if np.abs(f_i) < tolerance:  # достигли точность
        return phi_i, counter
    elif f_a * f_i > 0:  # нет нулей
        return _dichotomy_method(epsilon=epsilon, phi0=phi_i, phi1=phi1, nd=nd, n_as=n_as, e_as=e_as, e_f=e_f, t=t,
                                 e_out=e_out, counter=counter, tolerance=tolerance)
    else:
        return _dichotomy_method(epsilon=epsilon, phi0=phi0, phi1=phi_i, nd=nd, n_as=n_as, e_as=e_as, e_f=e_f, t=t,
                                 e_out=e_out, counter=counter, tolerance=tolerance)


def _fixed_point_method(epsilon: float, phi: float, phi_fixed: float, nd: float, n_as: float, e_as: float, e_f: float, t: float,
                        e_out: float, count: int,  tolerance=1e-7):
    """
    :math: $\phi(x)=f(x)+x$ : $x_{n+1}=x_n - \frac{f(x_n)}{f'(\ksi)}$
    :param epsilon: dielectric constant
    :param phi: current bend value in eV
    :param phi_fixed: ksi - fixed starter value for a band bend in eV
    :param nd: donors' concentration
    :param n_as: surface acceptors' concentration
    :param e_as: field created by surface acceptors
    :param e_f: Fermi level in eV
    :param t: temperature in Kelvin
    :param e_out: external field value in V/m
    :param count: counts iteration amount needed to calculate a solution
    :param tolerance: tolerant error
    :return: band bend value in eV and iterations steps amount
    """

    if np.abs(_bend_function(epsilon=epsilon, phi=phi, Nd=nd, Nas=n_as, Eas=e_as, Ef=e_f, t=t, Eout=e_out)) <= tolerance:
        return phi, count
    else:
        diff = _diff_funcrtion(epsilon=epsilon, phi=phi_fixed, Nd=nd, Nas=n_as, Eas=e_as, Ef=e_f, t=t)
        _lambda = 1 / diff
        phi_i = -_lambda * _bend_function(epsilon=epsilon, phi=phi, Nd=nd, Nas=n_as, Eas=e_as, Ef=e_f, t=t, Eout=e_out)
        return _fixed_point_method(epsilon=epsilon, phi=phi+phi_i, phi_fixed=phi_fixed, nd=nd, n_as=n_as, e_as=e_as,
                                   e_f=e_f, t=t, e_out=e_out, count=count+1, tolerance=tolerance)


def _newtown_method(epsilon: float, phi: float, nd: float, n_as: float, e_as: float, e_f: float, t: float,
                    e_out: float, count: int, tolerance=1e-7):
    """
    :math: $x_{n+1} = x_n - \frac{f(x_n)}{f'(x_n)}$
    :param epsilon: dielectric constant
    :param phi: starter band bend value in eV
    :param nd: donors' concentration
    :param n_as: surface acceptors' concentration
    :param e_as: field created by surface acceptors
    :param e_f: Fermi level in eV
    :param t: temperature in Kelvin
    :param e_out: external field value in V/m
    :param count: counts iteration amount needed to calculate a solution
    :param tolerance: tolerant error
    :return: band bend value in eV and iterations steps amount
    """
    if np.abs(_bend_function(epsilon=epsilon, phi=phi, Nd=nd, Nas=n_as, Eas=e_as, Ef=e_f, t=t, Eout=e_out)) <= tolerance:
        return phi, count
    else:
        _lambda = -_bend_function(epsilon=epsilon, phi=phi, Nd=nd, Nas=n_as, Eas=e_as, Ef=e_f, t=t, Eout=e_out) / \
                  _diff_funcrtion(epsilon=epsilon, phi=phi, Nd=nd, Nas=n_as, Eas=e_as, Ef=e_f, t=t)
        return _newtown_method(epsilon=epsilon, phi=phi+_lambda, nd=nd, n_as=n_as, e_as=e_as, e_f=e_f, t=t, e_out=e_out,
                               count=count + 1, tolerance=tolerance)


def _secant_method(epsilon: float, phi_0: float, nd: float, n_as: float, e_as: float, e_f: float, t: float,
                   e_out: float, counter: int, tolerance=1e-7, delta=1e-3, phi_1=None):
    """
    Метод секущих требует два последовательных значения, для этого на первом этапе мы задаем delta - стартовая разница
    двух измерений. Далее получаем x_i+1 из x_i и x_{i-1}

    :math: $x_{n+1}=x_n - f(x_n)\frac{x_n-x_{n-1}{f(x_n)-f(x_{n-1})}$
    :param epsilon: dielectric constant
    :param phi_0: band bend value in eV for {n-1}-th iteration
    :param nd: donors' concentration
    :param n_as: surface acceptors' concentration
    :param e_as: field created by surface acceptors
    :param e_f: Fermi level in eV
    :param t: temperature in Kelvin
    :param e_out: external field value in V/m
    :param counter: counts iteration amount needed to calculate a solution
    :param tolerance: tolerant error
    :param delta: x_1 - x_0 difference for a first iteration
    :param phi_1: band bend value in eV for n-th iteration
    :return: band bend value in eV and iterations steps amount
    """
    if phi_1 is not None:
        if (_bend_function(epsilon=epsilon, phi=phi_1, Nd=nd, Nas=n_as, Eas=e_as,
                                         Ef=e_f, t=t, Eout=e_out)) < tolerance:
            return phi_1, counter
        else:
            factor = delta / (_bend_function(epsilon=epsilon, phi=phi_1,
                                             Nd=nd, Nas=n_as, Eas=e_as,
                                             Ef=e_f, t=t, Eout=e_out) - _bend_function(epsilon=epsilon, phi=phi_0,
                                                                                       Nd=nd, Nas=n_as, Eas=e_as,
                                                                                       Ef=e_f, t=t, Eout=e_out))
            phi_ = phi_1 - _bend_function(epsilon=epsilon,
                                          phi=phi_1, Nd=nd,
                                          Nas=n_as, Eas=e_as,
                                          Ef=e_f, t=t, Eout=e_out) * factor
            return _secant_method(epsilon=epsilon, phi_0=phi_0, nd=nd, n_as=n_as, e_as=e_as, e_f=e_f, t=t, e_out=e_out,
                                  counter=counter+1, tolerance=tolerance, phi_1=phi_)
    else:
        if delta != 0:
            factor = delta / (_bend_function(epsilon=epsilon, phi=phi_0,
                                             Nd=nd, Nas=n_as, Eas=e_as,
                                             Ef=e_f, t=t, Eout=e_out) - _bend_function(epsilon=epsilon,
                                                                                       phi=phi_0 - delta, Nd=nd,
                                                                                       Nas=n_as, Eas=e_as, Ef=e_f,
                                                                                       t=t, Eout=e_out))
            phi_ = phi_0 - _bend_function(epsilon=epsilon,
                                          phi=phi_0, Nd=nd,
                                          Nas=n_as, Eas=e_as,
                                          Ef=e_f, t=t, Eout=e_out) * factor
            return _secant_method(epsilon=epsilon, phi_0=phi_0, nd=nd, n_as=n_as, e_as=e_as, e_f=e_f, t=t, e_out=e_out,
                                  counter=counter+1, tolerance=tolerance, phi_1=phi_)
        else:
            raise CatchZeroDelta


def bend_methods() -> list:
    return ['dichotomy', 'newtown', 'fixed-point', 'secant']


def calculate_band_bend(epsilon: float, Nd: float, t: float, Nas: float, Eas: float, Eout: float, method: str,
                        Ef: float, phi0: float, tolerance=1e-7, phi1=None, delta=1e-3) -> tuple:

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
            bend, counter = _secant_method(epsilon=epsilon, phi_0=phi0, nd=Nd, n_as=Nas, e_as=Eas, e_f=Ef, t=t,
                                           e_out=Eout, counter=counter, tolerance=tolerance, phi_1=phi1, delta=delta)
        else:
            raise CantMatchMethod(message=method, methods=methods)

        print(f'method: [{method}]\t bend width = {bend}, \t needed {counter} iterations')
        return bend, counter

    except CantMatchMethod as e:
        print(e.args)
    except CantRunDichotomyMethod as e:
        print(e.args)
    except CatchZeroDelta as e:
        print('Delta equals zero')
