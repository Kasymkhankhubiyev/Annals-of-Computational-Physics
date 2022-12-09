import numpy as np
from task_13.exceptions import CantMatchMethod, CatchZeroDelta, CantRunDichotomyMethod
from task_13.CalcParticles import *

me_effective = float
mh_effective = float
Kelvin = float
Nc = float
Nd = float
eV = float


def _count_Q(n: float, p: float, Nd: float) -> float:
    """
    :param n: кол-во негативных носителей
    :param p: кол-во положительных носителей
    :param Nd: кол-во ионизированных атомов
    :return: Q - заряд полупроводника
    """

    Q = n - p - Nd
    return Q


def _balance_function(nc: float, nv: float, nd: float, t: Kelvin, e_f: eV, e_c: eV, e_v: eV, e_d: eV) -> float:
    """
    :math: $Q=n - (N_d^+ + p)$
    :math: $Q = N_c \factor exp(\frac{E_f - E_c}{kT}) - N_v \factor exp(\frac{E_v-E_f}{kT}) - N_d \factor
            \frac{1}{1+0.5exp(\frac{E_f - E_d}{kT})}$
    :param nc: concentration of electrons
    :param nv: concentration of holes
    :param nd: concentration of donors
    :param t: temperature in Kelvin
    :param e_f: Fermi energy level in eV
    :param e_c: Conduction band energy level in eV
    :param e_v: Valence band energy level in eV
    :param e_d: Ionization energy in eV
    :return: difference in positively and negatively charged particles number
    """
    k = 1.38e-16  # эрг/К

    n = nc * np.exp((e_f - e_c) / (k * 6.24e11 * t))
    p = nv * np.exp((e_v - e_f) / (k * 6.24e11 * t))
    nd_plus = nd / (1. + 0.5 * np.exp((e_f - e_d) / (k * 6.24e11 * t)))
    Q = n - nd_plus - p
    return Q/(p + nd_plus)


def _diff_balance_function(nc: float, nv: float, nd: float, t: Kelvin, e_f: eV, e_c: eV, e_v: eV, e_d: eV) -> float:
    """
    $\frac{Q}{p + N_d^+}$
    :param nc: concentration of electrons
    :param nv: concentration of holes
    :param nd: concentration of donors
    :param t: temperature in Kelvin
    :param e_f: Fermi energy level in eV
    :param e_c: Conduction band energy level in eV
    :param e_v: Valence band energy level in eV
    :param e_d: Ionized Donors energy level in eV
    :return: differential of difference in positively and
             negatively charged particles number for a specific fermi level
    """
    k = 1.38e-16  # эрг/К
    n = nc * np.exp((e_f - e_c) / (k * 6.24e11 * t))
    p = nv * np.exp((e_v - e_f) / (k * 6.24e11 * t))
    nd_plus = nd / (1. + 0.5 * np.exp((e_f - e_d) / (k * 6.24e11 * t)))
    dn = nc * np.exp((e_f - e_c) / (k * 6.24e11 * t)) / (k * 6.24e11 * t)
    dp = -nv * np.exp((e_v - e_f) / (k * 6.24e11 * t)) / (k * 6.24e11 * t)
    dnd_plus = - nd / (1. + 0.5 * np.exp((e_f - e_d) / (k * 6.24e11 * t))) ** 2 * \
               0.5 * np.exp((e_f - e_d) / (k * 6.24e11 * t)) / (k * 6.24e11 * t)

    return (dn*(p + nd_plus) - (dp + dnd_plus)*n)/(p + nd_plus)**2


def _dichotomy_method(nc: float, nv: float, t: Kelvin, Jd: eV, Ef_low: eV, Ef_upper: eV, Ec: eV,
                      Ev: eV, Nd: float, counter: int, tolerance=1e-7):
    """

    :param nc: concentration of electrons
    :param nv: concentration of holes
    :param t: temperature in Kelvin
    :param Jd: ionization energy in eV
    :param Ef_low: lowest Fermi energy level in eV
    :param Ef_upper: upper Fermi energy level in eV
    :param Ec: conduction band energy level in eV
    :param Ev: valence band energy level in eV
    :param Nd: concentration of donors
    :param counter: counts iteration steps
    :param tolerance: acceptable error value
    :return: Fermi level  in eV and iterations steps amount
    """
    a, b = Ef_low, Ef_upper
    counter += 1

    Ef = (a + b) / 2.

    n = calc_n(nc=nc, Ef=Ef, Ec=Ec, t=t)
    p = calc_p(nv=nv, Ef=Ef, Ev=Ev, t=t)
    ndpl = calc_Ndplus(Nd=Nd, Ef=Ef, Ed=Ec - Jd, t=t)
    q = _count_Q(n=n, p=p, Nd=ndpl)

    if np.abs(q / (p + ndpl)) < tolerance:
        # print(f'Nd={Nd}     nc={nc}')
        return Ef, counter
    else:
        if q > 0:
            return _dichotomy_method(nc=nc, nv=nv, t=t, Jd=Jd, Ec=Ec, Ev=Ev, Nd=Nd, Ef_upper=Ef, Ef_low=Ef_low,
                                     counter=counter)
        elif q < 0:
            return _dichotomy_method(nc=nc, nv=nv, t=t, Jd=Jd, Ec=Ec, Ev=Ev, Nd=Nd, Ef_upper=Ef_upper, Ef_low=Ef,
                                     counter=counter)


def _fixed_point_method(nc: float, nv: float, nd: float, t: Kelvin, e_d: eV, e_f: eV, e_c: eV, e_v: eV,
                        e_f0: eV, counter: int, tolerance=1e-7) -> tuple:
    """
    :math: $\phi(x)=f(x)+x$ : $x_{n+1}=x_n - \frac{f(x_n)}{f'(\ksi)}$
    :param nc: concentration of electrons
    :param nv: concentration of holes
    :param nd: concentration of donors
    :param t: temperature in Kelvin
    :param e_d: ionization energy in eV
    :param e_f: current Fermi energy level in eV
    :param e_c: conduction band energy level in eV
    :param e_v: valence band energy level in eV
    :param e_f0: fixed level of Fermi energy in eV
    :param counter: counts iteration steps
    :param tolerance: acceptable error value
    :return: Fermi level  in eV and iterations steps amount
    """
    if np.abs(_balance_function(nc=nc, nv=nv, nd=nd, t=t, e_f=e_f, e_c=e_c, e_v=e_v, e_d=e_d)) < tolerance:
        return e_f, counter
    else:
        diff = _diff_balance_function(nc=nc, nv=nv, nd=nd, t=t, e_f=e_f0, e_c=e_c, e_v=e_v, e_d=e_d)
        _lambda = 1/diff
        ef_i = -_lambda * _balance_function(nc=nc, nv=nv, nd=nd, t=t, e_f=e_f, e_c=e_c, e_v=e_v, e_d=e_d)
        return _fixed_point_method(nc=nc, nv=nv, nd=nd, t=t, e_d=e_d, e_f=e_f+ef_i, e_c=e_c, e_v=e_v, counter=counter+1,
                                   e_f0=e_f0, tolerance=tolerance)


def _newtown_method(nc: float, nv: float, nd: float, t: Kelvin, e_d: eV, e_f: eV, e_c: eV, e_v: eV,
                    counter: int, tolerance=1e-7):
    """
    :math: $x_{n+1} = x_n - \frac{f(x_n)}{f'(x_n)}$
    :param nc:  concentration of electrons
    :param nv: concentration of holes
    :param nd: concentration of donors
    :param t: temperature in Kelvin
    :param e_d: ionization energy in eV
    :param e_f: current Fermi energy level in eV
    :param e_c: conduction band energy level in eV
    :param e_v: valence band energy level in eV
    :param counter: counts iteration steps
    :param tolerance: acceptable error value
    :return: Fermi level  in eV and iterations steps amount
    """
    if np.abs(_balance_function(nc=nc, nv=nv, nd=nd, t=t, e_f=e_f, e_c=e_c, e_v=e_v, e_d=e_d)) < tolerance:
        return e_f, counter
    else:
        diff = _diff_balance_function(nc=nc, nv=nv, nd=nd, t=t, e_f=e_f, e_c=e_c, e_v=e_v, e_d=e_d)
        _lambda = _balance_function(nc=nc, nv=nv, nd=nd, t=t, e_f=e_f, e_c=e_c, e_v=e_v, e_d=e_d) / diff
        return _newtown_method(nc=nc, nv=nv, nd=nd, t=t, e_d=e_d, e_f=e_f-_lambda, e_c=e_c, e_v=e_v, counter=counter+1,
                               tolerance=tolerance)


def _secant_method(nc: float, nv: float, nd: float, t: Kelvin, e_d: eV, e_f: eV, e_c: eV, e_v: float,
                   counter: int, tolerance=1e-7, delta=0.1, e_fi=None):
    """

    :param nc:  concentration of electrons
    :param nv: concentration of holes
    :param nd: concentration of donors
    :param t: temperature in Kelvin
    :param e_d: ionization energy in eV
    :param e_f: fermi level value in eV for {n-1}-th iteration
    :param e_c: conduction band energy level in eV
    :param e_v: valence band energy level in eV
    :param counter: counts iteration steps
    :param tolerance: Fermi level  in eV and iterations steps amount
    :param delta: x_1 - x_0 difference for a first iteration
    :param e_fi: fermi level value in eV for n-th iteration
    :return: Fermi level  in eV and iterations steps amount
    """
    if e_fi is not None:
        if (_balance_function(nc=nc, nv=nv, nd=nd, t=t, e_f=e_fi, e_c=e_c, e_v=e_v, e_d=e_d)) < tolerance:
            return e_fi, counter
        else:
            factor = delta / (_balance_function(nc=nc, nv=nv, nd=nd,
                                                t=t, e_f=e_fi, e_c=e_c,
                                                e_v=e_v, e_d=e_d) - _balance_function(nc=nc, nv=nv, nd=nd, t=t,
                                                                                      e_f=e_f, e_c=e_c, e_v=e_v,
                                                                                      e_d=e_d))

            ef_ = e_fi - _balance_function(nc=nc, nv=nv, nd=nd, t=t, e_f=e_fi, e_c=e_c, e_v=e_v, e_d=e_d) * factor
            return _secant_method(nc=nc, nv=nv, nd=nd, t=t, e_d=e_d, e_f=e_fi, e_c=e_c, e_v=e_v, counter=counter+1,
                                  tolerance=tolerance, e_fi=ef_)
    else:
        if delta != 0:
            factor = delta / (_balance_function(nc=nc, nv=nv, nd=nd,
                                                t=t, e_f=e_f, e_c=e_c,
                                                e_v=e_v, e_d=e_d) - _balance_function(nc=nc, nv=nv, nd=nd,
                                                                                      t=t, e_f=e_f-delta, e_c=e_c,
                                                                                      e_v=e_v, e_d=e_d))
            ef_ = e_f - _balance_function(nc=nc, nv=nv, nd=nd,
                                          t=t, e_f=e_f, e_c=e_c,
                                          e_v=e_v, e_d=e_d) * factor
            return _secant_method(nc=nc, nv=nv, nd=nd, t=t, e_f=e_f, e_fi=ef_, e_c=e_c, e_v=e_v, e_d=e_d,
                                  counter=counter+1, tolerance=tolerance)
        else:
            raise CatchZeroDelta


def fermi_methods() -> list:
    return ['dichotomy', 'newtown', 'fixed-point', 'secant']


def calculate_fermi_level(me: me_effective, mh: mh_effective, t: Kelvin, Jd: eV, Ef0: eV, Ec: eV,
                          Nd: float, method: str, Ev=0., tolerance=1e-7, Ef1=None, delta=1e-3) -> tuple:
    """
    :param me: эффективная масса электрона
    :param mh: эффективная масса дырки
    :param t: температура
    :param Jd: энергия ионизации
    :param Ef_low: $E_{f}^{+}$ - нижняя граница
    :param Ef_upper: $E_{f}^{-}$ - верхняя граница
    :param Ec: дно зоны проводимости
    :param Ev: потолок валентной зоны
    :param Nd: концентрация доноров
    :param method: метод поиска уровня ферми
    :return: Fermi level in eV and iterations steps amount
    """
    methods = ['dichotomy', 'newtown', 'fixed-point', 'secant']
    fermi_level, counter = None, 0

    nc = calc_Nc(me, t)
    nv = calc_Nv(mh, t)
    try:
        if method == 'dichotomy':
            if Ef1 is not None and Ef1 > Ef0:
                fermi_level, counter = _dichotomy_method(nc=nc, nv=nv, t=t, Jd=Jd, Ef_low=Ef0, Ef_upper=Ef1,
                                                         Ec=Ec, Ev=Ev, Nd=Nd, counter=counter, tolerance=tolerance)
            else:
                raise CantRunDichotomyMethod(phi0=Ef0, phi1=Ef1)
        elif method == 'fixed-point':
            fermi_level, counter = _fixed_point_method(nc=nc, nv=nv, nd=Nd, t=t, e_d=Ec-Jd, e_f=Ef0, e_f0=Ef0, e_c=Ec,
                                                       e_v=Ev, counter=counter, tolerance=tolerance)
        elif method == 'newtown':
            fermi_level, counter = _newtown_method(nc=nc, nv=nv, nd=Nd, t=t, e_d=Ec-Jd, e_f=Ef0, e_c=Ec, e_v=Ev,
                                                   counter=counter, tolerance=tolerance)
        elif method == 'secant':
            fermi_level, counter = _secant_method(nc=nc, nv=nv, nd=Nd, t=t, e_d=Ec-Jd, e_f=Ef0, e_c=Ec, e_v=Ev,
                                                  counter=counter, tolerance=tolerance, delta=delta)
        else:
            raise CantMatchMethod(message=method, methods=methods)

        print(f'method: [{method}]\t Fermi level = {fermi_level}, \t needed {counter} iterations')
        return fermi_level, counter

    except CantMatchMethod as e:
        print(e.args)
    except CantRunDichotomyMethod as e:
        print(e.args)
    except CatchZeroDelta as e:
        print('Delta equals zero')
