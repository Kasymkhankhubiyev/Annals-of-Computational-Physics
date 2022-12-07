import numpy as np
from typing import NamedTuple
from exceptions import CantMatchMethod
from CalcParticles import *

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


def _dichotomy_method(nc: float, nv: float, t: Kelvin, Jd: eV, Efpl: eV, Efneg: eV, Ec: eV,
                      Ev: eV, Nd: float, counter: int):
    a, b = Efpl, Efneg
    counter += 1

    Ef = (a + b) / 2.

    n = calc_n(nc=nc, Ef=Ef, Ec=Ec, t=t)
    p = calc_p(nv=nv, Ef=Ef, Ev=Ev, t=t)
    ndpl = calc_Ndplus(Nd=Nd, Ef=Ef, Ed=Ec - Jd, t=t)
    q = _count_Q(n=n, p=p, Nd=ndpl)

    if np.abs(q / (p + ndpl)) < 0.0001:
        print(f'Nd={Nd}     nc={nc}')
        return Ef, counter
    else:
        if q > 0:
            return _dichotomy_method(nc=nc, nv=nv, t=t, Jd=Jd, Ec=Ec, Ev=Ev, Nd=Nd, Efneg=Ef, Efpl=Efpl, counter=counter)
        elif q < 0:
            return _dichotomy_method(nc=nc, nv=nv, t=t, Jd=Jd, Ec=Ec, Ev=Ev, Nd=Nd, Efneg=Efneg, Efpl=Ef, counter=counter)


def _fixed_point_method():
    pass


def _newtown_method():
    pass


def _secant_method():
    pass


def calculate_fermi_level(me: me_effective, mh: mh_effective, t: Kelvin, Jd: eV, Efpl: eV, Efneg: eV, Ec: eV,
                      Nd: float, method: str, Ev = 0.) -> tuple:
    """

    :param me: эффективная масса электрона
    :param mh: эффективная масса дырки
    :param t: температура
    :param Jd: энергия ионизации
    :param Efpl: $E_{f}^{+}$ - нижняя граница
    :param Efneg: $E_{f}^{-}$ - верхняя граница
    :param Ec: дно зоны проводимости
    :param Ev: потолок валентной зоны
    :param Nd: концентрация доноров
    :param method: метод поиска уровня ферми
    :return:
    """
    methods = ['dichotomy', 'newtown', 'fixed-point', 'secant']
    fermi_level, counter = None, 0

    nc = calc_Nc(me, t)
    nv = calc_Nv(mh, t)
    try:
        if method == 'dichotomy':
            fermi_level, counter = _dichotomy_method(nc=nc, nv=nv, t=t, Jd=Jd, Efpl=Efpl, Efneg=Efneg,
                                                     Ec=Ec, Ev=Ev, Nd=Nd, counter=counter)
        elif method == 'newtown':
            pass
        elif method == 'fixed-point':
            pass
        elif method == 'secant':
            pass
        else:
            raise CantMatchMethod(message=method, methods=methods)

        print(f'Fermi level = {fermi_level}, \t needed {counter} iterations')
        return fermi_level, counter

    except CantMatchMethod as e:
        print(e.args)
