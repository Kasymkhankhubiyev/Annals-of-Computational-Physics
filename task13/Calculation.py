import numpy as np
from typing import NamedTuple

me_effective = float
mh_effective = float
Kelvin = float
Nc = float
Nd = float
eV = float


class Result(NamedTuple):
    Ef: float
    n: float
    p: float
    Ndpl: float
    Q: float
    ratio: float


def _count_Q(n: float, p: float, Nd: float) -> float:
    """
    :param n: кол-во негативных носителей
    :param p: кол-во положительных носителей
    :param Nd: кол-во ионизированных атомов
    :return: Q - заряд полупроводника
    """

    Q = n - p - Nd
    return Q


def calc_n(nc: float, Ef: float, Ec: float, t: Kelvin) -> float:  # Nparticle:
    """
    n = Nc * exp(- (Ec - Ef)/kT)
    """
    k = 1.38e-16  # эрг/К

    expl = np.exp((Ef - Ec)/(k * 6.24e11 * t))
    n = nc * expl
    return n


def calc_p(nv: float, Ef: float, Ev: float, t: Kelvin) -> float:  # Nparticle:
    """
    p = Nv * exp((Ev - Ef)/kT)
    """
    k = 1.38e-16  # эрг/К

    expl = np.exp((Ev - Ef) / (k * 6.24e11 * t))
    p = nv * expl
    return p


def calc_Ndplus(Nd: float, Ef: eV, Ed: eV, t: Kelvin):
    k = 1.38e-16  # эрг/К

    ndpl = Nd / (1. + 0.5 * np.exp((Ef - Ed) / (k * 6.24e11 * t)))
    return ndpl


def _calc_Nc(me: me_effective, t: Kelvin) -> float:  # Nparticle:
    k = 1.38e-023  # J/K
    h = 1.054e-034  # kg * m /sec^2
    m0 = 9.109e-031  # kg ~ 0.511MeV

    Nc = 2 * ((2 * np.pi * me * m0 * k * t)/((2 * np.pi * h)**2)) ** 1.5  # 1/m^3
    Nc /= 10**6  # 1/cm^3
    return Nc


def _calc_Nv(mh: mh_effective, t: Kelvin) -> float:  # Nparticle:
    k = 1.38e-023  # J/K
    h = 1.054e-034  # kg * m /sec^2
    m0 = 9.109e-031  # kg ~ 0.511MeV

    Nv = 2 * ((2 * np.pi * mh * m0 * k * t)/((2 * np.pi * h)**2)) ** 1.5  # 1/m^3
    Nv /= 10 ** 6  # 1/cm^3
    return Nv


def dihatomi_method():
    pass


def division_method():
    pass


def newtown_method():
    pass


def calculate_charges(me: me_effective, mh: mh_effective, t: Kelvin, Jd: eV, Efpl: eV, Efneg: eV, Ec: eV, Ev: eV, Nd: float, method: str):
    """
    Расчет делаем методом дихотомии
    :param me: эффективная масса плотности состояний электронов
    :param mh: эффективная масса плотности состояний дырок
    :param t: температура
    :param Efpl: уровень Ферми совпадает с дном зоны проводимости,
    все электроны в зоне проводимости, а дырки в валентной
    :param Efng: уровень Ферми близок к потолку валентной зоны
    :return:
    """

    # Jd = 0.05  # eV

    nc = _calc_Nc(me, t)
    nv = _calc_Nv(mh, t)

    a, b = Efpl, Efneg

    Ef = (a + b) / 2.

    n = calc_n(nc=nc, Ef=Ef, Ec=Ec, t=t)
    p = calc_p(nv=nv, Ef=Ef, Ev=Ev, t=t)
    ndpl = calc_Ndplus(Nd=Nd, Ef=Ef, Ed=Ec - Jd, t=t)
    q = _count_Q(n=n, p=p, Nd=ndpl)


    if np.abs(q/(p + ndpl)) < 0.0001:
        print(f'Nd={Nd}     nc={nc}')
        return Result(Ef=Ef, n=n, p=p, Ndpl=ndpl, Q=q, ratio=(q/(p + ndpl)))
    else:
        if q > 0:
            return calculate_charges(me=me, mh=mh, t=t, Jd=Jd, Ec=Ec, Ev=Ev, Nd=Nd, Efneg=Ef, Efpl=Efpl)
        elif q < 0:
            return calculate_charges(me=me, mh=mh, t=t, Jd=Jd, Ec=Ec, Ev=Ev, Nd=Nd, Efneg=Efneg, Efpl=Ef)