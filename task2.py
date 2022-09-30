"""
Используя метод дихотомии, простых итераций и Ньютона,
найти уровень энергии Е основного состояния квантовой
частицы в прямоугольной потенциальной яме:

Потенциал:

U(x)={ -U, |x| <= a; 0, |x| > a}

Уравнение:

-h^2/2m * phi''(x) + U(x)*phi(x) = E*phi(x)

Решение:

phi(x) = {Acos(kx), |x| < a; B*e^(-qx), |x| > a}

k^2 = 2m/h^2(U + E), q^2 = -2m/h^2 * E

Трансцендентное уравнение:

ctg((2ma^2U(1 + eps)/h^2)^1/2) - (1/eps - 1)^1/2 = 0

Метод дихотомии - деление отрезка пополам

f in [a, b]; if f(a)*f(1/2(a+b)) <= 0 --> f in [a, 1/2(a+b)];
            else f in [1/2(a+b), b]

            ai - bi <= n - точность

Метод простых итераций
"""

import math
import numpy as np

m, a, U0 = 1, 1, 1

alfa = 2 * m * pow(a, 2) * U0
tolerance = 1e-4
fhalf, half = [], []


def run():
    pass


def f(x):
    """
    ctg((2ma^2U(1 + eps)/h^2)^1/2) - (1/eps - 1)^1/2 = 0
    eps = -E/U
    alfa = 2ma^2U
    h = 1
    :param x: eps
    :return:
    """
    return 1/np.tan(np.sqrt(alfa * (1 - x))) - np.sqrt(1/x - 1)


def df(x):
    """
    2ma^2U/sin^2((2ma^2U(1 + eps)/h^2)^1/2) / (2(2ma^2U(1 + eps)/h^2)^1/2 + 1/(1 * (1/eps - 1)^1/2)eps^2))
    eps = -E/U
    alfa = 2ma^2U
    h = 1
    :param x: eps
    :return:
    """
    return alfa * 1 / pow(np.sin(np.sqrt(alfa * (1 - x))), 2) / (2 * np.sqrt(alfa * (1 - x))) + 1 / (
            2 * np.sqrt(1 / x - 1) * pow(x, 2))


def dichotomy_method(a, b):
    num = len(half)
    f_a = f(a)
    f_b = f(b)
    half.append((a + b) / 2)
    fhalf.append(f(half[num]))

    if np.abs(fhalf[num]) < tolerance: #достигли точность
        return fhalf, half
    elif np.sign(f_a) == np.sign(fhalf[num]): #нет нулей
        return dichotomy_method(half[num], b)
    elif np.sign(f_b) == np.sign(fhalf[num]):
        return dichotomy_method(a, half[num]) #Есть нуль


def simple_iter_method():
    pass


def newtown_method():
    pass
