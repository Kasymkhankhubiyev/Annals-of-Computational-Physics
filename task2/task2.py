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
"""

import matplotlib.pyplot as plt
import numpy as np
from typing import NamedTuple


class Result(NamedTuple):
    function: float
    energy: float
    iter_depth: int


def f(x: float, m: float, a: float, u0: float) -> float:
    """
    ctg((2ma^2U(1 + eps)/h^2)^1/2) - (1/eps - 1)^1/2 = 0
    eps = -E/U
    alfa = 2ma^2U/h^2
    """
    # h = 4.135e-015  # eV * s
    h = 1.

    alfa = 2 * m * a**2 * u0 / (h * h)  # 2ma^2U0/h^2 == const
    if np.tan(np.sqrt(alfa * (1 + x/u0))) == 0:
        return np.inf
    elif x == 0:
        return -np.inf
    else:
        return 1/np.tan(np.sqrt(alfa * (1 + x/u0))) - np.sqrt(-u0/x - 1)


def df(x: float, m: float, u0: float, a: float) -> float:
    # h = 6.582 * 1e-16  # eV * s
    h = 1.

    alfa = 2 * m * a**2 * u0 / (h * h)  # 2ma^2U0/h^2 == const
    sum1 = -1 * alfa / u0 / pow(np.sin(np.sqrt(alfa * (1 + x/u0))), 2) / (2 * np.sqrt(alfa * (1 + x/u0)))
    sum2 = - 0.5 * u0 / x**2 / np.sqrt(-u0/x - 1)
    return sum1 + sum2


def dichotomy_method(left: float, right: float, m: float, a: float, u0: float, tolerance: float, counter: int) -> Result:
    f_a = f(x=left, m=m, a=a, u0=u0)
    # f_b = f(x=right, m=m, a=a, u0=u0)
    x_i = (left + right) / 2
    f_i = f(x=x_i, m=m, a=a, u0=u0)

    if np.abs(f_i) < tolerance: #достигли точность
        return Result(function=f_i, energy=x_i, iter_depth=counter)
    elif f_a * f_i > 0: #нет нулей
        counter += 1
        return dichotomy_method(left=x_i, right=right, m=m, a=a, u0=u0, tolerance=tolerance, counter=counter)
    else:
        counter += 1
        return dichotomy_method(left=left, right=x_i, m=m, a=a, u0=u0, tolerance=tolerance, counter=counter)


def simple_iter_method(x0: float, m: float, a: float, u0: float, tolerance: float, counter :int) -> Result:
    """
    Метод простых итераций
    phi(x) := f(x) + x; f(x) --> phi(x) = x

    x_{n+1} = phi(x_{n}); |phi'(x)| < q < 1

    Чтобы схожимость была:
    x_{n+1} = x_{n} -lambda*f(x_n)

    sign(lambda):= sign(f'(x))

    верно для любой гладкой f(x)
    :param X0:
    :return:
    """
    diff = df(x=x0, m=m, a=a, u0=u0)
    _lambda = 1/diff  # * np.sign(diff)

    x = -_lambda * f(x=x0, m=m, a=a, u0=u0)

    if np.abs(x) <= tolerance:
        return Result(function=f(x=x0, m=m, a=a, u0=u0), energy=x0, iter_depth=counter)
    else:
        counter += 1
        return simple_iter_method(x0=x0+x, m=m, a=a, u0=u0, tolerance=tolerance, counter=counter)


def newtown_method(X0, tolerance):
    """
    x_{n+1} = x_{n} - f(x_{n})/f'(x_{n})
    :return:
    """
    nextX = []
    Xvals = []
    X1 = X0
    div = 1. / df(X0)
    X0 = X0 - f(X0) * div
    print(div * f(X0))
    Xvals.append(X0)
    nextX.append(X1)
    while np.abs(X1 - X0) > tolerance:
        X1 = X0
        div = 1. / df(X0)
        X0 -= f(X0) * div
        Xvals.append(X0)
        nextX.append(X1)
    return nextX, Xvals


def predictX(left: float, right: float, tolerance: float, a: float, m: float, u0: float) -> float:
    root = 0
    x, y = [], []
    for i in range(int((right-left)/tolerance)+1):
        if np.abs(f(x=(left + i * tolerance), m=m, a=a, u0=u0)) != np.inf:
            y.append(f(x=(left + i * tolerance), m=m, a=a, u0=u0))
            x.append((left + i * tolerance))
        if np.abs(f(x=(left + i * tolerance), m=m, a=a, u0=u0)) <= tolerance:
            root = left + i * tolerance

    plt.plot(x, y, color='green')
    # plt.plot(root, color='red')
    plt.ylim(-10, 10)

    plt.grid()
    plt.savefig('task2/original_function.png')
    return root


def draw_graph() -> None:
    pass


def run():
    # m, a, u0 = 0.5e6, np.float64(1), np.float64(1)  # ev, cm, eV
    m, a, u0 = 1., np.float64(1), np.float64(3)
    tolerance = np.float64(1e-12)

    left = np.float64(-u0)
    right = np.float64(0.0)

    predicted_root = predictX(left=left, right=right, tolerance=1e-4, a=a, m=m, u0=u0)
    print(predicted_root)

    # fhalf = []
    # half = []

    print(f'Предположительный ответ: {predicted_root}')

    counter = 0
    dichotomy = dichotomy_method(left=left, right=right, m=m, a=a, u0=u0, tolerance=tolerance, counter=counter)
    print(f' Ответ методом дихотомии: {dichotomy}')

    #за начальную точку возьмем значение из метода Дихотомии

    counter = 0
    iteration_X = simple_iter_method(x0=predicted_root, m=m, a=a, u0=u0, tolerance=tolerance, counter=counter)
    print(f' Ответ методом интераций: {iteration_X}')

    # Newtown_X= newtown_method(a)
    # print(f' Ответ методом Ньютона:   {Newtown_X}')
