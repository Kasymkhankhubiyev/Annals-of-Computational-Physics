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

m, A, U0 = 1, 1, 1
h = 6.582 * 1e-16 #eV * s

alfa = 2 * m * pow(A, 2) * U0 / (h*h)
min_val = np.float64(1e-4)


def run():
    a = np.float64(0.99)
    b = np.float64(1) - min_val
    x = predictX(a, b)[-1]

    fhalf = []
    half = []

    print(f'Предположительный ответ: {x}')

    dichotomy_X = dichotomy_method(a, b, fhalf, half)
    print(f' Ответ методом дихотомии: {half[-1]}')

    iteration_X = simple_iter_method(a)
    print(f' Ответ методом интераций: {iteration_X[1][-1]}')

    Newtown_X= newtown_method(a)
    print(f' Ответ методом Ньютона:   {Newtown_X}')



def f(x):
    """
    ctg((2ma^2U(1 + eps)/h^2)^1/2) - (1/eps - 1)^1/2 = 0
    eps = -E/U
    alfa = 2ma^2U/h^2
    :param x: eps
    :return:
    """
    return 1/np.tan(np.sqrt(alfa * (1 - x))) - np.sqrt(1/x - 1)


def df(x):
    """
    2ma^2U/sin^2((2ma^2U(1 + eps)/h^2)^1/2) / (2(2ma^2U(1 + eps)/h^2)^1/2 + 1/(1 * (1/eps - 1)^1/2)eps^2))
    eps = -E/U
    alfa = 2ma^2U/h^2
    :param x: eps
    :return:
    """
    return alfa * 1 / pow(np.sin(np.sqrt(alfa * (1 - x))), 2) / (2 * np.sqrt(alfa * (1 - x))) + 1 / (
            2 * np.sqrt(1 / x - 1) * pow(x, 2))


def dichotomy_method(a, b, farr, xarr):
    num = len(farr)
    f_a = f(a)
    f_b = f(b)
    xarr.append((a + b) / 2)
    farr.append(f(xarr[num]))

    if np.abs(farr[num]) < min_val: #достигли точность
        return farr, xarr
    elif np.sign(f_a) == np.sign(farr[num]): #нет нулей
        return dichotomy_method(xarr[num], b, farr, xarr)
    elif np.sign(f_b) == np.sign(farr[num]):
        return dichotomy_method(a, xarr[num], farr, xarr) #Есть нуль


def simple_iter_method(X0):
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
    nextX = []
    Xvals = []
    lam = np.float64(0.0001)
    X1 = X0
    X0 -= lam * np.sign(df(X0)) * f(X0)
    while np.abs(X1 - X0) > min_val:
        X1 = X0
        X0 -= lam * np.sign(df(X0)) * f(X0)
        Xvals.append(X0)
        nextX.append(X1)
    return Xvals, nextX


def newtown_method(X0):
    """
    x_{n+1} = x_{n} - f(x_{n})/f'(x_{n})
    :return:
    """
    nextX = []
    Xvals = []
    X1 = X0
    div = 1. / df(X0)
    X0 = X0 - f(X0) * div
    Xvals.append(X0)
    nextX.append(X1)
    while np.abs(X1 - X0) > min_val:
        X1 = X0
        div = 1. / df(X0)
        X0 -= f(X0) * div
        Xvals.append(X0)
        nextX.append(X1)
    return nextX, Xvals


def predictX(a, b):
    x_arr = []
    for x in np.arange(a, b, min_val):
        if f(x) > 0:
            x_arr.append(x)
    return x_arr


def draw_graph(val0, val1, color, name, plt, y):
    plt.plot(val0, val1, color=color)
    plt.scatter(val0, val1, color='black')
    plt.set_title(name)
    plt.set_xlabel('x')
    plt.set_ylabel(y)
    plt.plot(val0[-1], val1[-1], '*', color='orange', markersize=10)

