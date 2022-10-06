"""
Используя интегральное представление функции Бесселя целого индекса m

J_m(x) = 1/pi * Integrate(0, pi, cos(mt-xsint)dt)

и вычисляя производную с помощью конечной разности в тех же точках,
что и сам интеграл, продемонстрировать выполнение равенства

J'_0(x) + J_1(X) = 0

с точностью не хуже 1е-10 на отрезке [0, 2pi]

Для вычислений интегралов использовать метод трапеций и Симпсона
"""

import numpy as np
import matplotlib.pyplot as plt


def J(x, m, t):
    res = np.cos(m*t - x*np.sin(t))
    return res


def dJ0simpson(a, b, N, x, delta):
    res = (simpson_method(a, b, N, x+delta, 0) - simpson_method(a, b, N, x-delta, 0))/(2 * delta)
    return res

def dJ0trapeze(a, b, N, x, delta):
    res = (trapeze_method(a, b, N, x+delta, 0) - trapeze_method(a, b, N, x-delta, 0))/(2 * delta)
    return res


def simpson_method(a, b, N, x, m):
    """
    S[a,b] = (b - a)/6 * [f(a) +4f((a+b)/2) + f(b)]

    S[a,b] = h/3 [f_0_ + 4f_1_ + 2f_2_ + 4f_3_ + ... + 2f_2k-1_ + f_2k_]
    :return:
    """

    tarr = np.linspace(a, b, N)
    h = tarr[1] - tarr[0]
    yarr = [J(x, m, t) for t in tarr]
    S = 0
    for i in range(len(yarr)):
        if i % 2 == 0:
            factor = 2.
        else:
            factor = 4.
        S += factor * yarr[i]
    S = (S - yarr[0] - 3. * yarr[N - 1]) * h/3.

    return S


def trapeze_method(a, b, N, x, m):
    """
    S[a,b]=(f(a) + f(b))/2 * (b - a)
    :return:
    """

    tarr = np.linspace(a, b, N)
    h = tarr[1] - tarr[0]
    yarr = [J(x, m, t) for t in tarr]
    S = 0
    for y in yarr:
        S += y
    S = (S - 0.5*(yarr[0] + yarr[len(tarr)-1]))*h

    return S


def run():
    N = 100
    a = 0
    b = np.pi
    points = np.linspace(a, b, N)
    delta = points[1] - points[0]
    step = 2 * np.pi / N
    xarr, J1s,J0s, dJ0s = [], [], [], []
    J1t, J0t, dJ0t = [], [], []
    deltaJs, deltaJt = [], []
    for i in range(N + 1):
        x = i * step
        xarr.append(x)
        J0s.append(simpson_method(a, b, N, x, 0))
        J1s.append(simpson_method(a, b, N, x, 1))
        dJ0s.append(dJ0simpson(a, b, N, x, delta))
        deltaJs.append(dJ0s[i] + J1s[i])
        J0t.append(trapeze_method(a, b, N, x, 0))
        J1t.append(trapeze_method(a, b, N, x, 0))
        dJ0t.append(dJ0trapeze(a, b, N, x, delta))
        deltaJt.append(dJ0t[i] + J1t[i])

    fig, axs = plt.subplots(nrows=2, ncols=3)


    colors = ['green', 'red', 'blue']

    axs[0, 0].plot(xarr, J1s, label='J1 simps', color=colors[0])
    axs[0, 0].plot(xarr, J0s, label='J0 simps', color=colors[1])
    axs[0, 0].plot(xarr, dJ0s, label='dJ0 simps', color=colors[2])
    axs[0, 0].set(xlabel='X')
    axs[0, 0].set(ylabel='Jm(x)')
    axs[0, 0].set(title='Ф-ции Бесселя')

    axs[1, 0].plot(xarr, J1t, label='J1 trapeze', color=colors[0])
    axs[1, 0].plot(xarr, J0t, label='J0 trapeze', color=colors[1])
    axs[1, 0].plot(xarr, dJ0t, label='dJ0 trapeze', color=colors[2])
    axs[1, 0].set(xlabel='X')
    axs[1, 0].set(ylabel='Jm(x)')
    axs[1, 0].set(title='Ф-ции Бесселя')

    axs[0, 1].plot(xarr, J1s, label='J1 simps', color=colors[0])
    axs[0, 1].plot(xarr, dJ0s, label='dJ0 simps', color=colors[2])
    axs[0, 1].set(title='J1 and dJ0 simps')

    axs[1, 1].plot(xarr, J1t, label='J1 trapeze', color=colors[0])
    axs[1, 1].plot(xarr, dJ0t, label='dJ0 trapeze', color=colors[2])
    axs[1, 1].set(title='J1 and dJ0 trapeze')

    axs[0, 2].plot(xarr, deltaJs, label='Simpson', color=colors[1])
    axs[0, 2].set(title="J'0(x) + J1(x)")

    axs[1, 2].plot(xarr, deltaJt, label='Trapeze', color=colors[1])
    axs[1, 2].set(title="J'0(x) + J1(x)")

    plt.legend(loc='best')
    plt.savefig('task4.png')


