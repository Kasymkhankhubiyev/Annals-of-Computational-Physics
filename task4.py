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
    res = np.float64(np.cos(np.float64(m)*np.float64(t) - np.float64(x)*np.sin(np.float64(t))))
    return res


def dJ0simpson(a, b, N, x, delta):
    factor = np.float64(1000.)
    j0l = simpson_method(a, b, N, x+(np.float64(delta)/factor), 0)
    j0r = simpson_method(a, b, N, x-(np.float64(delta)/factor), 0)
    res = np.float64((j0l - j0r)/(np.float64(2.) * (np.float64(delta)/factor)))
    return res


def dJ0trapeze(a, b, N, x, delta):
    factor = np.float64(1000.)
    j0r = trapeze_method(a, b, N, x+(np.float64(delta)/factor), 0)
    j0l = trapeze_method(a, b, N, x-(np.float64(delta)/factor), 0)
    res = np.float64((j0r - j0l)/(np.float64(2.) * (np.float64(delta)/factor)))
    return res


def simpson_method(a, b, N, x, m):
    """
    S[a,b] = (b - a)/6 * [f(a) +4f((a+b)/2) + f(b)]

    S[a,b] = h/3 [f_0_ + 4f_1_ + 2f_2_ + 4f_3_ + ... + 2f_2k-1_ + f_2k_]
    :return:
    """

    tarr = np.linspace(np.float64(a), np.float64(b), 2**9)
    h = tarr[1] - tarr[0]
    yarr = [J(x, m, t) for t in tarr]
    S = np.float64(0.)
    for i in range(len(yarr)):
        if i % 2 == 0:
            factor = np.float64(2.)
        else:
            factor = np.float64(4.)
        S += np.float64(factor * yarr[i])
    S = np.float64((S - yarr[0] - np.float64(3.) * yarr[len(tarr)-1]) * h/np.float64(3.))

    return S


def trapeze_method(a, b, N, x, m):
    """
    S[a,b]=(f(a) + f(b))/2 * (b - a)
    :return:
    """

    tarr = np.linspace(a, b, 2**9)
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
        J1t.append(trapeze_method(a, b, N, x, 1))
        dJ0t.append(dJ0trapeze(a, b, N, x, delta))
        deltaJt.append(dJ0t[i] + J1t[i])

    fig, axs = plt.subplots(nrows=1, ncols=3)
    plt.tight_layout(pad=3, h_pad=3, w_pad=3)


    colors = ['green', 'red', 'blue']

    axs[0].plot(xarr, J1s, label='J1 simps', color=colors[0])
    axs[0].plot(xarr, J0s, label='J0 simps', color=colors[1])
    axs[0].plot(xarr, dJ0s, label='dJ0 simps', color=colors[2])
    axs[0].set(xlabel='X')
    axs[0].set(ylabel='Jm(x)')
    axs[0].set(title='Ф-ции Бесселя')
    axs[0].legend(fontsize=7,
                  ncol=1,
                  facecolor='oldlace',
                  edgecolor='r')

    axs[1].plot(xarr, J1s, label='J1 simps', color=colors[0])
    axs[1].plot(xarr, dJ0s, label='dJ0 simps', color=colors[2])
    axs[1].set(title='J1 and dJ0 simps')
    axs[1].legend(fontsize=7,
                  ncol=1,
                  facecolor='oldlace',
                  edgecolor='r')

    axs[2].plot(xarr, deltaJs, label='Simpson', color=colors[1])
    #axs[2].set(title="      J'0(x) + J1(x)")

    fig.set_size_inches(11.5, 6.5)
    plt.legend(loc='best')
    plt.savefig('task4simps.png')

    plt.close()
    fig, axs = plt.subplots(nrows=1, ncols=3)
    plt.tight_layout(pad=3, h_pad=3, w_pad=3)

    axs[0].plot(xarr, J1t, label='J1 trapeze', color=colors[0])
    axs[0].plot(xarr, J0t, label='J0 trapeze', color=colors[1])
    axs[0].plot(xarr, dJ0t, label='dJ0 trapeze', color=colors[2])
    axs[0].set(xlabel='X')
    axs[0].set(ylabel='Jm(x)')
    axs[0].set(title='Ф-ции Бесселя')
    axs[0].legend(fontsize=7,
                  ncol=1,
                  facecolor='oldlace',
                  edgecolor='r')

    axs[1].plot(xarr, J1t, label='J1 trapeze', color=colors[0])
    axs[1].plot(xarr, dJ0t, label='dJ0 trapeze', color=colors[2])
    axs[1].set(title='J1 and dJ0 trapeze')
    axs[1].legend(fontsize=7,
                  ncol=1,
                  facecolor='oldlace',
                  edgecolor='r')

    axs[2].plot(xarr, deltaJt, label='Trapeze', color=colors[1])
    axs[2].set(title="J'0(x) + J1(x)")

    fig.set_size_inches(11.5, 6.5)
    plt.legend(loc='best')
    plt.savefig('task4trapeze.png')


