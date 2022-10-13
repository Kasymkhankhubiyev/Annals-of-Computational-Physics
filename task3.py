"""
Вычислить инетегралы методами трапеций и Симпсона,
разделив отрезок интегрирования на 4, 8, 16, ... 2^n интервалов

Как убывает погрешность численного интегрирования с ростом числа интервалов?
"""
import numpy as np
import math
import matplotlib.pyplot as plt
from scipy.integrate import quad
a1, b1, a2, b2 = 1., -1., 0., 1.
N = 50


def f1(x) -> float:
    y = 1. / (1. + x * x)
    return y


def f2(x) -> float:
    y = x ** (1./3.) * np.exp(np.sin(x))
    return y


def trapeze_method(a, b, N, f, I) -> list:
    """
    S[a,b]=(f(a) + f(b))/2 * (b - a)
    :return:
    """
    errarr = []

    for n in range(2, N+1, 1):
        xarr = np.linspace(a, b, n)
        h = (b - a) / (n)  # xarr[1] - xarr[0]
        yarr = [f(x) for x in xarr]
        S = 0
        for y in yarr:
            S += y
        S = (S - 0.5*(yarr[0] + yarr[len(xarr)-1]))*h
        errarr.append(np.abs((I-S)/I))
        print(f' Кол-во разбиений:\t {n} \t\t Значение интеграла:\t{S}')

    return errarr


def trapeze(a, b, N, f, I) -> list:

    err = []
    for n in range(2, N+1, 1):
        x = np.linspace(a, b, n**2 + 1)
        y = np.array([f(i) for i in x])#np.vectorize(f(x))

        s = (b - a) / (2 * (n**2)) * np.sum(y[1:] + y[:-1])
        err.append(np.abs(I-s)/I)

    return err


def simpson(a, b, N, f, I) -> list:

    err = []
    for n in range(1, N+1, 1):
        x = np.linspace(a, b, n**2 + 1)
        y = np.array([f(i) for i in x])  # np.vectorize(f(x))

        s = (b - a)/((n ** 2) * 3) * np.sum(y[0:-1:2] + 4*y[1::2] + y[2::2])
        err.append(np.abs(I-s)/I)

    return err


def simpson_method(a, b, N, f, I) -> list:
    """
    S[a,b] = (b - a)/6 * [f(a) +4f((a+b)/2) + f(b)]

    S[a,b] = h/3 [f_0_ + 4f_1_ + 2f_2_ + 4f_3_ + ... + 2f_2k-1_ + f_2k_]
    :return:
    """

    errarr = []
    for n in range(2, N+1, 1):
        xarr = np.linspace(a, b, 2 * n)
        h = (b - a) / (2 * n)  #xarr[1] - xarr[0]
        yarr = [f(x) for x in xarr]
        S = 0
        for i in range(1, len(yarr)-1, 1):
            if i % 2 == 0:
                    factor = 2.
            else:
                    factor = 4.
            S = S + factor * yarr[i]
        #S = (S - yarr[0] - 3. * yarr[2 ** n - 1]) * h/3.
        S = (S + yarr[0] + yarr[len(yarr) - 1])*h/3.
        errarr.append(np.abs((I-S)/I))
        print(f' Кол-во разбиений:\t {n} \t\t Значение интеграла:\t{S}')

    return errarr


def run():
    Integral1 = quad(func=f1, a=a1, b=b1)[0]
    Integral2 = quad(func=f2, a=a2, b=b2)[0]
    print(Integral1)
    print(Integral2)
    print('Считаем первый интеграл:')
    print('trapeze method: \n')
    x = [n for n in range(2, N+1, 1)]
    y1arr = trapeze_method(a1, b1, N, f1, Integral1)
    # y1arr = trapeze(a1, b1, N, f1, Integral1)
    print('\n')
    print('Simpson method: \n')
    colors = ['blue', 'red']
    y2arr = simpson_method(a1, b1, N, f1, Integral1)
    # y2arr = simpson(a1, b1, N, f1, Integral1)
    #plt.ylim(0, 0.1)
    plt.plot(x, y1arr, label='трапеции', color=colors[0])
    plt.plot(x, y2arr, label='Симпсон', color=colors[1])
    plt.xlabel('Кол-во делений')
    plt.ylabel('Относительная ошибка')
    plt.title('1/(1+x^2)')
    plt.legend(loc='best')
    plt.savefig('page1.png')
    print('\n\nСчитаем второй интеграл:')
    print('trapeze method: \n')
    y1arr = trapeze_method(a2, b2, N, f2, Integral2)
    # y1arr = trapeze(a2, b2, N, f2, Integral2)
    print('\n')
    print('Simpson method: \n')
    y2arr = simpson_method(a2, b2, N, f2, Integral2)
    # y2arr = simpson(a2, b2, N, f2, Integral2)
    plt.close()
    #plt.ylim(0, 0.1)
    plt.plot(x, y1arr, label='трапеции', color=colors[0])
    plt.plot(x, y2arr, label='Симпсон', color=colors[1])
    plt.xlabel('Кол-во делений')
    plt.ylabel('Относительная ошибка')
    plt.legend(loc='best')
    plt.title('x^(1/3)*exp(sin(x))')
    plt.savefig('page2.png')
