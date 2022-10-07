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
N = 5



def run():
    Integral1 = quad(func=f1, a=a1, b=b1)[0]
    Integral2 = quad(func=f2, a=a2, b=b2)[0]
    print(Integral2)
    print('Считаем первый интеграл:')
    print('trapeze method: \n')
    x = [2**n for n in range(1, N)]
    y1arr = trapeze_method(a1, b1, N, f1, Integral1)
    print('\n')
    print('Simpson method: \n')
    colors = ['blue', 'red']
    y2arr = simpson_method(a1, b1, N, f1, Integral1)
    #plt.ylim(0, 0.1)
    plt.plot(x, y1arr, label='трапеции', color=colors[0])
    plt.plot(x, y2arr, label='Симпсон', color=colors[1])
    plt.xlabel('Кол-во делений')
    plt.ylabel('Относительная ошибка')
    plt.title('1/(1+x^2')
    plt.legend(loc='best')
    plt.savefig('page1.png')
    print('\n\nСчитаем второй интеграл:')
    print('trapeze method: \n')
    y1arr = trapeze_method(a2, b2, N, f2, Integral2)
    print('\n')
    print('Simpson method: \n')
    y2arr = simpson_method(a2, b2, N, f2, Integral2)
    plt.close()
    #plt.ylim(0, 0.1)
    plt.plot(x, y1arr, label='трапеции', color=colors[0])
    plt.plot(x, y2arr, label='Симпсон', color=colors[1])
    plt.xlabel('Кол-во делений')
    plt.ylabel('Относительная ошибка')
    plt.legend(loc='best')
    plt.title('x^(1/3)*exp(sin(x))')
    plt.savefig('page2.png')


def f1(x):
    y = 1. / (1. + x * x)
    return y

def f2(x):
    y = x ** (1./3.) * np.exp(np.sin(x))
    return y

def trapeze_method(a, b, N, f, I):
    """
    S[a,b]=(f(a) + f(b))/2 * (b - a)
    :return:
    """
    errarr = []

    for n in range(1, N):
        xarr = np.linspace(a, b, 2 ** n)
        h = xarr[1] - xarr[0]
        yarr = [f(x) for x in xarr]
        S = 0
        for y in yarr:
            S += y
        S = (S - 0.5*(yarr[0] + yarr[len(xarr)-1]))*h
        errarr.append(np.abs((I-S)/I))
        print(f' Кол-во разбиений:\t {2**n} \t\t Значение интеграла:\t{S}')

    return errarr

def simpson_method(a, b, N, f, I):
    """
    S[a,b] = (b - a)/6 * [f(a) +4f((a+b)/2) + f(b)]

    S[a,b] = h/3 [f_0_ + 4f_1_ + 2f_2_ + 4f_3_ + ... + 2f_2k-1_ + f_2k_]
    :return:
    """

    errarr = []
    for n in range(1, N):
        xarr = np.linspace(a, b, 4 ** n)
        h = xarr[1] - xarr[0]
        yarr = [f(x) for x in xarr]
        S = 0
        for i in range(len(yarr)):
            if i%2 == 0:
                if i == 0 or i == len(yarr) - 1:
                    factor = 1.
                else:
                    factor = 2.
            else:
                if i == 0 or i == len(yarr) - 1:
                    factor = 1.
                else:
                    factor = 4.
            S = S + factor * yarr[i]
        #S = (S - yarr[0] - 3. * yarr[2 ** n - 1]) * h/3.
        S = S*h/3.
        errarr.append(np.abs((I-S)/I))
        print(f' Кол-во разбиений:\t {2**n} \t\t Значение интеграла:\t{S}')

    return errarr
