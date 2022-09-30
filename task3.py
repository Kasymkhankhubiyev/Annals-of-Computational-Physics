"""
Вычислить инетегралы методами трапеций и Симпсона,
разделив отрезок интегрирования на 4, 8, 16, ... 2^n интервалов

Как убывает погрешность численного интегрирования с ростом числа интервалов?
"""
import numpy as np
a1, b1, a2, b2 = 1, -1, 0, 1
N = 10

def run():
    print('trapeze method: \n')
    trapeze_method(a1, b1, N)
    print('\n\n')
    print('Simpson method: \n')
    simpson_method(a1, b1, N)

def f1(x):
    y = 1. / (1. + x * x)
    return y

def f2(x):
    y = pow(x, 1/3) * np.exp(np.sin(x))
    return y

def trapeze_method(a, b, N):
    """
    S[a,b]=(f(a) + f(b))/2 * (b - a)
    :return:
    """
    for n in range(1, N):
        xarr = np.linspace(a, b, 2 ** n)
        h = xarr[1] - xarr[0]
        yarr = [f1(x) for x in xarr]
        S = 0
        for y in yarr:
            S += y
        S = (S - 0.5*(yarr[0] + yarr[2**n-1]))*h
        print(f' Кол-во разбиений:\t {2**n} \t\t Значение интеграла:\t{S}')

def simpson_method(a, b, N):
    """
    S[a,b] = (b - a)/6 * [f(a) +4f((a+b)/2) + f(b)]

    S[a,b] = h/3 [f_0_ + 4f_1_ + 2f_2_ + 4f_3_ + ... + 2f_2k-1_ + f_2k_]
    :return:
    """
    for n in range(1, N):
        xarr = np.linspace(a, b, 2 ** n)
        h = xarr[1] - xarr[0]
        yarr = [f1(x) for x in xarr]
        S = 0
        for i in range(len(yarr)):
            if i%2 == 0:
                factor = 2.
            else:
                factor = 4.
            S += factor * yarr[i]
        S = (S - yarr[0] - 3. * yarr[2 ** n - 1]) * h/3.
        print(f' Кол-во разбиений:\t {2**n} \t\t Значение интеграла:\t{S}')
