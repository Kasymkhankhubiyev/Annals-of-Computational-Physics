"""
Методом погонки решить разностный аналог граничной задачи
для уравнения y''=cosx на промежутке -pi/2 < x < pi/2.

Рассмотреть различный варианты граничных условий (га функцию и ее производную).
"""

import numpy as np
import matplotlib.pyplot as plt

Constant = float

# задаем граничные условия
y0, yn = 0, 1
x0, xn = -np.pi/2., np.pi/2.

N = 100
h = np.float64(xn - x0)/np.float64(N)

a0, b0, c0 = 1, 1, 1


def second_diff(x: float) -> float:
    res = np.cos(x)
    return res


def first_diff(x: float, c0: Constant):
    """
    y' = sinx + c0
    """
    res = np.sin(x) + c0
    return res


def func(x: float, c0: Constant, c1: Constant):
    """
    y = -cosx + c0*x + c1
    """
    res = -1 * np.cos(x) + c0 * x + c1
    return res


def calculation_first_order_const(y_0, y_n):
    """
    правая прогонка

    граничные условия первого рода, константные

    Fi = -cos(x)
    yi'' = y_i+1 /h^2 - y_i * 2 / h^2 + y_i-1/h^2
    Ai = 1/h^2, Bi = 1/h^2, Ci = 2/h^2
    n1 = y0, n2 = yn
    alfa1 = x0, betta1 = y0
    """
    h = np.float64(xn - x0)/np.float64(N)

    A, B, C = [], [], [2*np.float64(1)/np.float64(h*h)]
    F = [-np.cos(x0)]

    alfas, bettas = [1], [1]
    y = [y_n]

    for i in range(1, N):
        A.append(np.float64(1)/np.float64(h*h))
        C.append(C[0])
        B.append(np.float64(1)/np.float64(h*h))
        F.append(- np.cos(x0 + i * h))

    for i in range(0, N-2):
        alfas.append(B/(C-A*alfas[i]))
        bettas.append((A*B - np.cos(x0 + (i+1) * h))/(C - A * alfas[i]))

    for i in range(0, N-1):
        y.append(alfas[i] * y[i] + bettas[i])

    return y



def run() -> None:

    fig, axs = plt.subplots(nrows=1, ncols=2)

    c0 = np.float64(yn + np.cos(xn) - y0 - np.cos(x0)) / np.float64(xn - x0)
    c1 = y0 + np.cos(x0) - c0 * x0

    h = np.float64(xn - x0)/np.float64(N)

    x, gt = [], []

    for i in range(N):
        x.append(x0 + i*h)
        gt.append(func(x=x0+i*h, c0=c0, c1=c1))

    axs[0].plot(x, gt, label='analytic', color='green')
    axs[0].legend(fontsize=7, ncol=1, facecolor='oldlace', edgecolor='r')

    calc = calculation_first_order_const(y_0=y0, y_n=yn)
    axs[1].plot(x, calc, color='red', label='calculated')
    axs[1].legend(fontsize=7, ncol=1, facecolor='oldlace', edgecolor='r')

    plt.savefig('task9/task9.png')
    plt.close()



