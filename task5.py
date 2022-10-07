"""
Провести интерполяционный полином Pn(x) через точки
x_k = -5+k*10/n, y_k = 1/(1 + x_k^2),
k = {0, 1, ..., n} при n = 4, ..., 15.

Нарисовать графики Pn(x) - 1/(1 + x_k^2)
Объяснить результат
"""

import numpy as np
import matplotlib.pyplot as plt


def f_x(k, n):
    res = np.float64(-5) + np.float64(10.) * np.float64(k) / np.float64(n)
    return res


def f_y(x):
    res = np.float64(1.)/(np.float64(1.) + x * x)
    return res


def divided_dif(x_points, y_points, n):
    """
    f(x_0, x_1) = [f(x_1) - f(x_0)]/[x_1 - x_0]
    f(x_0, x_1, x_2)= {f(x_2) - f(x_2)]/[x_2 - x_1] - f(x_1) - f(x_0)]/[x_1 - x_0]}/[x_2 - x_0]

    """
    for i in range(1, n):
        y_points[i:] = (y_points[i:] - y_points[i-1])/(x_points[i:] - x_points[i-1])
    return y_points


def newtown(x_points, y_points, x, n):
    factors = divided_dif(x_points, y_points, n)
    pn = factors[n-1]
    for i in range(1, n):
        pn = factors[n-1-i] + (x - x_points[n-1-i]) * pn
    return pn

# def newtown(x_points, y_points, x, n):
#     factors = divided_dif(x_points, y_points, n+1)
#     pn = factors[n]
#     for i in range(1, n+1):
#         pn = factors[n-i] + (x - x_points[n-i]) * pn
#     return pn


def run():
    fig, axs = plt.subplots(nrows=4, ncols=3)
    plt.tight_layout(pad=1, h_pad=1, w_pad=1)
    for n in range(4, 16, 1):  #до 16-ти
        x_points, y_points = [], []
        for k in range(n):
            x_points.append(f_x(k, n))
            y_points.append(f_y(x_points[k]))
            #print(n, len(x_points), len(y_points))
        x = np.arange(-5., 5., 0.01)
        polynom = []
        for i in range(len(x)):
            f = f_y(x[i])
            p = newtown(x_points, y_points, x[i], n)
            polynom.append(p-f)
        axs[(n-4)//3, (n-4) % 3].plot(x, polynom, color='blue')
        axs[(n-4)//3, (n-4) % 3].set(xlabel='X')
        axs[(n-4)//3, (n-4) % 3].set(ylabel='P-f')
        axs[(n - 4) // 3, (n - 4) % 3].set(title=str(n))

    fig.set_size_inches(18.5, 10.5)
    plt.savefig('task5.png')
