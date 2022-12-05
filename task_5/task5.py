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
    for i in range(1, n+1):
        y_points[i:] = (y_points[i:] - y_points[i-1])/(x_points[i:] - x_points[i-1])
    return y_points


def newtown(x_points, y_points, x, n):
    factors = divided_dif(x_points, y_points, n)
    pn = factors[n-1]
    for i in range(1, n-1):
        pn = factors[n-i] + (x - x_points[n-i]) * pn
    return pn


def div_dif(x_points, y_points, n):

    """
    y(x_1, x_2, ..., x_n) = SUM(i=1; n){y(x_i)/П(x_i-x_j)}
    """
    div_difs = [y_points[0]]
    for i in range(1, n):
        summ = 0.
        for j in range(i):
            frac = 1.
            for k in range(i):
                if k != j:
                    frac *= (x_points[j] - x_points[k])
            summ += y_points[j] / frac
        div_difs.append(summ)

    return div_difs


def newtown_ext(x_points, y_points, x, n):
    factors = div_dif(x_points, y_points, n)
    pn = factors[n-1]
    for i in range(n):
        pn = factors[n-1-i] + (x - x_points[n-1-i]) * pn
    return pn


def lagrange_factors(x_points: list, n: int) -> list:
    factors = []
    for j in range(n):
        factor = 1.
        for i in range(n):
            if j != i:
                factor *= (x_points[j] - x_points[i])
        factors.append(factor)

    return factors


def Lagrange(x_points, y_points, x, n) -> float:

    factors = lagrange_factors(x_points, n)
    polynom = 0
    for i in range(n):
        phi = 1.

        for j in range(n):
            phi *= (x - x_points[i])

        polynom += phi * y_points[i] / factors[i]
    return polynom


def run():
    fig, axs = plt.subplots(nrows=5, ncols=3)
    plt.tight_layout(pad=1, h_pad=1, w_pad=1)
    for n in range(4, 16, 1):  #до 16-ти
        x_points, y_points = [], []
        for k in range(n+1):
            x_points.append(f_x(k, n))
            y_points.append(f_y(x_points[k]))
            #print(n, len(x_points), len(y_points))
        x = np.arange(-5., 5., 0.01)
        polynom = []
        f_s = []
        for i in range(len(x)):
            f = f_y(x[i])
            f_s.append(f)
            p = newtown_ext(x_points, y_points, x[i], n+1)
            polynom.append(p)
        axs[(n - 4) // 3, (n - 4) % 3].plot(x, polynom, color='blue')
        axs[(n - 4) // 3, (n - 4) % 3].plot(x, f_s, color='blue')
        axs[(n - 4) // 3, (n - 4) % 3].set(xlabel='X')
        axs[(n - 4) // 3, (n - 4) % 3].set(ylabel='f')
        axs[(n - 4) // 3, (n - 4) % 3].set(title=str(n))
        # print(f_s)

    fig.set_size_inches(18.5, 10.5)
    plt.savefig('task_5/task5_1.png')

def run_lagrange():
    fig, axs = plt.subplots(nrows=4, ncols=3)
    plt.tight_layout(pad=1, h_pad=1, w_pad=1)
    for n in range(4, 16, 1):  #до 16-ти
        x_points, y_points = [], []
        for k in range(n+1):
            x_points.append(f_x(k, n))
            y_points.append(f_y(x_points[k]))
            #print(n, len(x_points), len(y_points))
        x = np.arange(-5., 5., 0.01)
        polynom = []
        f_points, p_points = [], []
        for i in range(len(x)):
            f = f_y(x[i])

            p = Lagrange(x_points, y_points, x[i], n+1)

            polynom.append(p - f)
        axs[(n-4)//3, (n-4) % 3].plot(x, polynom, color='blue')
        axs[(n-4)//3, (n-4) % 3].set(xlabel='X')
        axs[(n-4)//3, (n-4) % 3].set(ylabel='P-f')
        axs[(n - 4) // 3, (n - 4) % 3].set(title=str(n))

    fig.set_size_inches(18.5, 10.5)
    plt.savefig('task_5/task5_lagrange.png')


def run_lagrange_two_plots():
    fig, axs = plt.subplots(nrows=4, ncols=3)
    plt.tight_layout(pad=1, h_pad=1, w_pad=1)
    for n in range(4, 16, 1):  # до 16-ти
        x_points, y_points = [], []
        for k in range(n):
            x_points.append(f_x(k, n))
            y_points.append(f_y(x_points[k]))
            # print(n, len(x_points), len(y_points))
        x = np.arange(-5., 5., 0.01)
        polynom = []
        f_points, p_points = [], []
        for i in range(len(x)):
            f = f_y(x[i])
            f_points.append(f)

            p = Lagrange(x_points, y_points, x[i], n+1)
            p_points.append(p)

            # print(f, p)

        axs[(n - 4) // 3, (n - 4) % 3].plot(x, f_points, color='blue')
        axs[(n - 4) // 3, (n - 4) % 3].plot(x, p_points, color='red')
        axs[(n - 4) // 3, (n - 4) % 3].set(xlabel='X')
        axs[(n - 4) // 3, (n - 4) % 3].set(ylabel='P-f')
        axs[(n - 4) // 3, (n - 4) % 3].set(title=str(n))

    fig.set_size_inches(18.5, 10.5)
    plt.savefig('task_5/task5_two_func.png')
