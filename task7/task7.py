"""
Решить систему диф уравнений хищник-жертва

1) x' = ax - bxy
2) y' = cxy - dy

Методом Рунге_кутты второго порядка точности
"""

import matplotlib.pyplot as plt
import numpy as np

a, b, c, d = 10, 2, 2, 10
mint, maxt = 0., 1.
alpha = 3./4.


def dif_x(x: float, y: float) -> float:
    res = float(a) * x - float(b) * x * y
    return res


def dif_y(x: float, y: float) -> float:
    res = float(c) * x * y - float(d) * y
    return res


def runge_kutta(m, x0, y0):
    """
    y_n+1 = y_n + h*((1-alpha)*f(x,y) + alpha*f(x+h/(2*alpha), y+f(x_n, y_n)*h/(2*alpha)))
    :return:
    """
    x, y = [x0], [y0]
    h = (maxt-mint)/m
    for i in range(m):
        x.append(x[i] + h * ((1 - alpha) * dif_x(x[i], y[i]) + alpha * dif_x(x[i] + h * dif_x(x[i], y[i]) / (2 * alpha),
                                                                             y[i] + h * dif_y(x[i], y[i]) / (2 * alpha))))
        y.append(y[i] + h * ((1 - alpha) * dif_y(x[i], y[i]) + alpha * dif_y(x[i] + h * dif_x(x[i], y[i]) / (2 * alpha),
                                                                             y[i] + h * dif_y(x[i], y[i]) / (2 * alpha))))
    return x, y


def run() -> None:
    N = 40
    plt.title('phase trajectory')

    colors = ['red', 'blue']

    h = (maxt - mint) / N
    t = np.arange(mint, maxt, h)

    for x0 in range(1, 6, 5):
        for y0 in range(1, 6, 5):
            results = runge_kutta(N, x0, y0)
            color = (np.random.random(), np.random.random(), np.random.random())
            plt.plot(results[0], results[1], label=str(x0) + 'vs' + str(y0), color=color)
            # plt.plot(t, results[0][1:], label='prey', color=colors[0])
            # plt.plot(t, results[1][1:], label='predator', color=colors[1])

    plt.legend(fontsize=7, ncol=1, facecolor='oldlace', edgecolor='r')
    plt.xlabel('t')
    plt.ylabel('population')
    plt.savefig(f'task7/task7_11.png')

    plt.close()

