import numpy as np
import matplotlib.pyplot as plt

from typing import NamedTuple


class Matrix(NamedTuple):
    a: np.array
    b: np.array
    c: np.array
    d: np.array


def uxy0(x, y, L):
    return (1 - (x/L)**2), (1 - (y/L)**2)


def _matrix(N, tau, h):
    a, b, c, d = np.zeros(N) - tau/h**2, np.ones(N) + 2*tau/h**2, np.zeros(N) - tau/h**2, np.zeros(N)
    a[0], b[0], c[0] = 0, -1, 1
    a[N-1], b[N-1], c[N-1] = 1, -1, 0
    return Matrix(a=a, b=b, c=c, d=d)


def _step_along_y(matrix, N, t_i, tau, ux, uy, h):
    for j in range(1, N - 1):
        uy[j][t_i + 1] = uy[j][t_i] + tau / 2 * (
                    uy[j + 1][t_i] - 2 * uy[j][t_i] + uy[j - 1][t_i]) / h ** 2  # шаг по неявной схеме
        matrix.d[j] = uy[j][t_i + 1]

    for i in range(2, N):  # прямой ход
        ksi = matrix.a[i] / matrix.b[i - 1]
        matrix.a[i] = 0
        matrix.b[i] = matrix.b[i] - ksi * matrix.c[i - 1]
        matrix.d[i] = matrix.d[i] - ksi * matrix.d[i - 1]

    ux[N - 2][t_i + 1] = matrix.d[N - 2] / matrix.b[N - 2]  # обратный ход по x

    for i in range(N - 2, 0, -1):  # обратный ход по x
        ux[i][t_i + 1] = (1 / matrix.b[i] * (matrix.d[i] - matrix.c[i] * ux[i + 1][t_i + 1]))

    return ux, uy


def _step_along_x(matrix, N, t_i, tau, ux, uy, h):
    for j in range(1, N - 1):  # неявный шаг по х
        ux[j][t_i + 1] = ux[j][t_i] + tau / 2 * (ux[j + 1][t_i] - 2 * ux[j][t_i] + ux[j - 1][t_i]) / h ** 2
        matrix.d[j] = ux[j][t_i + 1]

    for i in range(2, N):  # прямой ход по х
        ksi = matrix.a[i] / matrix.b[i - 1]
        matrix.a[i] = 0
        matrix.b[i] = matrix.b[i] - ksi * matrix.c[i - 1]
        matrix.d[i] = matrix.d[i] - ksi * matrix.d[i - 1]

    uy[N - 2][t_i + 1] = matrix.d[N - 2] / matrix.b[N - 2]  # обратный ход по у

    for i in range(N - 2, 0, -1):  # обратный ход по у
        uy[i][t_i + 1] = (1 / matrix.b[i] * (matrix.d[i] - matrix.c[i] * uy[i + 1][t_i + 1]))

    return ux, uy


def localy_1d_method(N, tau, h, t_steps, L):
    # matrix = _matrix(N, tau, h)
    ux, uy = np.zeros((N, t_steps)), np.zeros((N, t_steps))
    x, y = np.linspace(-L, L, N), np.linspace(-L, L, N)
    for i in range(N):
        ux[i][0], uy[i][0] = uxy0(x=x[i], y=y[i], L=L)  #заполняем граничные условия при t=0

    temp_origin = np.zeros(t_steps)
    temp_origin[0] = ux[int(N / 2)][0] * uy[int(N / 2)][0]

    for t_i in range(t_steps-1):  # шаг по временной сетке
        matrix = _matrix(N, tau, h)

        if t_i % 2 == 0:  # делаем шаг по y
            ux, uy = _step_along_y(matrix=matrix, N=N, t_i=t_i, tau=tau, ux=ux, uy=uy, h=h)

        if t_i % 2 == 1:  # делаем шаг по x
            ux, uy = _step_along_x(matrix=matrix, N=N, t_i=t_i, tau=tau, ux=ux, uy=uy, h=h)

        temp_origin[t_i + 1] = ux[int(N / 2), t_i + 1] * uy[int(N / 2), t_i + 1]

    return temp_origin


def run():
    L, N, t = 1, 200, 1
    t_steps = 101
    x, t = np.linspace(-L, L, N), np.linspace(0, t, t_steps)
    temps = localy_1d_method(N=N, tau=t[1] - t[0], h=x[1] - x[0], t_steps=t_steps, L=L)
    plt.plot(t, temps, color='blue', label='temp in the origin')
    plt.plot(t, np.exp(-6*t), color='red', label='exp(-6x)')
    plt.plot(t, np.exp(-7*t), color='green', label='exp(-7x)')
    plt.legend(fontsize=7, ncol=1, facecolor='oldlace', edgecolor='r')
    plt.savefig('task13/temp_vs_time.png')
