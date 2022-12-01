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
    b[0], c[0] = -1, 1
    c[N-1], b[N-1] = -1, 0
    return Matrix(a=a, b=b, c=c, d=d)


def inverse_iteration(N, tau, h, t_steps, L):
    matrix = _matrix(N, tau, h)
    ux, uy = np.zeros((N, t_steps)), np.zeros((N, t_steps))
    x, y = np.linspace(-L, L, N), np.linspace(-L, L, N)
    for i in range(N):
        ux[i][0], uy[i][0] = uxy0(x=x[i], y=y[i], L=L)  #заполняем граничные условия при t=0

    t_origin = np.zeros(t_steps)
    t_origin[0] = ux[N / 2][0] * uy[N / 2][0]

    for t_i in range(t_steps-1):

        if t_i % 2 == 0:  # делаем шаг по х
            for j in range(1, N-1):
                uy[j][t_i + 1] = uy[j][t_i] + tau/2 * (uy[j+1][t_i] - 2*uy[j][t_i] + uy[j-1][t_i]) / h**2 # шаг по неявной схеме
                matrix.d[j] = uy[j][t_i + 1]

            for i in range(2, N):  # прямой ход
                ksi = matrix.a[i] / matrix.b[i - 1]
                matrix.a[i] = 0
                matrix.b[i] = matrix.b[i] - ksi * matrix.c[i - 1]
                matrix.d[i] = matrix.d[i] - ksi * matrix.d[i - 1]
            ux[N - 2][t_i + 1] = matrix.d[N - 2] / matrix.b[N - 2]
            for i in range(N - 2, 0, -1):
                ux[i][t_i + 1] = (1 / matrix.b[i] * (matrix.d[i] - matrix.c[i] * ux[i + 1][t_i + 1]))


def solve(N: int, L: int, tau: float, t: int, t_steps: int, x: np.array, y:np.array):
    uxy0j = np.zeros((N, N))  # граничные условия для разных х, y при t=0
    for i in range(N):
        for j in range(N):
            uxy0j[i][j] = uxy0(x[i], y[j])

    v = [uxy0j]
    h = x[1] - x[0]  # np.float64(xn - x0) / np.float64(N)

    for t in range(t_steps):
        # d = [v[t][i] + tau / 2. * (v[t][i + 1] - 2 * v[t][i] + v[t][i - 1]) / (h ** 2) for i in range(1, len(x) - 1)]

        v.append(inverse_iteration(d, tau=tau, h=h, N=len(x) - 1))
    return v




def run():
    L, N, t = 1, 100, 10
    t_steps = 1000
    # tau_step = 1e-2
    x0, y0, xn, yn = -L, -L, L, L
    x, y, t = np.linspace(-L, L, N), np.linspace(-L, L, N), np.linspace(0, t, t_steps)
    ux, uy = np.zeros((N, t_steps)), np.zeros((N, t_steps))  # создаем матрицу N - кол-во узлов пространственной сетки
    u = np.zeros((N, N, t_steps))


