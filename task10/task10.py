"""
Решить задачу Коши для одномерного уравнения
теплопроводности по схеме Кранкера-Николсона

u_t = u_xx, 0 < x < L, L = 1

u(0, t)=0, u(L, t)=0,
u(x, 0)=x*(1 - x/L)^2

как зависит от шага
влияние L на сходимость, корректность метода
"""
import matplotlib.pyplot as plt
import numpy as np
from typing import NamedTuple


class Matrix(NamedTuple):
    a: np.array
    b: np.array
    c: np.array


def ux0(x: float, L: int) -> float:
    """
    Граничное условие в начальный момент времени
    """
    return x * (1 - x/np.float64(L))**2


def matrix_coef(tau: float, h: float, N: int) -> Matrix:
    """
    Рассчет коэфициентов матрицы для З.Дирихле
    с нулевыми граничными условиями: u(0, t) = u(L, t) = 0
    a_i = 0, i = 1; a_i = -1/2 * tau/h^2, i = 2, ..., N-1
    b_i = 1 + tau/h^2
    c_i = - 1/2 * tau/h^2, i = 1, ..., N-2; c_i = 0, i = N-1
    """
    a, b, c = [], [], []

    for i in range(N-1):
        a.append(-0.5 * tau/(h*h))
        b.append(1. + tau/(h*h))
        c.append(-0.5 * tau/(h*h))

    a[0] = 0
    c[N-2] = 0

    return Matrix(a=np.array(a), b=np.array(b), c=np.array(c))


def triagonal(d: list, tau: float, h: float, N: int) -> list:

    matrix = matrix_coef(tau=tau, h=h, N=N)
    a, b, c = matrix.a, matrix.b, matrix.c

    for i in range(0, N - 1):
        ksi = a[i] / b[i - 1]
        b[i] -= ksi * c[i - 1]
        d[i] -= ksi * d[i - 1]

    y = [i for i in range(N - 1)]

    y[N - 2] = d[N - 2] / b[N - 2]

    for i in range(N - 3, -1, -1):
        y[i] = (d[i] - c[i] * y[i + 1]) / b[i]

    y.insert(0, 0)  # u(0, t) = 0
    y.append(0)  # u(L, t) = 0
    return y


def solution(t: float, x: float, L: float):
    """
    $$u(x,t) = \sum \frac{4L(2\pi k + \pi k(-1)^k)}{\pi ^4 k^4}exp(-\mu _{k}^2 t) sin(\mu _{k} x)$$
    $$\mu _{k} = \frac{k\pi}{L}$$
    """
    K = 1000
    res = []
    for k in (1, K):
        mu = (k * np.pi)/ L
        c_k = 4 * L * (2*np.pi*k + np.pi * k * np.cos(np.pi * k)) / (np.pi * k)**4
        exp = np.exp((-1) * mu**2 * t)
        sin = np.sin(mu * x)
        res.append(c_k * exp * sin)
    sol = np.sum(np.array(res))
    return sol


def run() -> None:

    N = 10
    L = 30
    T = 1
    x0, xn = 0, L
    t0, tn = 0, T

    fig, axs = plt.subplots(nrows=2, ncols=2)
    plt.tight_layout(pad=1, h_pad=1, w_pad=1)
    # plt.yscale('log')

    for n in range(1, 5):
        t_num = 10**n
        x = np.linspace(x0, xn, N+1)
        # это можно перевести в numpy если объединить в одну матрицу и обрабатывать поэлементно
        ux0j = [ux0(x=x[i], L=L) for i in range(N+1)]  # граничные условия для разных х при t=0

        v = [ux0j]
        h = np.float64(xn - x0) / np.float64(N)
        tau = np.float64(tn - t0) / np.float64(t_num)

        for t in range(t_num):
            d = [v[t][i] + tau/2. * (v[t][i+1] - 2*v[t][i] + v[t][i-1]) / (h**2) for i in range(1, N)]
            v.append(triagonal(d, tau=tau, h=h, N=N))

        # находим макс температру по всему стержню в каждый момент времени
        max_temp = np.max(np.array(v), axis=1)

        time = np.linspace(t0, tn, t_num+1)
        analytic = []
        for t in time:
            temps = []
            for step in x:
                temps.append(solution(t=t, x=step, L=L))
            analytic.append(np.max(np.array(temps)))

        error = np.abs(np.array(max_temp) / np.array(analytic) - 1)
        axs[(n-1)//2, (n-1)%2].plot(time, error, color='red', label=f'error N={N}, L={L}, t={t_num}')
        axs[(n-1)//2, (n-1)%2].legend(fontsize=7, ncol=1, facecolor='oldlace', edgecolor='r')

    # plt.plot(time, error, color='red', label=f'error N={N}, L={L}, T={T}')
    # plt.legend(fontsize=7, ncol=1, facecolor='oldlace', edgecolor='r')
    plt.savefig(f'task10/task10_L_{L}_errror.png')
    plt.close()
