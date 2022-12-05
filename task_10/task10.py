"""
Решить задачу Коши для одномерного уравнения
теплопроводности по схеме Кранкера-Николсона

u_t = u_xx, 0 < x < L, L = 1

u(0, t)=0, u(L, t)=0,
u(x, 0)=x*(1 - x/L)^2

как ывыставлять шаг в зависимости от L?
видим, что достаточно, чтобы шах по Х и по времени был на порядок меньше длины.
"""
import matplotlib.pyplot as plt
import numpy as np
from typing import NamedTuple


class Matrix(NamedTuple):
    a: np.array
    b: np.array
    c: np.array


def ux0(x: float, L: float) -> float:
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


def solve(N, x0, xn, tn, t0, L):
    t_num = 10000  # round(N/10)
    x = np.linspace(x0, xn, N + 1)
    x = x[::int(np.log10(N))+1]
    # это можно перевести в numpy если объединить в одну матрицу и обрабатывать поэлементно
    ux0j = [ux0(x=x[i], L=L) for i in range(len(x))]  # граничные условия для разных х при t=0

    v = [ux0j]
    h = x[1] - x[0]#np.float64(xn - x0) / np.float64(N)
    tau = np.float64(tn - t0) / np.float64(t_num)

    for t in range(t_num):
        d = [v[t][i] + tau / 2. * (v[t][i + 1] - 2 * v[t][i] + v[t][i - 1]) / (h ** 2) for i in range(1, len(x)-1)]
        v.append(triagonal(d, tau=tau, h=h, N=len(x)-1))
    return v

def run() -> None:

    N = 10
    L = 1e-3
    T = 1
    x0, xn = 0, L
    t0, tn = 0, 0.1

    fig, axs = plt.subplots(nrows=3, ncols=2)
    plt.tight_layout(pad=1, h_pad=1, w_pad=1)
    plt.yscale('log')

    time_steps = 10001

    time = np.linspace(t0, tn, time_steps)

    # находим макс температру по всему стержню в каждый момент времени
    print('step1')
    h = 1e-4
    N = round(L/h)
    max_temp = np.max(np.array(solve(N=N, x0=x0, xn=xn, tn=tn, t0=t0, L=L)), axis=1)
    axs[0, 0].plot(time, max_temp, color='red', label=f'error N={N}, L={L}, pow={L/(xn-x0) * N}')
    axs[0, 0].legend(fontsize=7, ncol=1, facecolor='oldlace', edgecolor='r')

    # находим макс температру по всему стержню в каждый момент времени
    print('step2')
    h = 1e-4
    N = round(L / h)
    max_temp = np.max(np.array(solve(N=N, x0=x0, xn=xn, tn=tn, t0=t0, L=L)), axis=1)
    axs[0, 1].plot(time, max_temp, color='red', label=f'error N={N}, L={L}, pow={L / (xn - x0) * N}')
    axs[0, 1].legend(fontsize=7, ncol=1, facecolor='oldlace', edgecolor='r')

    # находим макс температру по всему стержню в каждый момент времени
    print('step3')
    h = 1e-5
    N = round(L / h)
    max_temp = np.max(np.array(solve(N=N, x0=x0, xn=xn, tn=tn, t0=t0, L=L)), axis=1)
    axs[1, 0].plot(time, max_temp, color='red', label=f'error N={N}, L={L}, pow={L / (xn - x0) * N}')
    axs[1, 0].legend(fontsize=7, ncol=1, facecolor='oldlace', edgecolor='r')

    # находим макс температру по всему стержню в каждый момент времени
    print('step4')
    h = 1e-6
    N = round(L / h)
    max_temp = np.max(np.array(solve(N=N, x0=x0, xn=xn, tn=tn, t0=t0, L=L)), axis=1)
    axs[1, 1].plot(time, max_temp, color='red', label=f'error N={N}, L={L}, pow={L / (xn - x0) * N}')
    axs[1, 1].legend(fontsize=7, ncol=1, facecolor='oldlace', edgecolor='r')

    # находим макс температру по всему стержню в каждый момент времени
    print('step5')
    h = 1e-7
    N = round(L / h)
    max_temp = np.max(np.array(solve(N=N, x0=x0, xn=xn, tn=tn, t0=t0, L=L)), axis=1)
    axs[2, 0].plot(time, max_temp, color='red', label=f'error step_x={h}, L={L}, time_step={tn/(time_steps-1)}')
    axs[2, 0].legend(fontsize=7, ncol=1, facecolor='oldlace', edgecolor='r')

    # находим макс температру по всему стержню в каждый момент времени
    # print('step6')
    # h = 1e-8
    # N = round(L / h)
    # max_temp = np.max(np.array(solve(N=N, x0=x0, xn=xn, tn=tn, t0=t0, L=L)), axis=1)
    # axs[2, 1].plot(time, max_temp, color='red', label=f'error N={N}, L={L}, pow={L / (xn - x0) * N}')
    # axs[2, 1].legend(fontsize=7, ncol=1, facecolor='oldlace', edgecolor='r')

    plt.savefig(f'task10/task10_temps_L_{L}_l_dependence_time_steps={time_steps}.png')
    plt.close()

    #там где L=30 слишком маленький интервал времени, поэтому линейно спадает
