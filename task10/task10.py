"""
Решить задачу Коши для одномерного уравнения
теплопроводности по схеме Кранкера-Николсона

u_t = u_xx, 0 < x < L, L = 1

u(0, t)=0, u(L, t)=0,
u(x, 0)=x*(1 - x/L)^2
"""

import matplotlib.pyplot as plt
import numpy as np
from typing import NamedTuple

N = 10
L = 1
T = 10
x0, xn = 0, L
t0, tn = 0, 1

h = np.float64(xn - x0)/np.float64(N)
tau = np.float64(tn - t0)/np.float64(T)


class Matrix(NamedTuple):
    a: np.array
    b: np.array
    c: np.array

def ux0(x: float) -> float:
    return x * (1 - x/np.float64(L))**2


def matrix_coef() -> Matrix:
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


def triagonal(d):

    matrix = matrix_coef()
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


def run() -> None:
    x = [x0 + i * h for i in range(N+1)]
    ux0j = [ux0(x[i]) for i in range(N+1)]  # граничные условия для разных х при t=0

    v = [ux0j]

    """
    f(x, t) = 0
    """
    for t in range(T):
        d = [v[t][i] + tau/2. * (v[t][i+1] - 2*v[t][i] + v[t][i-1]) / (h**2) for i in range(1, N)]
        v.append(triagonal(d))

    max_temp = [max(v[i]) for i in range(len(v))]
    time = [t0 + tau * i for i in range(len(v))]

    plt.plot(time, max_temp, color='green', label='max temp')
    plt.title('Зависимость максимальной температуры от времени')
    plt.xlabel('t')
    plt.ylabel('temp')
    plt.savefig('task10/task10.png')