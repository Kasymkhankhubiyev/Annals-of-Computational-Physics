"""
решить жесткую систему уравнений

1) u' = 998*u + 1998*v
2) v' = -999*u - 1999*v

lambda1 = -1, lambda2 = -1000

lambda2/lambda1 = 1000 -> система жесткая

общее решение:
u(x) = a * 2*exp(-x) + b * exp(-1000*x)
v(x) = a * -exp(-x) + b * -exp(-1000*x)

пусть при x = 0 u(0) = v(0) = 1
Тогда а = 2, b = -3

"""

import matplotlib.pyplot as plt
import numpy as np
import typing

a, b, c, d = 998, 1998, -999, -1999
x0, u0, v0 = 0, 1, 1
lambda_max = 1000  # -1000

x_min, x_max = 0, 6


class Solution(typing.NamedTuple):
    name: str
    x: list
    u: list
    v: list


def dif_u(u: float, v: float) -> float:

    res = np.float(a) * u + np.float(b) * v
    return res


def dif_v(u: float, v: float) -> float:

    res = np.float(c) * u + np.float(d) * v
    return res


def linear_solution(x: list) -> Solution:

    u, v = [], []

    for i in range(len(x)):
        u.append(4.0 * np.exp(-1 * x[i]) - 3.0 * np.exp(-1000*x[i]))
        v.append(-2.0 * np.exp(-1 * x[i]) + 3.0 * np.exp(-1000*x[i]))

    return Solution(name='gt', u=u, v=v, x=x)


def explicit_first_scheme() -> Solution:
    """
    y_n+1 = y_n + h*f(x_n, y_n)
    """
    x, u, v = [x0], [u0], [v0]

    n = lambda_max * x_max
    h = np.float64(x_max) / np.float64(n)

    for i in range(n+1):
        u.append(u[i] + h * dif_u(u=u[i], v=v[i]))
        v.append(v[i] + h * dif_v(u=u[i], v=v[i]))
        x.append(x[0] + i * h)

    return Solution(name='exp1', u=u, v=v, x=x)


def implicit_first_scheme() -> Solution:
    """
    y_n+1 = y_n + h/2 * f(x_n, y_n) + h/2 * f(x_n+1, y_n+1)
    """
    n = lambda_max * x_max
    h = np.float64(x_max) / np.float64(n)

    A = np.matrix([[1 - h * a / 2, -h * b / 2],
                   [-h * c / 2, 1 - h * d / 2]])
    B = np.matrix([[1 + h * a / 2, h * b / 2],
                   [h * c / 2, 1 + h * d / 2]])
    C = np.linalg.inv(A) * B

    y = [np.array([[u0], [v0]])]

    num_steps = n
    x = [x0 + i * h for i in range(num_steps + 1)]
    for i in range(num_steps):
        y.append(C * y[i])

    # y = [(y[i][0].item(), y[i][1].item()) for i in range(len(y))]
    # print(y[40][0])

    u = [np.float64(y[i][0]) for i in range(len(y))]
    v = [np.float64(y[i][1]) for i in range(len(y))]

    return Solution(name='imp1', u=u, v=v, x=x)


def run() -> None:
    fig, axs = plt.subplots(nrows=1, ncols=2)

    solution = explicit_first_scheme()
    axs[0].plot(solution.x, solution.u, label='u(x)', color='green')
    axs[0].plot(solution.x, solution.v, label='v(x)', color='blue')
    gt = linear_solution(solution.x)
    axs[0].plot(gt.x, gt.u, label='gt u(x)', color='red', linestyle='dashed')
    axs[0].plot(gt.x, gt.v, label='gt v(x)', color='black', linestyle='dashed')
    axs[0].set(title=solution.name)
    axs[0].legend(fontsize=7, ncol=1, facecolor='oldlace', edgecolor='r')
    axs[1].plot(solution.x, np.array(gt.u)-np.array(solution.u), label='diff u(x)', color='green')
    axs[1].plot(solution.x, np.array(gt.v) - np.array(solution.v), label='diff v(x)', color='green')
    axs[1].set(title=solution.name)
    axs[1].legend(fontsize=7, ncol=1, facecolor='oldlace', edgecolor='r')
    plt.savefig('task8/task8_exp1.png')
    plt.close()

    fig, axs = plt.subplots(nrows=1, ncols=2)
    solution = implicit_first_scheme()
    axs[0].plot(solution.x, solution.u, label='u(x)', color='green')
    axs[0].plot(solution.x, solution.v, label='v(x)', color='blue')
    gt = linear_solution(solution.x)
    axs[0].plot(gt.x, gt.u, label='gt u(x)', color='red', linestyle='dashed')
    axs[0].plot(gt.x, gt.v, label='gt v(x)', color='black', linestyle='dashed')
    axs[0].set(title=solution.name)
    axs[0].legend(fontsize=7, ncol=1, facecolor='oldlace', edgecolor='r')
    axs[1].plot(solution.x, np.array(gt.u) - np.array(solution.u), label='diff u(x)', color='green')
    axs[1].plot(solution.x, np.array(gt.v) - np.array(solution.v), label='diff v(x)', color='green')
    axs[1].set(title=solution.name)
    axs[1].legend(fontsize=7, ncol=1, facecolor='oldlace', edgecolor='r')
    plt.savefig('task8/task8_imp1.png')
    plt.close()

