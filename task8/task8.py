"""
решить жесткую систему уравнений

1) u' = 998*u + 1998*v
2) v' = -999*u - 1999*v

lambda1 = -1, lambda2 = -1000

lambda2/lambda1 = 1000 -> система жесткая

общее решение:
u(x) = a * 2*exp(-x) + b * exp(-1000*x)
v(x) = a * -exp(-x) + b * -exp(-1000*x)
"""

import matplotlib.pyplot as plt
import numpy as np
import typing

a, b, c, d = 998, 1998, -999, -1999
x0, u0, v0 = 0, 1, 1
lambda_max = 1000

x_min, x_max = 0, 1


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
    a = u0 + v0
    b = -(u0 + 2*v0)
    for i in range(len(x)):
        u.append(2.0 * a * np.exp(-1 * x[i]) + b * np.exp(-1000*x[i]))
        v.append(-1.0 * a * np.exp(-1 * x[i]) - b * np.exp(-1000*x[i]))
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

    u = [np.float64(y[i][0]) for i in range(len(y))]
    v = [np.float64(y[i][1]) for i in range(len(y))]

    return Solution(name='imp1', u=u, v=v, x=x)


def newtown():
    """
    f = 998u + 1998v
    g = -999u - 1999v
    """
    f_u, f_v, g_u, g_v = 998, 1998, -999, -1999
    n = lambda_max * x_max
    h = np.float64(x_max) / np.float64(n)
    x = np.arange(x_min, x_max, h)
    E = np.eye(2)
    J = np.matrix([[f_u, f_v], [g_u, g_v]])   # constant
    A = np.linalg.inv(E - h*J)

    u, v = [u0], [v0]

    for i in range(len(x)):
        res = np.matrix([[u[i]], [v[i]]]) + A * np.matrix([[h * dif_u(u=u[i], v=v[i])], [h * dif_v(u=u[i], v=v[i])]])
        u.append(res[0, 0])
        v.append(res[1, 0])
    return Solution(name='Newtown', x=x, u=u, v=v)


def run() -> None:
    fig, axs = plt.subplots(nrows=1, ncols=3)

    solution_exp = explicit_first_scheme()
    axs[0].plot(solution_exp.x, solution_exp.u, label='u(x)', color='green')
    axs[0].plot(solution_exp.x, solution_exp.v, label='v(x)', color='blue')
    gt = linear_solution(solution_exp.x)
    axs[0].plot(gt.x, gt.u, label='gt u(x)', color='red', linestyle='dashed')
    axs[0].plot(gt.x, gt.v, label='gt v(x)', color='black', linestyle='dashed')
    axs[0].set(title=solution_exp.name)
    axs[0].legend(fontsize=7, ncol=1, facecolor='oldlace', edgecolor='r')

    solution_imp = implicit_first_scheme()
    axs[1].plot(solution_imp.x, solution_imp.u, label='u(x)', color='green')
    axs[1].plot(solution_imp.x, solution_imp.v, label='v(x)', color='blue')
    gt = linear_solution(solution_imp.x)
    axs[1].plot(gt.x, gt.u, label='gt u(x)', color='red', linestyle='dashed')
    axs[1].plot(gt.x, gt.v, label='gt v(x)', color='black', linestyle='dashed')
    axs[1].set(title=solution_imp.name)
    axs[1].legend(fontsize=7, ncol=1, facecolor='oldlace', edgecolor='r')

    solution_newtown = newtown()
    axs[2].plot(solution_newtown.x, solution_newtown.u[1:], label='u(x)', color='green')
    axs[2].plot(solution_newtown.x, solution_newtown.v[1:], label='v(x)', color='blue')
    gt = linear_solution(solution_newtown.x)
    axs[2].plot(gt.x, gt.u, label='gt u(x)', color='red', linestyle='dashed')
    axs[2].plot(gt.x, gt.v, label='gt v(x)', color='black', linestyle='dashed')
    # plt.set(title=solution_newtown.name)
    axs[2].legend(fontsize=7, ncol=1, facecolor='oldlace', edgecolor='r')
    # plt.savefig('task8/task8_newtown.png')

    plt.savefig('task8/task8_exp1.png')
    plt.close()

    # results = newtown()
    # print(results)

