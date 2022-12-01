import matplotlib.pyplot as plt
import numpy as np

"""
зависимость ошибки от шага сетки
вид ошибки
"""

N = 30
iter_number = 10

x0, xn = -15, 15
h = (xn - x0) / N

x_data = np.linspace(x0, xn, N)


def U(x):
    """
    math: $U(x) = \frac{(x^2)}{2}$
    """
    return pow(x, 2) / 2


def solution(x):
    """
    math: $\pi^{-\frac{1}{4}}exp(-\frac{x^2}{2})$
    """
    return np.exp(- x**2 / 2) * pow(np.pi, -1 / 4)


def matrix():
    a, b, c = [], [], []
    for xi in x_data:
        a.append(-0.5 / (h ** 2))
        b.append(1 / (h ** 2) + U(xi))
        c.append(-0.5 / (h ** 2))

    a[0], c[N - 1] = 0, 0

    return a, b, c


def tridiagonal(d):
    a, b, c = matrix()
    d = d.copy()

    for i in range(1, N):
        ksi = a[i] / b[i - 1]
        b[i] -= ksi * c[i - 1]
        d[i] -= ksi * d[i - 1]

    y = [i for i in range(N)]
    y[N - 1] = d[N - 1] / b[N - 1]
    for i in range(N - 2, -1, -1):
        y[i] = (d[i] - c[i] * y[i + 1]) / b[i]

    return y


def inverse_iteration(psi_prev):  # на первом шаге получаем вектор приближенных собственных значений
    for i in range(0, iter_number - 1):
        psi_prev = tridiagonal(psi_prev)
    psi = tridiagonal(psi_prev)

    energy = max(psi_prev) / max(psi)
    psi /= (max(psi) * pow(np.pi, 1 / 4))  # нормализуем волновую функцию

    return energy, psi

def error_steps():
    e_error, psi_error, nums = [], [], []
    for n in range(30, 200, 10):
        N = n
        x_data = np.linspace(x0, xn, N)
        print(n)
        nums.append(n)
        energy, psi = inverse_iteration(np.array([0.7 for _ in range(N)]))
        analytic_solution = solution(x_data)
        e_error.append(np.abs(energy-0.5))
        psi_error.append(np.sum(np.abs(np.array(analytic_solution) - np.array(psi)))/len(analytic_solution))

    fig, axs = plt.subplots(nrows=1, ncols=2)
    axs[0].plot(nums, e_error, label='energy error')
    axs[1].plot(nums, psi_error, label='psi error')
    axs[0].legend(fontsize=7, ncol=1, facecolor='oldlace', edgecolor='r')
    axs[1].legend(fontsize=7, ncol=1, facecolor='oldlace', edgecolor='r')

    plt.savefig('task11/task11_error_dependence.png')



def run():
    # energy, psi = inverse_iteration(np.array([0.7 for _ in range(N)]))
    # print(f'Calculated E0 = {energy}')
    # print(f'Real E0 = 0.5')
    #
    # fig, axs = plt.subplots(nrows=1, ncols=2)
    #
    # x = np.linspace(x0, xn, N)
    # analytic_solution = solution(x)
    #
    # axs[0].plot(x, analytic_solution, color='pink', label='analytic')
    # axs[0].plot(x, psi, color='purple', label=f'inverse iter')
    # axs[1].plot(x, np.abs(analytic_solution - np.array(psi)), color='blue', label='error')
    # axs[0].legend(fontsize=7, ncol=1, facecolor='oldlace', edgecolor='r')
    # axs[1].legend(fontsize=7, ncol=1, facecolor='oldlace', edgecolor='r')
    #
    # plt.savefig(f'task11/task11_[{x0},{xn}].png')
    error_steps()
