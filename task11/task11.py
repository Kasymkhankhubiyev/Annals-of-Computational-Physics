import matplotlib.pyplot as plt
import numpy as np

"""
зависимость ошибки от шага сетки
вид ошибки
"""

# N = 30
iter_number = 10

x0, xn = -15, 15


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


def matrix(x, n):
    """
    когда n - нечетное, у нас всегда есть х==0, тогда вклад от потенциала в нуле -> 0, только от $1/h^2$
    когда n - четное, у нас нет точки в нуле и мы как бы перепрыгиваем эту точку.
    """

    h = x[1]-x[0]
    a, c = np.zeros(len(x)) - 0.5 / (h ** 2), np.zeros(len(x)) - 0.5 / (h ** 2)
    b = np.zeros(len(x)) + 1 / (h ** 2) + x**2/2
    a[0], c[n - 1] = 0, 0

    return a, b, c


def tridiagonal(d, x, n):
    a, b, c = matrix(x, n)
    d = d.copy()

    for i in range(1, n - 1):
        ksi = a[i] / b[i - 1]
        b[i] -= ksi * c[i - 1]
        d[i] -= ksi * d[i - 1]

    y = np.zeros(n)  # [i for i in range(n)]
    y[n - 1] = d[n - 1] / b[n - 1]
    for i in range(n - 2, -1, -1):
        y[i] = (d[i] - c[i] * y[i + 1]) / b[i]

    return y


def inverse_iteration(psi_prev, x, n):  # на первом шаге получаем вектор приближенных собственных значений
    for i in range(0, iter_number - 1):
        psi_prev = tridiagonal(psi_prev, x, n)
    psi = tridiagonal(psi_prev, x, n)

    energy = max(psi_prev) / max(psi)
    # $H\psi (x)= E\psi (x)$
    psi /= (max(psi) * pow(np.pi, 1 / 4))  # нормализуем волновую функцию

    return energy, psi


def error_steps():
    e_error, psi_error, nums = [], [], []
    for n in range(10, 200, 10):
        x = np.linspace(x0, xn, n)
        print(n)
        nums.append(n)
        energy, psi = inverse_iteration(np.zeros(n) + 0.7, x, n)
        analytic_solution = solution(x)
        e_error.append(np.abs(energy-0.5))
        psi_error.append(np.sum(np.abs(np.array(analytic_solution) - np.array(psi)))/len(analytic_solution))

    fig, axs = plt.subplots(nrows=1, ncols=2)
    axs[0].plot(nums, e_error, label='energy error')
    axs[1].plot(nums, psi_error, label='psi error')
    axs[0].legend(fontsize=7, ncol=1, facecolor='oldlace', edgecolor='r')
    axs[1].legend(fontsize=7, ncol=1, facecolor='oldlace', edgecolor='r')

    plt.savefig('task11/task11_error_dependence.png')



def run():
    N = 10**6  # интересно получается с четным и нечетным числом разбиений.
    x = np.linspace(x0, xn, N)
    energy, psi = inverse_iteration(np.zeros(N) + 0.7, x, N)
    print(f'Calculated E0 = {energy}')
    print(f'Real E0 = 0.5')

    fig, axs = plt.subplots(nrows=1, ncols=2)

    analytic_solution = solution(x)

    axs[0].plot(x, analytic_solution, color='pink', label='analytic')
    axs[0].plot(x, psi, color='purple', label=f'inverse iter')
    axs[1].plot(x, np.abs(analytic_solution - np.array(psi)), color='blue', label='error')
    axs[0].legend(fontsize=7, ncol=1, facecolor='oldlace', edgecolor='r')
    axs[1].legend(fontsize=7, ncol=1, facecolor='oldlace', edgecolor='r')

    plt.savefig(f'task11/task11_[{x0},{xn}].png')
    # error_steps()

    """
    Чем больше шаг сетки, тем лучше приближается
    если сетка четная - удин пик, если нечетная - два пика.
    """
