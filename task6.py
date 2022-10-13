"""
Решить задачу Коши:

    dx/dt = -x, x(0)=1, 0<t<3

методом Эйлера первого порядка точности и
методами Рунге-Кутты второго и четвертого порядка точности.

 dx/dt = -x <-> dx/x = -t <-> x = exp(-t)

"""
import numpy as np
import matplotlib.pyplot as plt

mint, maxt = 0, 3
x0, t0 = 1, 0
N = 10


def difur(x: float) -> float:
    res = -x
    return res


def func(t: float) -> float:
    x = np.float64(1 / np.exp(t))
    return x


def euler(M:int):
    # t0 == mint

    h = (np.float64(maxt) - np.float64(mint))/np.float64(M)

    x = [x0]
    for i in range(M):
        x.append(x[i] + h*difur(x[i]))
    return x


def runge_kutta_2order():
    pass


def runge_kutta_4order():
    pass


def error_trand():
    M = 500
    errors = []
    for K in range(10, M+10, 10):
        t = np.linspace(mint, maxt, M)

        gt = [func(i) for i in t]

        euler_sol = euler(K)

        s = 0
        for i in range(len(gt)):
            s += (gt[i] - euler_sol[i+1])/gt[i]

        errors.append(s/K)

    num = [K for K in range(10, M+10, 10)]

    fig, axs = plt.subplots(nrows=1, ncols=3)
    plt.tight_layout(pad=1, h_pad=1, w_pad=1)
    axs[0].plot(num, errors, label='euler_error', color='blue')

    axs[0].legend(fontsize=7,
                  ncol=1,
                  facecolor='oldlace',
                  edgecolor='r')

    fig.set_size_inches(10.5, 6.5)
    plt.savefig('task6_errors.png')


def run():

    t = np.linspace(mint, maxt, N)

    gt = [func(i) for i in t]

    fig, axs = plt.subplots(nrows=1, ncols=3)
    plt.tight_layout(pad=1, h_pad=1, w_pad=1)
    euler_sol = euler(N)
    euler_dif = [gt[i] - euler_sol[i+1] for i in range(len(gt))]

    axs[0].plot(t, euler_sol[1:], label='euler', color='red')
    axs[0].plot(t, gt, label='ground_truth', color='blue')
    axs[0].plot(t, euler_dif, label='difference', color='green')

    axs[0].legend(fontsize=7,
                  ncol=1,
                  facecolor='oldlace',
                  edgecolor='r')

    fig.set_size_inches(10.5, 6.5)
    plt.savefig('task6.png')
