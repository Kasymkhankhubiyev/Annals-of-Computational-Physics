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


def difur(x: float) -> float:
    res = -1. * x
    return res


def func(t: float) -> float:
    # x = np.exp(-1 * t)
    x = np.float64(1 / np.exp(t))
    return x


def euler(M:int):
    # t0 == mint

    h = (np.float64(maxt) - np.float64(mint))/np.float64(M)

    x = [x0]
    for i in range(M):
        x.append(x[i] + h*difur(x[i]))
    return x


def runge_kutta_2order(M: int):
    """
    alpha = 3/4 - optimum
    alpha = 1/2 - corrected Euler
    alpha = 1 - modified Euler
    y_n+1 = y_n + h*[(1 - alpha)*f(x, y) + alpha*f(x + h/2*alpha, y + h/[2*alpha] * f(x_n, y_n))]
    :return: array of points
    """

    alpha = np.float64(1./2.)
    h = (np.float64(maxt) - np.float64(mint))/np.float64(M)
    x = [x0]

    for i in range(M):
        x.append(x[i] + h * ((1. - alpha)*difur(x[i]) + alpha * difur(x[i] + h * difur(x[i]) / (2 * alpha))))

    return x




def runge_kutta_4order(M:int):
    """
    y_n+1 = y_n + h/6 * [k1 + 2*k2 + 2*k3 + k4]
    k1 = f(x_n, y_n)
    k2 = f(x_n + h/2, y_n + h/2 * k1)
    k3 = f(x_n + h/2, y_n + h/2 * k3)
    k4 = f(x_n + h, y_n + h * k3)
    """

    h = (np.float64(maxt) - np.float64(mint))/np.float64(M)
    x = [x0]

    for i in range(M):
        k1 = difur(x[i])
        k2 = difur(x[i] + k1 * h / 2.)
        k3 = difur(x[i] + k2 * h / 2.)
        k4 = difur(x[i] + k3 * h)
        x.append(x[i] + h / 6. * (k1 + 2. * k2 + 2. * k3 + k4))

    return x


def one_shot_errors(num, errors_euler, errors_runge2, errors_runge4):
    plt.yscale('log')
    # plt.plot(num, errors_euler, label='euler_error', color='red')
    plt.plot(num, errors_runge2, label='runge2_error', color='blue')
    plt.plot(num, errors_runge4, label='runge4_error', color='green')

    plt.legend(fontsize=7, ncol=1, facecolor='oldlace', edgecolor='r')
    plt.savefig('task6_one_shot_errors.png')


def error_trand():
    M = 50
    errors_euler = []
    errors_runge2 = []
    errors_runge4 = []
    for K in range(10, M+10, 10):

        h = (np.float64(maxt) - np.float64(mint))/np.float64(K)

        gt = [func(mint + i * h) for i in range(K)]

        euler_sol = euler(K)
        runge2 = runge_kutta_2order(K)
        runge4 = runge_kutta_4order(K)

        euler_s = 0
        runge2_s = 0
        runge4_s = 0
        for i in range(len(gt)):
            # euler_s += (gt[i] - euler_sol[i+1])/gt[i]
            runge2_s += np.abs(gt[i] - runge2[i])/gt[i]
            runge4_s += np.abs(gt[i] - runge4[i])/gt[i]

        # errors_euler.append(euler_s/float(len(gt)))
        errors_runge2.append(runge2_s/float(len(gt)))
        errors_runge4.append(runge4_s/float(len(gt)))

    num = [K for K in range(10, M+10, 10)]

    fig, axs = plt.subplots(nrows=1, ncols=3)
    plt.tight_layout(pad=1, h_pad=1, w_pad=1)
    # axs[0].plot(num, errors_euler, label='euler_error', color='blue')
    axs[1].plot(num, errors_runge2, label='runge2_error', color='blue')
    axs[2].plot(num, errors_runge4, label='runge4_error', color='blue')

    # axs[0].legend(fontsize=7, ncol=1, facecolor='oldlace', edgecolor='r')
    axs[1].legend(fontsize=7, ncol=1, facecolor='oldlace', edgecolor='r')
    axs[2].legend(fontsize=7, ncol=1, facecolor='oldlace', edgecolor='r')

    fig.set_size_inches(10.5, 6.5)
    plt.savefig('task6_errors.png')

    fig.clear()
    plt.close(fig)

    return num, errors_euler, errors_runge2, errors_runge4



def run():
    N = 15

    h = (np.float64(maxt) - np.float64(mint)) / np.float64(N)

    gt = [func(mint + i * h) for i in range(N)]

    t = [i * h for i in range(N)]

    fig, axs = plt.subplots(nrows=2, ncols=3)
    plt.tight_layout(pad=1, h_pad=1, w_pad=1)
    # plt.yscale('log')

    # euler_sol = euler(N)
    # euler_dif = [gt[i] - euler_sol[i+1] for i in range(len(gt))]

    runge2 = runge_kutta_2order(N)
    runge2_dif = [gt[i] - runge2[i] for i in range(len(gt))]

    runge4 = runge_kutta_4order(N)
    runge4_dif = [np.abs(gt[i] - runge4[i]) for i in range(len(gt))]

    # axs[0, 0].plot(t, euler_sol[1:], label='euler', color='red')
    # axs[0, 0].plot(t, gt, label='ground_truth', color='blue')
    # axs[1, 0].plot(t, euler_dif, label='difference', color='green')

    # axs[0, 0].legend(fontsize=7, ncol=1, facecolor='oldlace', edgecolor='r')

    axs[0, 1].plot(t, runge2[1:], label='Runge 2 order', color='red')
    axs[0, 1].plot(t, gt, label='ground_truth', color='blue')

    axs[0, 1].legend(fontsize=7, ncol=1, facecolor='oldlace', edgecolor='r')

    # axs[1, 1].plot(t, euler_dif, label='euler_dif', color='green')
    axs[1, 1].plot(t, runge2_dif, label='runge2_dif', color='red')
    axs[1, 1].plot(t, runge4_dif, label='runge4_dif', color='blue')

    axs[1, 1].legend(fontsize=7, ncol=1, facecolor='oldlace', edgecolor='r')

    axs[1, 2].plot(t, runge2_dif, label='runge2_dif', color='red')
    axs[1, 2].plot(t, runge4_dif, label='runge4_dif', color='blue')

    axs[1, 2].legend(fontsize=7, ncol=1, facecolor='oldlace', edgecolor='r')

    axs[0, 2].plot(t, runge4[1:], label='Runge 4 order', color='red')
    axs[0, 2].plot(t, gt, label='ground_truth', color='blue')
    # axs[1, 2].plot(t, runge4_dif, label='difference', color='green')

    axs[0, 2].legend(fontsize=7, ncol=1, facecolor='oldlace', edgecolor='r')

    fig.set_size_inches(10.5, 6.5)
    plt.savefig('task6.png')

    fig.clear()
    plt.close(fig)

    items = error_trand()
    one_shot_errors(items[0], items[1], items[2], items[3])
