import matplotlib.pyplot as plt
import numpy as np

N = 2
iter_number = 2

x0, xn = -10, 10
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
    return pow(np.pi, -1 / 4) * np.exp(- pow(x, 2) / 2)


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

    for i in range(0, N):
        ksi = a[i] / b[i - 1]
        b[i] -= ksi * c[i - 1]
        d[i] -= ksi * d[i - 1]

    y = np.array(0, N, 1, dtype=int)  # [i for i in range(N)]
    y[N - 1] = d[N - 1] / b[N - 1]

    for i in range(N - 2, -1, -1):
        y[i] = (d[i] - c[i] * y[i + 1]) / b[i]

    return y


def inverse_iteration(psi_prev):
    for i in range(0, iter_number - 1):
        psi_prev = tridiagonal(psi_prev)
    psi = tridiagonal(psi_prev)

    energy = max(psi_prev) / max(psi)
    psi /= (max(psi) * pow(np.pi, 1 / 4))  # нормализуем волновую функцию

    return energy, psi


def run():
    energy, psi = inverse_iteration(np.array([0.7 for _ in range(N)]))
    print(f'Calculated E0 = {energy}')
    print(f'Real E0 = 0.5')

    fig = plt.figure(figsize=(12, 4))
    plt_wave = fig.add_subplot(121)
    plt_err = fig.add_subplot(122)

    plt_wave.set_title(r'Wave-function')
    plt_wave.plot(x_data, solution(x_data), color='pink', label='real solution')
    plt_wave.plot(x_data, psi, color='purple', label=f'inverse iteration,\nN = {N},\niteration = {iter_number}')

    plt_err.set_title('Error value')
    plt_err.plot(x_data, solution(x_data) - np.array(psi), color='blue')

    plt_wave.legend()
    plt.savefig('task11/task11.png')