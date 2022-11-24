"""
Методом прогонки решить разностный аналог граничной задачи
для уравнения y''=cosx на промежутке -pi/2 < x < pi/2.

y = -cos(x) + bx + c
"""

import numpy as np
import matplotlib.pyplot as plt


def _matrix(y0: float, yn: float, x0: float, xn: float, N: int):
    d, x = [y0], [x0]

    a = np.ones(N+1)
    b = np.zeros(N+1) - 2
    c = np.ones(N+1)

    a[0], a[N] = 0, 0
    b[0], b[N] = 1, 1
    c[0], c[N] = 0, 0

    h = (xn - x0) / N  # step

    for i in range(1, N):
        xi = x0 + h * i
        d.append(h**2 * np.cos(xi))
        x.append(xi)

    d.append(yn)
    x.append(xn)
    return a, b, c, d, x


def _tridiagonal(a, b, c, d, N):
    # прямой ход метода Гаусса - исключение поддиагональных элементов ai, i =2,..,n
    for i in range(1, N + 1):
        xi = a[i] / b[i - 1]  # a_0/b_1
        b[i] -= xi * c[i - 1]
        d[i] -= xi * d[i - 1]

    y = [i for i in range(N + 1)]
    y[N] = d[N] / b[N]
    """
    $y_{n} = \frac{d_{n}}{b_{n}}; x_{i}=\frac{(d_{i} - c_{i}x_{i+1}}{b_{i}})$
    """
    for i in range(N - 1, -1, -1):
        y[i] = (d[i] - c[i] * y[i + 1]) / b[i]
    return y


def _solution(x, y0, yn, x0, xn):
    c1 = (y0 - yn + np.cos(x0) - np.cos(xn)) / (x0 - xn)
    c2 = y0 + np.cos(x0) - c1 * x0
    return -np.cos(x) + c1 * x + c2


def _draw_solution(x, y, sol, error, y0, yn):
    fig = plt.figure(figsize=(10, 4))

    plt_gr = fig.add_subplot(121)
    plt_err = fig.add_subplot(122)

    plt_gr.set_xlabel('x')
    plt_err.set_ylabel('y')
    plt_gr.plot(x, y, color='red', label='calculated')
    plt_gr.plot(x, sol, color='blue', label='analytic', linestyle='dashed')

    plt_err.set_title('error')
    plt_err.plot(x, error, color='green', label='error')

    fig.subplots_adjust(hspace=0.3, wspace=0.3)
    plt_gr.legend()
    plt_err.legend()
    plt.savefig(f'task9/task9_y0={y0}_yn={yn}.png')


def _solve_dif_equation(y0: float, yn: float, x0: float, xn: float, N: int):
    a, b, c, d, x = _matrix(y0=y0, yn=yn, x0=x0, xn=xn, N=N)
    y = _tridiagonal(a=a, b=b, c=c, d=d, N=N)
    # sol = [_solution(i) for i in x]
    sol = np.zeros_like(x)
    sol += _solution(x=np.array(x), x0=x0, xn=xn, y0=y0, yn=yn)
    error = np.abs(np.array(y) - np.array(sol))
    _draw_solution(x, y, sol, error, y0, yn)


def run():
    # y0, yn = 0, 1  # граничные условия ->
    x0, xn = -np.pi / 2, np.pi / 2  # промежуток

    N = 1000  # кол-во интервалов

    _solve_dif_equation(y0=0, yn=1, x0=x0, xn=xn, N=N)
    _solve_dif_equation(y0=1, yn=0, x0=x0, xn=xn, N=N)
    _solve_dif_equation(y0=1, yn=1, x0=x0, xn=xn, N=N)
    _solve_dif_equation(y0=0, yn=0, x0=x0, xn=xn, N=N)
    _solve_dif_equation(y0=-1, yn=1, x0=x0, xn=xn, N=N)
