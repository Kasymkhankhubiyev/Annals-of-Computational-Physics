import numpy as np
import matplotlib.pyplot as plt

min_x, max_x = 0, 6
x0, u0, v0 = 0, 1, 1

# h < 2/lambda // из решения по явной схеме эйлера, условие для устойчивости системы
lamb_max = -1000
N_ex = abs(lamb_max) * max_x
h_beg, h_end = (max_x - min_x) / N_ex, (max_x - min_x) / N_ex

N_im = 20
h_beg_im, h_end_im = (max_x - min_x) / N_im, (max_x - min_x) / N_im

# две h нужны: тк в самом начале слишком малые изменения по x,
# поэтому необходим более маленький h, далее уже можно взять покрупнее

a, b, c, d = 998, 1998, -999, -1999


def solution(x):
    return 4.0 * np.exp(-x) - 3.0 * np.exp(-1000 * x), -2.0 * np.exp(-x) + 3.0 * np.exp(-1000 * x)


def diff_u(u, v):
    return a * u + b * v


def diff_v(u, v):
    return c * u + d * v


def det(h):
    return (1 - h * d) * (1 - h * a) - b * c * pow(h, 2)


def diff_imp_u(u, v, h):
    return ((1 - h * d) * u + b * h * v) / det(h)


def diff_imp_v(u, v, h):
    return ((1 - h * a) * v + c * h * u) / det(h)


def euler():
    x, u, v = [x0], [u0], [v0]
    for i in range(0, 101):
        u.append(u[i] + h_beg * diff_u(u[i], v[i]))
        v.append(v[i] + h_beg * diff_v(u[i], v[i]))
        x.append(x[0] + i * h_beg)
    for i in range(101, N_ex + 1):
        u.append(u[i] + h_end * diff_u(u[i], v[i]))
        v.append(v[i] + h_end * diff_v(u[i], v[i]))
        x.append(x[0] + i * h_end)
    return u, v, x


def euler_im():
    x, u, v = [x0], [u0], [v0]
    for i in range(0, 101):
        u.append(diff_imp_u(u[i], v[i], h_beg_im))
        v.append(diff_imp_v(u[i], v[i], h_beg_im))
        x.append(x[0] + i * h_beg_im)
    for i in range(101, N_im + 1):
        u.append(diff_imp_u(u[i], v[i], h_end_im))
        v.append(diff_imp_v(u[i], v[i], h_end_im))
        x.append(x[0] + i * h_end_im)
    return u, v, x


def run():
    x_data = np.arange(min_x, max_x, h_beg)
    y = solution(x_data)
    u_solution, v_solution = y[0], y[1]

    u_ex, v_ex, x_ex = euler()
    u_im, v_im, x_im = euler_im()

    fig = plt.figure(figsize=(10, 4))

    plt_ex = fig.add_subplot(121)
    plt_im = fig.add_subplot(122)

    plt_ex.set_title('explicit euler')
    plt_ex.set_xlabel('x')
    plt_ex.plot(x_ex, u_ex, color='purple', label='u(x)')
    plt_ex.plot(x_ex, v_ex, color='blue', label='v(x)')
    plt_ex.plot(x_data, u_solution, color='purple', label='u(x) sol', linestyle='dashed')
    plt_ex.plot(x_data, v_solution, color='blue', label='v(x) sol', linestyle='dashed')

    plt_im.set_title('implicit euler')
    plt_im.set_xlabel('x')
    plt_im.axis([-0.1, 6, -2, 4])
    plt_im.plot(x_im, u_im, color='purple', label='u(x)')
    plt_im.plot(x_im, v_im, color='blue', label='v(x)')
    plt_im.plot(x_data, u_solution, color='purple', label='u(x) sol', linestyle='dashed')
    plt_im.plot(x_data, v_solution, color='blue', label='v(x) sol', linestyle='dashed')

    fig.subplots_adjust(hspace=0.3, wspace=0.3)
    plt_ex.legend()
    plt_im.legend()
    plt.savefig('task8/test8_exp1.png')