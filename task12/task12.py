"""
Сигнал имеет вид:

f(t) = a0*sin(w0*t) + a1*sin(w1*t)
"""

import numpy as np
import matplotlib.pyplot as plt

a0, a1 = 1., 0.002
w0, w1 = 5.1, 25.5
T = 2 * np.pi


def func(t: np.array) -> np.array:
    f = a0 * np.sin(w0 * t) + a1*np.sin(w1 * t)
    return f


def run() -> None:
    time = np.linspace(0, T, 1000)
    signal = func(time)
    plt.plot(time, signal, color='green')
    plt.legend('best')
    plt.savefig('task12/signal_graph.png')
    plt.close()

