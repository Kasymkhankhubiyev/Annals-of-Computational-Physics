"""
Сигнал имеет вид:

f(t) = a0*sin(w0*t) + a1*sin(w1*t)
"""

import numpy as np
import matplotlib.pyplot as plt

a0, a1 = 1., 0.002
w0, w1 = 5.1, 25.5
T = 2 * np.pi

time = float
signal = float


def func(t: time) -> signal:
    f = a0 * np.sin(w0 * t) + a1*np.sin(w1 * t)
    return f


def run() -> None:
    pass
