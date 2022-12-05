"""
Сигнал имеет вид:

f(t) = a0*sin(w0*t) + a1*sin(w1*t)
"""

import numpy as np
import matplotlib.pyplot as plt
from typing import NamedTuple

a0, a1 = 1., 0.002
w0, w1 = 5.1, 25.5/2.
T, N = 2 * np.pi, 50


class Result(NamedTuple):
    w: np.array
    spectrum: np.array


def func(t: np.array) -> np.array:
    f = a0 * np.sin(w0 * t) + a1*np.sin(w1 * t)
    return f


def square_window_func(k) -> int:
    if N > k >= 0:
        return 1
    return 0


def hann_window_func(k) -> float:
    if N > k >= 0:
        return (1 - np.cos(2 * np.pi * k / N)) / 2


def py_fftw(t: np.array, h):
    fr = 0
    for i in range(N):
        fr = np.fft.fft(func(t)) * h(i)

    return fr * fr.conjugate()


def FFTW(t: np.array, h):
    w, pow_spect = [], []
    for i in range(N):
        fi = complex(0, 0)
        for k in range(N):
            fi += func(t[k]) * np.exp(2 * np.pi * 1j * i * k / N) * h(k)
        pow_spect.append((fi * fi.conjugate()))
        w.append(2 * np.pi * i / T)
    return Result(w=np.array(w), spectrum=np.array(pow_spect))


def run() -> None:
    time = np.linspace(0, T, N)

    # signal = func(time)
    # fig, axs = plt.subplots(nrows=1, ncols=2)
    plt.yscale('log')
    sq_window, hann_window = FFTW(time, square_window_func), FFTW(time, hann_window_func)
    plt.plot(sq_window.w, sq_window.spectrum, color='green', label='square window')
    plt.plot(hann_window.w, hann_window.spectrum, color='red', label='hann window')
    # plt.plot(hann_window.w, py_fftw(time, hann_window_func), color='blue', label='gt')
    plt.legend(fontsize=7, ncol=1, facecolor='oldlace', edgecolor='r')

    # error_sq = np.abs(np.array(py_fftw(time, square_window_func) - sq_window.spectrum))
    # axs[1].plot(sq_window.w, error_sq, color='green', label='gt VS sq_win')
    # error_hann = np.abs(np.array(py_fftw(time, square_window_func) - hann_window.spectrum))
    # axs[1].plot(hann_window.w, error_hann, color='blue', label='gt VS hann_win')
    # axs[1].legend(fontsize=7, ncol=1, facecolor='oldlace', edgecolor='r')
    plt.title(f'w0={w0},  w1={w1}')
    plt.savefig('task_12/signal_spectr.png')
    plt.close()

