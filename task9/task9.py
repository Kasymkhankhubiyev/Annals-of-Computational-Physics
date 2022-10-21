"""
Методом погонки решить разностный аналог граничной задачи
для уравнения y''=cosx на промежутке -pi/2 < x < pi/2.

Рассмотреть различный варианты граничных условий (га функцию и ее производную).
"""

import numpy as np
import matplotlib.pyplot as plt

Constant = float


def second_diff(x: float) -> float:
    res = np.cos(x)
    return res


def first_diff(x: float, c0: Constant):
    """
    y' = sinx + c0
    """
    res = np.sin(x) + c0
    return res


def func(x: float, c0: Constant, c1: Constant):
    """
    y = -cosx + c0*x + c1
    """
    res = -1 * np.cos(x) + c0 * x + c1 * x
    return res


def run() -> None:
    pass



