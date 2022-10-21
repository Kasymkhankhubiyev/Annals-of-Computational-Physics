"""
решить жесткую систему уравнений

1) u' = 998*u + 1998*v
2) v' = -999*u - 1999*v

"""

import matplotlib.pyplot as plt
import numpy as np


def dif_u(u: float, v:float) -> float:

    res = np.float(998) * u + np.float(1998) * v
    return res


def dif_v(u: float, v: float) -> float:

    res = np.float(-999) * u - np.float(1999) * v
    return res


def explicit_first_scheme() -> None:
    """
    y_n+1 = y_n + h*f(x_n, y_n)
    """
    pass


def implicit_first_scheme() -> None:
    """
    y_n+1 = y_n + h/2 * f(x_n, y_n) + h/2 * f(x_n+1, y_n+1)
    """
    pass


def explicit_second_scheme() -> None:
    """
    y_n+3 = y_n+2 + h/12*[23*f(x_n+2, y_n+2) - 16*f(x_n+1, y_n+1) + 5*f(x_n, y_n)]
    """
    pass


def implicit_second_scheme() -> None:
    """
    y_n+2 = y_n+1 + h/12*[5*f(x_n+2, y_n+2) + 8*f(x_n+1, y_n+1) - f(x_n, y_n)]
    """
    pass


def explicit_third_scheme() -> None:
    """
    y_n+4 = y_n+3 + h/24*[55*f(x_n+3, y_n+3) - 59*f(x_n+2, y_n+2) + 37*f(x_n+1, y_n+1) - 9*f(x_n, y_n)]
    """
    pass


def implicit_third_scheme() -> None:
    """
    y_n+3 = y_n+2 + h/24*[9*f(x_n+3, y_n+3) + 19*f(x_n+2, y_n+2) - 5*f(x_n+1, y_n+1) + f(x_n, y_n)]
    """
    pass


def run() -> None:
    pass
