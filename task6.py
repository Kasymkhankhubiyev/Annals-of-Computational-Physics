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


def difur(x):
    res = -x
    return res


def func(t):
    x = np.float64(np.exp(-t))
    return x


def euler():
    x=[]
    y=[]

def runge_kutta():
    pass

def run():
    pass