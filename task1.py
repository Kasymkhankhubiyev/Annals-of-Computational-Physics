import numpy as np
from math import inf


def f_epsilon():
    """
    eps :
        1 + eps/2 = 1
        1 + eps != 1

    eps = 1/(2^count) и eps < 1,
    тогда мантисса = count + неявная единица
    :return:
    """
    eps = np.float32(1)
    runnable = 1
    count = 0
    while(runnable):
        new = eps/np.float32(2)
        count += 1
        if ((np.float32(1) + new - np.float32(1) != np.float32(0))
                and (np.float32(1) + new/np.float32(2) - np.float32(1) == np.float32(0))):
            eps = new
            runnable = 0
        else:
            eps = new

    return eps, count


def d_epsilon():
    eps = np.float64(1)
    runnable = 1
    count = 0
    while(runnable):
        new = eps/np.float64(2)
        count += 1
        if ((np.float64(1) + new - np.float64(1) != np.float64(0))
                and (np.float64(1) + new/np.float64(2) - np.float64(1) == np.float64(0))):
            eps = new
            runnable = 0
        else:
            eps = new

    return eps, count


def f_power(state):
    num = np.float32(1)
    count = 0

    if state == 'MAX':
        base = np.float32(2)
        while num != inf:
            num *= base
            count += 1
        return count - 1
    elif state == 'MIN':
        base = np.float32(1/2)
        while num != np.float32(0):
            num *= base
            count += 1
        return count - 1



def d_power(state):
    num = np.float64(1)
    count = 0

    if state == 'MAX':
        base = np.float64(2)
        while num != inf:
            num *= base
            count += 1
        return count - 1
    elif state == 'MIN':
        base = np.float64(1 / 2)
        while num != np.float32(0):
            num *= base
            count += 1
        return count - 1

def f_compare(eps):
    values= {'1': (np.float32(1)), '1+eps/2': np.float32(1) + eps / np.float32(2),
                  '1+eps': np.float32(1) + eps,
                  '1+eps+eps/2': np.float32(1) + eps + eps / np.float32(2)}

    print(values.items())
    for name1, num1 in values.items():
        for name2, num2 in values.items():
            if num2 > num1:
                print(f'({num2}) > ({num1})   ::  {name2} > {name1}')
            elif num1 > num2:
                print(f'({num2}) < ({num1})   ::  {name2} < {name1}')
            else:
                print(f'({num2}) == ({num1})   ::  {name2} == {name1}')


def d_compare(eps):
    values= {'1': (np.float64(1)), '1+eps/2': np.float64(1) + eps / np.float64(2),
                  '1+eps': np.float64(1) + eps,
                  '1+eps+eps/2': np.float64(1) + eps + eps / np.float64(2)}

    print(values.items())
    for name1, num1 in values.items():
        for name2, num2 in values.items():
            if num2 > num1:
                print(f'({num2}) > ({num1})   ::  {name2} > {name1}')
            elif num1 > num2:
                print(f'({num2}) < ({num1})   ::  {name2} < {name1}')
            else:
                print(f'({num2}) == ({num1})   ::  {name2} == {name1}')


def run():
    #print(np.finfo(np.float32))
    # print(np.finfo(np.float64))

    f_eps, f_m = f_epsilon()
    d_eps, d_m = d_epsilon()
    f_pow_max, f_pow_min = f_power('MAX'), f_power('MIN')
    d_pow_max, d_pow_min = d_power('MAX'), d_power('MIN')

    print(f'Эпсилон для 32 бит:     {f_eps}')
    print(f'Эпсилон для 64 бит:      {d_eps} \n\n')

    print(f'Мантисса для 32 бит:    {f_m}')
    print(f'Мастисса для 64 бит:    {d_m}\n\n')

    print(f'Порядок для 32 бит:     MAX: {f_pow_max}    MIN: {f_pow_min}')
    print(f'Порядок для 64 бит:     MAX: {d_pow_max}    MIN: {d_pow_min}')

    f_compare(f_eps)
    d_compare(d_eps)
