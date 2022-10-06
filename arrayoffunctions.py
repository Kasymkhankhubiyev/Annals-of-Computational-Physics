
def func1(x, y):
    z = x + y
    return z

def hello_world():
    return 0


def func2(x, y):
    z = x ** 2 + y
    return z


def func3(x, y):
    z = x ** 3 + y
    return z


def run():
    functions = [func1, func2, func3]

    for i in range(10):
        result = []
        for func in functions:
            result.append(func(i, 1))
        print(result)
