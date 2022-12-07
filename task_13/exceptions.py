
class CantMatchMethod(Exception):
    def __init__(self, message, methods):
        super().__init__(
            f'Cannot math a method: {message} '
            f'choose one from {methods}'
        )


class CantRunDichotomyMethod(Exception):
    def __init__(self, phi0, phi1):
        super().__init__(
            f'Cannot run the dichotomy method \t'
            f'with wrong borders: [{phi0}, {phi1}]\t'
            f'req: [a, b]: a < b'
        )
