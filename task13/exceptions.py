
class CantMatchMethod(Exception):
    def __init__(self, message, methods):
        super().__init__(
            f'Cannot math a method: {message}'
            f'choose one from {methods}'
        )
