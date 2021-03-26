

from time import time

class Timer():

    def __init__(self):
        self.t0=time()

    def reset(self):
        self.t0=time()

    def check(self):
        t = time() - self.t0
        self.reset()
        return t