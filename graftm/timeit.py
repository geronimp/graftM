import time

class Timer:
    
    def __init__(self): pass

    def timeit(self, method):
        def timed(*args, **kw):
            ts = time.time()
            result = method(*args, **kw)
            te = time.time()
            return round(te-ts, 2), result
        return timed