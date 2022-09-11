import numpy as np

class SmartSlice:
    """slice-like object that support basic airthmatic operation +. -, *.
    start, stop and step should be non-negative integers and cannot be None.
    """
    def __init__(self, start, stop, step=1):
        assert start is not None
        assert stop is not None
        assert step is not None

        self._slice = slice(start, stop, step)

    @property
    def start(self):
        return self._slice.start

    @property
    def stop(self):
        return self._slice.stop

    @property
    def step(self):
        return self._slice.step

    def __mul__(self, i):
        start = self.start * i
        stop = self.stop * i
        step = self.step * i

        return type(self)(start, stop, step)

    def __rmul__(self, i):
        return self * i

    def __add__(self, i):
        start = self.start + i
        stop = self.stop + i
        step = self.step

        return type(self)(start, stop, step)

    def __radd__(self, i):
        return self + i

    def __sub__(self, i):
        return self.__add__(-i)

    def __rsub__(self, i):
        return self - i


def test():
    t = SmartSlice(1, 10)
    t2 = 2 * t - 1
    print(t2._slice)


class NDArraySmartSlice(np.ndarray):
    """ numpy NDArray that support indexing by SmartSlice.
    """
    def __getitem__(self, i):
        if isinstance(i, SmartSlice):
            return super().__getitem__(slice(i.start,i.stop,i.step))
        else:
            return super().__getitem__(i)

    def __setitem__(self, i, v):
        if isinstance(i, SmartSlice):
            return super().__setitem__(slice(i.start,i.stop,i.step), v)
        else:
            return super().__setitem__(i, v)

    @classmethod
    def cast(cls, x):
        return x.view(cls)


def cast(x):
    return NDArraySmartSlice.cast(x)


def test_ndarrayslice():
    # NDArraySlice([1, 1, 1])
    N = 100
    a = np.arange(N).view(NDArraySmartSlice)
    print(a)
    i = SmartSlice(1, N-1, 1)
    print(a[i])
