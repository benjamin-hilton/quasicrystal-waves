import scipy as sp
import matplotlib.pyplot as plt
import math

phi = (1 + math.sqrt(5)) / 2

class Shape(object):

    def __init__(self, base, head, left, right):
        self._base = base
        self._head = head
        self._left = left
        self._right = right

    def split(self):
        raise NotImplemented("Not implemented by sublcass")

    def should_split(self):
        raise NotImplemented("Not implemented by sublcass")

    @staticmethod
    def interpolate(p1, p2, ratio):
        dx = p1[0] - p2[0]
        dy = p1[1] - p2[1]

        # length = math.sqrt(dx ** 2 + dy ** 2)

        new_x = p1[0] + ratio * dx
        new_y = p1[1] + ratio * dy

        return new_x, new_y

class Kite(Shape):

    def __init__(self, base, head, left, right):
        Shape.__init__(self, base, head, left, right)

    def should_split(self):
        return True

    def split(self):
        base_right = Shape.interpolate(self._base, self._right, phi / (1 + phi))
        base_left = Shape.interpolate(self._base, self._left, phi / (1 + phi))
        dart_head = Shape.interpolate(self._base, self._head, 1 / (1 + phi))

        dart = Dart(base_head, self._base, base_right, base_left)
        kite_right = Kite(self._head, base_left, base_head, self._left)
        kite_left = Kite(self._head, base_right, self._right, base_head)

        return dart, kite_right, kite_left
        

class Dart(Shape):

    def __init__(self, base, head, left, right):
        Shape.__init__(self, base, head, left, right)

    def should_split(self):
        return True

    def split(self):
        pass


class HalfDart(Shape):

    def __init__():
        pass

    def should_split():
        return False

Kite()
