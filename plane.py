from decimal import Decimal, getcontext

from vector import Vector

getcontext().prec = 30


class Plane(object):

    NO_NONZERO_ELTS_FOUND_MSG = 'No nonzero elements found'

    def __init__(self, normal_vector=None, constant_term=None):
        self.dimension = 3

        if not normal_vector:
            all_zeros = ['0']*self.dimension
            normal_vector = Vector(all_zeros)
        self.normal_vector = normal_vector

        if not constant_term:
            constant_term = Decimal('0')
        self.constant_term = Decimal(constant_term)

        self.set_basepoint()


    def set_basepoint(self):
        try:
            n = self.normal_vector
            c = self.constant_term
            basepoint_coords = ['0']*self.dimension

            initial_index = Plane.first_nonzero_index(n)
            initial_coefficient = n[initial_index]

            basepoint_coords[initial_index] = c/initial_coefficient
            self.basepoint = Vector(basepoint_coords)

        except Exception as e:
            if str(e) == Plane.NO_NONZERO_ELTS_FOUND_MSG:
                self.basepoint = None
            else:
                raise e

    def __eq__(self, other):
        if self.normal_vector.is_zero():
            if not other.normal_vector.is_zero():
                return False
            else:
                diff = self.constant_term - other.constant_term
                return MyDecimal(diff).is_near_zero()
        elif other.normal_vector.is_zero():
            return False

        if not self.is_parallel_with(other):
            return False
        # calculate connecting vector
        cv = self.basepoint.minus(other.basepoint)
        if cv.is_orthogonal_to(self.normal_vector):
            # It is unnecessary to check that cv is
            # orthogonal to other.normal_vector
            # because we already verified that the
            # lines are parallel. Hence their normal
            # vectors are also parallel. So if cv
            # is orthogonal to one, it's orthogonal
            # to the other.
            return True
        return False

    def __str__(self):

        num_decimal_places = 3

        def write_coefficient(coefficient, is_initial_term=False):
            coefficient = round(coefficient, num_decimal_places)
            if coefficient % 1 == 0:
                coefficient = int(coefficient)

            output = ''

            if coefficient < 0:
                output += '-'
            if coefficient > 0 and not is_initial_term:
                output += '+'

            if not is_initial_term:
                output += ' '

            if abs(coefficient) != 1:
                output += '{}'.format(abs(coefficient))

            return output

        n = self.normal_vector

        try:
            initial_index = Plane.first_nonzero_index(n)
            terms = [write_coefficient(n[i], is_initial_term=(i==initial_index)) + 'x_{}'.format(i+1)
                     for i in range(self.dimension) if round(n[i], num_decimal_places) != 0]
            output = ' '.join(terms)

        except Exception as e:
            if str(e) == self.NO_NONZERO_ELTS_FOUND_MSG:
                output = '0'
            else:
                raise e

        constant = round(self.constant_term, num_decimal_places)
        if constant % 1 == 0:
            constant = int(constant)
        output += ' = {}'.format(constant)

        return output

    def is_parallel_with(self, other):
        return self.normal_vector.is_parallel_to(other.normal_vector, tolerance=1e-4)

    @staticmethod
    def first_nonzero_index(iterable):
        for k, item in enumerate(iterable):
            if not MyDecimal(item).is_near_zero():
                return k
        raise Exception(Plane.NO_NONZERO_ELTS_FOUND_MSG)


class MyDecimal(Decimal):
    def is_near_zero(self, eps=1e-10):
        return abs(self) < eps

def test_parallel_equal():
    v = Vector([-0.412, 3.806, 0.728])
    p1 = Plane(normal_vector=v, constant_term=-3.46)
    w = Vector([1.03, -9.515, -1.82])
    p2 = Plane(normal_vector=w, constant_term=8.65)
    print("parallel.1({}) equal({})".format(p1.is_parallel_with(p2), p1==p2))

    v = Vector([2.611, 5.528, 0.283])
    p1 = Plane(normal_vector=v, constant_term=4.6)
    w = Vector([7.715, 8.306, 5.342])
    p2 = Plane(normal_vector=w, constant_term=3.76)
    print("parallel.2({}) equal({})".format(p1.is_parallel_with(p2), p1==p2))

    v = Vector([-7.926, 8.625, -7.217])
    p1 = Plane(normal_vector=v, constant_term=-7.952)
    w = Vector([-2.642, 2.875, -2.405666])
    p2 = Plane(normal_vector=w, constant_term=-2.651)
    wayne = [round(c/3,3) for c in v]
    print("parallel.3({}) equal({})".format(p1.is_parallel_with(p2), p1==p2))

def test():
    test_parallel_equal()

if __name__ == '__main__':
    test()
