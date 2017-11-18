from decimal import Decimal, getcontext

from vector import Vector

getcontext().prec = 30


class Line(object):

    NO_NONZERO_ELTS_FOUND_MSG = 'No nonzero elements found'

    def __init__(self, normal_vector=None, constant_term=None):
        self.dimension = 2

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

            initial_index = Line.first_nonzero_index(n)
            initial_coefficient = n[initial_index]

            basepoint_coords[initial_index] = c/initial_coefficient
            self.basepoint = Vector(basepoint_coords)

        except Exception as e:
            if str(e) == Line.NO_NONZERO_ELTS_FOUND_MSG:
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
            initial_index = Line.first_nonzero_index(n)
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


    @staticmethod
    def first_nonzero_index(iterable):
        for k, item in enumerate(iterable):
            if not MyDecimal(item).is_near_zero():
                return k
        raise Exception(Line.NO_NONZERO_ELTS_FOUND_MSG)

    def is_parallel_with(self, other):
        return self.normal_vector.is_parallel_to(other.normal_vector)

    def intersection_point(self, other):
        if not self.is_parallel_with(other):
            # then we must find the single point of intersection
            A = self.normal_vector[0]
            B = self.normal_vector[1]
            C = other.normal_vector[0]
            D = other.normal_vector[1]
            k1, k2 = self.constant_term, other.constant_term
            det = A*D - B*C
            x = (D * k1 - B * k2) / det
            y = (-C * k1 + A * k2) / det
            return Vector([x, y])
        elif self == other:
            return self
        else:
            return None


class MyDecimal(Decimal):
    def is_near_zero(self, eps=1e-10):
        return abs(self) < eps


def test_intersection():
    v = Vector([4.046, 2.836])
    l1 = Line(normal_vector=v, constant_term = 1.21)
    w = Vector([10.115, 7.09])
    l2 = Line(normal_vector=w, constant_term = 3.025)
    intersect_pt = l1.intersection_point(l2)
    is_parallel = l1.is_parallel_with(l2)
    is_eq = l1 == l2
    print("intersect.1({}) is parallel({}) is equal({})".format(intersect_pt, is_parallel, is_eq))
    
    v = Vector([7.204, 3.182])
    l1 = Line(normal_vector=v, constant_term = 8.68)
    w = Vector([8.172, 4.114])
    l2 = Line(normal_vector=w, constant_term = 9.883)
    intersect_pt = l1.intersection_point(l2)
    is_parallel = l1.is_parallel_with(l2)
    is_eq = l1 == l2
    print("intersect.2({}) is parallel({}) is equal({})".format(intersect_pt, is_parallel, is_eq))

    v = Vector([1.182, 5.562])
    l1 = Line(normal_vector=v, constant_term = 6.744)
    w = Vector([1.773, 8.343])
    l2 = Line(normal_vector=w, constant_term = 9.525)
    intersect_pt = l1.intersection_point(l2)
    is_parallel = l1.is_parallel_with(l2)
    is_eq = l1 == l2
    print("intersect.3({}) is parallel({}) is equal({})".format(intersect_pt, is_parallel, is_eq))

def test():
    test_intersection()


if __name__ == '__main__':
    test()
