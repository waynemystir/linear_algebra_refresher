from math import sqrt, acos, pi, sin
from decimal import Decimal, getcontext

getcontext().prec = 30

class Vector(object):
    
    CANNOT_NORMALIZE_ZERO_VECTOR_MSG = 'Cannot normalize the zero vector'
    NO_UNIQUE_PARALLEL_COMPONENT_TO_MSG = 'Cannot find a unique component parallel to the given basis vector'
    NO_UNIQUE_ORTHOGONAL_COMPONENT_TO_MSG = 'Cannot find a unique component orthogonal to the given basis vector'
    CROSS_PRODUCT_WRONG_NUM_DIMS = 'Cross product wrong number of dimensions, should be 3'

    def __init__(self, coordinates):
        try:
            if not coordinates:
                raise ValueError
            self.coordinates = tuple([Decimal(x) for x in coordinates])
            self.dimension = len(coordinates)

        except ValueError:
            raise ValueError('The coordinates must be nonempty')

        except TypeError:
            raise TypeError('The coordinates must be an iterable')

    def __iter__(self):
        self.current = 0
        return self

    def __next__(self):
        if self.current >= len(self.coordinates):
            raise StopIteration
        else:
            current_value = self.coordinates[self.current]
            self.current += 1
            return current_value

    def __getitem__(self, i):
        return self.coordinates[i]

    def __str__(self):
        return 'Vector: {}'.format([round(c,3) for c in self.coordinates])

    def __eq__(self, v):
        return self.coordinates == v.coordinates

    def __len__(self):
        return len(self.coordinates)

    def plus(self, v):
        if type(v) == list: v = Vector(v)
        if len(self) != len(v):
            raise Exception("plus function: vector"
                   " lengths not equal (%d) != (%d)" % (len(self), len(b)))
        return Vector([s+v for s,v in zip(self.coordinates, v.coordinates)])

    def minus(self, v):
        if type(v) == list: v = Vector(v)
        if len(self) != len(v):
            raise Exception("minus function: vector"
                   " lengths not equal (%d) != (%d)" % (len(self), len(b)))
        return Vector([Decimal(s-v) for s,v in zip(self.coordinates, v.coordinates)])

    def scalar_multiply(self, sc):
        return Vector([Decimal(sc)*Decimal(s) for s in self.coordinates])

    def magnitude(self):
        return Decimal(sqrt(sum([Decimal(c**2) for c in self.coordinates])))

    def normalized(self):
        try:
            return self.scalar_multiply(Decimal('1.0')/self.magnitude())
        except ZeroDivisionError:
            raise Exception(self.CANNOT_NORMALIZE_ZERO_VECTOR_MSG)

    def dot(self, v):
        if type(v) == list: v = Vector(v)
        if len(self) != len(v):
            raise Exception("dot function: vector"
                   " lengths not equal (%d) != (%d)" % (len(self), len(b)))
        return Decimal(sum([Decimal(s*v) for s,v in zip(self.coordinates, v.coordinates)]))

    def angle_with(self, v, in_degrees=False, tolerance=1e-10):
        try:
            sn = self.normalized()
            vn = v.normalized()
            dp = sn.dot(vn)
            if dp > Decimal('1.0') and abs(dp - Decimal('1.0')) < tolerance:
                dp = Decimal('1.0')
            if dp < Decimal('-1.0') and abs(dp - Decimal('-1.0')) < tolerance:
                dp = Decimal('-1.0')
            air = Decimal(acos(dp))
        except Exception as e:
            if str(e) == self.CANNOT_NORMALIZE_ZERO_VECTOR_MSG:
                raise Exception('Cannot compute an angle with the zero vector')
            else: raise e
        if in_degrees: return air * Decimal(180. / pi)
        else: return air

    def is_parallel_to_original(self, v, tolerance=1e-10):
        if len(self) != len(v):
            raise Exception("is_parallel_to function: vector"
                   " lengths not equal (%d) != (%d)" % (len(self), len(b)))
        lst = [0 if v == 0 else Decimal(s/v) for s,v in zip(self.coordinates, v.coordinates)]
        for l in lst:
            if abs(l - lst[0]) > tolerance:
                return False
        return True

    def is_parallel_to(self, v, tolerance=1e-10):
        if len(self) != len(v):
            raise Exception("is_parallel_to function: vector"
                   " lengths not equal (%d) != (%d)" % (len(self), len(b)))
        if self.is_zero() or v.is_zero(): return True
        air = self.angle_with(v)
        return abs(air) < tolerance or abs(air - Decimal(pi)) < tolerance

    def is_zero(self, tolerance=1e-10):
        return self.magnitude() < tolerance

    def is_orthogonal_to(self, v, tolerance=1e-10):
        if len(self) != len(v):
            raise Exception("is_orthogonal_to function: vector"
                   " lengths not equal (%d) != (%d)" % (len(self), len(b)))
        return abs(self.dot(v)) < tolerance

    def projection_to(self, b):
        if len(self) != len(b):
            raise Exception("projection_to function: vector"
                   " lengths not equal (%d) != (%d)" % (len(self), len(b)))

        try:
            ub = b.normalized()
            parallellen = self.dot(ub)
            projection = ub.scalar_multiply(parallellen)
            return projection

        except Exception as e:
            if str(e) == self.CANNOT_NORMALIZE_ZERO_VECTOR_MSG:
                raise Exception(self.NO_UNIQUE_PARALLEL_COMPONENT_TO_MSG)
            else: raise e

    def orthogonal_component_to(self, b):
        if len(self) != len(b):
            raise Exception("orthogonal_component_to function: vector"
                   " lengths not equal (%d) != (%d)" % (len(self), len(b)))

        try:
            projection = self.projection_to(b)
            return self.minus(projection)

        except Exception as e:
            if str(e) == self.NO_UNIQUE_PARALLEL_COMPONENT_TO_MSG:
                raise Exception(self.NO_UNIQUE_ORTHOGONAL_COMPONENT_TO_MSG)
            else: raise e

    def cross_product(self, w):
        if len(self) != len(w):
            raise Exception("cross_product function: vector"
                   " lengths not equal (%d) != (%d)" % (len(self), len(b)))
        if self.dimension != 3 or w.dimension !=3:
            raise Exception(self.CROSS_PRODUCT_WRONG_NUM_DIMS)
        x = self.coordinates[1] * w.coordinates[2] - w.coordinates[1] * self.coordinates[2]
        y = -(self.coordinates[0] * w.coordinates[2] - w.coordinates[0] * self.coordinates[2])
        z = self.coordinates[0] * w.coordinates[1] - w.coordinates[0] * self.coordinates[1]
        return Vector([x,y,z])

    def area_of_parallelogram(self, w):
        return self.cross_product(w).magnitude()
#        base = self.magnitude()
#        θ = self.angle_with(w)
#        height = Decimal(sin(θ)) * w.magnitude()
#        return base * height

    def area_of_triangle(self, w):
        return Decimal(1/2) * self.area_of_parallelogram(w)


def test_plus_minus_scalar_mult():
    v1 = Vector([8.218, -9.341])
    v2 = Vector([-1.129, 2.111])
    vplus = v1.plus(v2)
    v1 = Vector([7.119, 8.215])
    v2 = Vector([-8.223, 0.878])
    vminus = v1.minus(v2)
    v = Vector([1.671, -1.012, -0.318])
    sm = v.scalar_multiply(7.41)
    print("plus {}".format(vplus))
    print("minus {}".format(vminus))
    print("scalar_multiply {}".format(sm))

def test_magnitude_direction():
    v = Vector([-0.221, 7.437])
    print("v.magnitude.1: {}".format(round(v.magnitude(), 3)))
    v = Vector([8.813, -1.331, -6.247])
    print("v.magnitude.2: {}".format(round(v.magnitude(), 3)))
    v = Vector([0, 0])
    try:
        print("v.direction.0: {} is it unit: {}".format(v.normalized(), v.normalized().magnitude()))
    except Exception as e:
        print("Error occurred with the zero vector: %s" % str(e))
    v = Vector([5.581, -2.136])
    print("v.direction.1: {} is it unit: {}".format(v.normalized(), v.normalized().magnitude()))
    v = Vector([1.996, 3.108, -4.554])
    print("v.direction.2: {} is it unit: {}".format(v.normalized(), v.normalized().magnitude()))

def test_dot_angle():
    v = Vector([7.887, 4.138])
    w = Vector([-8.802, 6.776])
    print("dot-prod-1 {}".format(round(v.dot(w), 3)))
    v = Vector([-5.955, -4.904, -1.874])
    w = Vector([-4.496, -8.755, 7.103])
    print("dot-prod-2 {}".format(round(v.dot(w), 3)))
    v = Vector([3.183, -7.627])
    w = Vector([-2.668, 5.319])
    print("angle_with-1 {}".format(round(v.angle_with(w), 3)))
    v = Vector([7.35, 0.221, 5.188])
    w = Vector([0.,0.,0.])
    try:
        v.angle_with(w)
    except Exception as e:
        print("We tried to get an angle with the zero vector".format(e))
    w = Vector([2.751, 8.259, 3.985])
    print("angle_with-2 {}".format(round(v.angle_with(w, in_degrees=True), 3)))

def test_parallel_orthogonal():
    v = Vector([-7.579, -7.88])
    w = Vector([22.737, 23.64])
    print("is_para_orth_1 ({}) ({})".format(v.is_parallel_to(w), v.is_orthogonal_to(w)))
    v = Vector([-2.029, 9.97, 4.172])
    w = Vector([-9.231, -6.639, -7.245])
    print("is_para_orth_2 ({}) ({})".format(v.is_parallel_to(w), v.is_orthogonal_to(w)))
    v = Vector([-2.328, -7.284, -1.214])
    w = Vector([-1.821, 1.072, -2.94])
    print("is_para_orth_3 ({}) ({})".format(v.is_parallel_to(w), v.is_orthogonal_to(w)))
    v = Vector([2.118, 4.827])
    w = Vector([0, 0])
    print("is_para_orth_4 ({}) ({})".format(v.is_parallel_to(w), v.is_orthogonal_to(w)))

def test_projections():
    v = Vector([3.039, 1.879])
    b = Vector([0, 0])
    try:
        print("projection_to.1 {}".format(v.projection_to(b)))
    except Exception as e:
        print("test_projections_zero_vector_exception: (%s)" % (e))
    try:
        print("orthogonal_component_to.1 {}".format(v.orthogonal_component_to(b)))
    except Exception as e:
        print("test_orthogonal_comps_zero_vector_exception: (%s)" % (e))

    v = Vector([3.039, 1.879])
    b = Vector([0.825, 2.036])
    print("projection_to.1 {}".format(v.projection_to(b)))

    v = Vector([-9.88, -3.264, -8.159])
    b = Vector([-2.155, -9.353, -9.473])
    print("orthogonal_component_to.1 {}".format(v.orthogonal_component_to(b)))

    v = Vector([3.009, -6.172, 3.692, -2.51])
    b = Vector([6.404, -9.144, 2.759, 8.718])
    vproj = v.projection_to(b)
    print("projection_to.2 {}".format(vproj))
    voc = v.orthogonal_component_to(b)
    print("orthogonal_component_to.2 {}".format(v.orthogonal_component_to(b)))
    vc = vproj.plus(voc)
    diff = v.minus(vc)
    print("diff={}".format(diff))

def test_cross_prod_area():
    v = Vector([8.462, 7.893, -8.187])
    w = Vector([6.984, -5.975, 4.778])
    print("cross_product.1 {}".format(v.cross_product(w)))

    v = Vector([-8.987, -9.838, 5.031])
    w = Vector([-4.268, -1.861, -8.866])
    print("area_of_parallelogram.1 {}".format(round(v.area_of_parallelogram(w), 3)))

    v = Vector([1.5, 9.547, 3.691])
    w = Vector([-6.007, 0.124, 5.772])
    print("area_of_triangle.1 {}".format(round(v.area_of_triangle(w), 3)))


def test_vector():
    test_plus_minus_scalar_mult()
    test_magnitude_direction()
    test_dot_angle()
    test_parallel_orthogonal()
    test_projections()
    test_cross_prod_area()

if __name__ == "__main__":
    test_vector()
