from decimal import Decimal, getcontext
from copy import deepcopy

from vector import Vector
from plane import Plane

getcontext().prec = 30


class LinearSystem(object):

    ALL_PLANES_MUST_BE_IN_SAME_DIM_MSG = 'All planes in the system should live in the same dimension'
    NO_SOLUTIONS_MSG = 'No solutions'
    INF_SOLUTIONS_MSG = 'Infinitely many solutions'

    def __init__(self, planes):
        try:
            d = planes[0].dimension
            for p in planes:
                assert p.dimension == d

            self.planes = planes
            self.dimension = d

        except AssertionError:
            raise Exception(self.ALL_PLANES_MUST_BE_IN_SAME_DIM_MSG)


    def swap_rows(self, row1, row2):
        self[row1], self[row2] = self[row2], self[row1]

    def multiply_coefficient_and_row(self, coefficient, row):
        nv = self[row].normal_vector.scalar_multiply(coefficient)
        ct = coefficient * self[row].constant_term
        self[row] = Plane(normal_vector=nv, constant_term=ct)

    def add_multiple_times_row_to_row(self, coefficient, row_to_add, row_to_be_added_to):
        nv = self[row_to_add].normal_vector.scalar_multiply(coefficient)
        ct = coefficient * self[row_to_add].constant_term
        w = Plane(normal_vector=nv, constant_term=ct)
        e = self[row_to_be_added_to]
        self[row_to_be_added_to] = w.add(e)

    def compute_triangular_form(self):
        system = deepcopy(self)
        still_working = True

        while still_working:
            n0x = self.indices_of_first_nonzero_terms_in_each_row()
            if 0 in n0x and sum(n0x) == 3:
                # so we have some combination of [0,1,2]
                # first let's make it exactly [0,1,2]
                if n0x[0] != 0:
                    if n0x[1] == 0:
                        self.swap_rows(0,1)
                    else:
                        self.swap_rows(0,2)
                if n0x[1] != 1:
                    self.swap_rows(1,2)
                # and now we have the desired triangular form
                return system, True
            elif 0 not in n0x:
                # none of the planes has all the
                # variables
                return system, False
            elif sum(n0x) > 3:
                # so the sum must be 4: [0, 2, 2]
                return system, False
            elif sum(n0x) == 2:
                # so some combination of [0, 1, 1] or [0, 2, 0]
                if 1 in n0x:
                    # first let's make it exactly [0, 1, 1]
                    if n0x[0] != 0:
                        self.swap_rows(0,1)
                        n0x[0], n0x[1] = n0x[1], n0x[0]
                    # now let's make [0, 1, 2]
                    p1 = self[0]
                    p3 = self[2]
                    x1,x3 = p1[0], p3[0]
                    coefficient = -x3/x1
                    self.add_multiple_times_row_to_row(coefficient, 0, 2)
                    # and the next loop will see [0, 1, 2]
                else:
                    # first let's make it exactly [0, 0, 2]
                    if n0x[2] != 2:
                        if n0x[0] == 2:
                            self.swap_rows(0,2)
                            n0x[0], n0x[2] = n0x[2], n0x[0]
                        else:
                            self.swap_rows(1,2)
                            n0x[1], n0x[2] = n0x[2], n0x[1]
                    # now let's make [0, 1, 2]
                    p1 = self[0]
                    p2 = self[1]
                    x1,x2 = p1[0], p2[0]
                    coefficient = -x2/x1
                    self.add_multiple_times_row_to_row(coefficient, 0, 1)
                    # and the next loop will see [0, 1, 2]
            elif sum(n0x) == 1:
                # some combination of [0, 0, 1]
                # first let's make it exactly [0, 0, 1]
                if n0x[0] != 0:
                    self.swap_rows(0,1)
                    n0x[0], n0x[1] = n0x[1], n0x[0]
                if n0x[1] == 1:
                    self.swap_rows(1,2)
                    n0x[1], n0x[2] = n0x[2], n0x[1]
                # now let's make it [0, 1, 1]
                print("06")
                print(self)
                p1 = self[0]
                p2 = self[1]
                x1,x2 = p1[0], p2[0]
                coefficient = -x2/x1
                print("108 - tp1({}) tp2({}) x1({}) x2({}) coeff({})".format(type(p1), type(p2), x1, x2, coefficient))
                self.add_multiple_times_row_to_row(coefficient, 0, 1)
                # and the next loop will see [0, 1, 1]
            else:
                # all zeros [0, 0, 0]
                p1 = self[0]
                p3 = self[2]
                x1,x3 = p1[0], p3[0]
                coefficient = -x3/x1
                self.add_multiple_times_row_to_row(coefficient, 0, 2)
                # and the next loop will see [0, 0, 1]

        return system, False

    def indices_of_first_nonzero_terms_in_each_row(self):
        num_equations = len(self)
        num_variables = self.dimension

        indices = [-1] * num_equations

        for i,p in enumerate(self.planes):
            try:
                indices[i] = p.first_nonzero_index(p.normal_vector)
            except Exception as e:
                if str(e) == Plane.NO_NONZERO_ELTS_FOUND_MSG:
                    continue
                else:
                    raise e

        return indices


    def __len__(self):
        return len(self.planes)


    def __getitem__(self, i):
        return self.planes[i]


    def __setitem__(self, i, x):
        try:
            assert x.dimension == self.dimension
            self.planes[i] = x

        except AssertionError:
            raise Exception(self.ALL_PLANES_MUST_BE_IN_SAME_DIM_MSG)


    def __str__(self):
        ret = 'Linear System:\n'
        temp = ['Equation {}: {}'.format(i+1,p) for i,p in enumerate(self.planes)]
        ret += '\n'.join(temp)
        return ret


class MyDecimal(Decimal):
    def is_near_zero(self, eps=1e-10):
        return abs(self) < eps


p0 = Plane(normal_vector=Vector(['1','1','1']), constant_term='1')
p1 = Plane(normal_vector=Vector(['0','1','0']), constant_term='2')
p2 = Plane(normal_vector=Vector(['1','1','-1']), constant_term='3')
p3 = Plane(normal_vector=Vector(['1','0','-2']), constant_term='2')

s = LinearSystem([p0,p1,p2,p3])

print(s.indices_of_first_nonzero_terms_in_each_row())
print('s0:{},s1:{},s2:{},s3:{}'.format(s[0],s[1],s[2],s[3]))
print(len(s))
print(s)

s[0] = p1
print(s)

print(MyDecimal('1e-9').is_near_zero())
print(MyDecimal('1e-11').is_near_zero())

def test_row_ops():
    p0 = Plane(normal_vector=Vector(['1','1','1']), constant_term='1')
    p1 = Plane(normal_vector=Vector(['0','1','0']), constant_term='2')
    p2 = Plane(normal_vector=Vector(['1','1','-1']), constant_term='3')
    p3 = Plane(normal_vector=Vector(['1','0','-2']), constant_term='2')
    
    s = LinearSystem([p0,p1,p2,p3])
    s.swap_rows(0,1)
    if not (s[0] == p1 and s[1] == p0 and s[2] == p2 and s[3] == p3):
        print('test case 1 failed! @#$%^&*!@#$%^&*!@#$%^&*!@#$%^&*')
    else: print('test case 1 PASSED')
    
    s.swap_rows(1,3)
    if not (s[0] == p1 and s[1] == p3 and s[2] == p2 and s[3] == p0):
        print('test case 2 failed! @#$%^&*!@#$%^&*!@#$%^&*!@#$%^&*')
    else: print('test case 2 PASSED')
    
    s.swap_rows(3,1)
    if not (s[0] == p1 and s[1] == p0 and s[2] == p2 and s[3] == p3):
        print('test case 3 failed! @#$%^&*!@#$%^&*!@#$%^&*!@#$%^&*')
    else: print('test case 3 PASSED')
    
    s.multiply_coefficient_and_row(1,0)
    if not (s[0] == p1 and s[1] == p0 and s[2] == p2 and s[3] == p3):
        print('test case 4 failed! @#$%^&*!@#$%^&*!@#$%^&*!@#$%^&*')
    else: print('test case 4 PASSED')
    
    s.multiply_coefficient_and_row(-1,2)
    if not (s[0] == p1 and
            s[1] == p0 and
            s[2] == Plane(normal_vector=Vector(['-1','-1','1']), constant_term='-3') and
            s[3] == p3):
        print('test case 5 failed! @#$%^&*!@#$%^&*!@#$%^&*!@#$%^&*')
    else: print('test case 5 PASSED')
    
    s.multiply_coefficient_and_row(10,1)
    if not (s[0] == p1 and
            s[1] == Plane(normal_vector=Vector(['10','10','10']), constant_term='10') and
            s[2] == Plane(normal_vector=Vector(['-1','-1','1']), constant_term='-3') and
            s[3] == p3):
        print('test case 6 failed! @#$%^&*!@#$%^&*!@#$%^&*!@#$%^&*')
    else: print('test case 6 PASSED')
    
    s.add_multiple_times_row_to_row(0,0,1)
    if not (s[0] == p1 and
            s[1] == Plane(normal_vector=Vector(['10','10','10']), constant_term='10') and
            s[2] == Plane(normal_vector=Vector(['-1','-1','1']), constant_term='-3') and
            s[3] == p3):
        print('test case 7 failed! @#$%^&*!@#$%^&*!@#$%^&*!@#$%^&*')
    else: print('test case 7 PASSED')
    
    s.add_multiple_times_row_to_row(1,0,1)
    if not (s[0] == p1 and
            s[1] == Plane(normal_vector=Vector(['10','11','10']), constant_term='12') and
            s[2] == Plane(normal_vector=Vector(['-1','-1','1']), constant_term='-3') and
            s[3] == p3):
        print('test case 8 failed! @#$%^&*!@#$%^&*!@#$%^&*!@#$%^&*')
    else: print('test case 8 PASSED')
    
    s.add_multiple_times_row_to_row(-1,1,0)
    if not (s[0] == Plane(normal_vector=Vector(['-10','-10','-10']), constant_term='-10') and
            s[1] == Plane(normal_vector=Vector(['10','11','10']), constant_term='12') and
            s[2] == Plane(normal_vector=Vector(['-1','-1','1']), constant_term='-3') and
            s[3] == p3):
        print('test case 9 failed! @#$%^&*!@#$%^&*!@#$%^&*!@#$%^&*')
    else: print('test case 9 PASSED')

    s.add_multiple_times_row_to_row(7,2,1)
    if not (s[0] == Plane(normal_vector=Vector(['-10','-10','-10']), constant_term='-10') and
            s[1] == Plane(normal_vector=Vector(['3','4','17']), constant_term='-9') and
            s[2] == Plane(normal_vector=Vector(['-1','-1','1']), constant_term='-3') and
            s[3] == p3):
        print('test case 10 failed! @#$%^&*!@#$%^&*!@#$%^&*!@#$%^&*')
    else: print('test case 10 PASSED')

def test_triangular_form():
    p1 = Plane(normal_vector=Vector(['0','1','1']), constant_term='1')
    p2 = Plane(normal_vector=Vector(['1','-1','1']), constant_term='2')
    p3 = Plane(normal_vector=Vector(['1','2','-5']), constant_term='3')
    s = LinearSystem([p1,p2,p3])
    t = s.compute_triangular_form()
    if not (t[0] == Plane(normal_vector=Vector(['1','-1','1']), constant_term='2') and
            t[1] == Plane(normal_vector=Vector(['0','1','1']), constant_term='1') and
            t[2] == Plane(normal_vector=Vector(['0','0','-9']), constant_term='-2')):
        print('test case 4 failed')

def test():
    test_row_ops()
    test_triangular_form()

if __name__ == '__main__':
    test()
