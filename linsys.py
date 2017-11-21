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
            self.wc = 0

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
        num_equations = len(system)
        num_variables = system.dimension

        col = 0
        for row in range(num_equations):
            while col < num_variables:
                c = MyDecimal(system[row].normal_vector[col])
                if c.is_near_zero():
                    swap_succeeded = system.swap_with_row_below_for_nonzero_coefficient_if_able(row, col)
                    if not swap_succeeded:
                        col += 1
                        continue

                system.clear_coefficients_below(row, col)
                col += 1
                break

        return system

    def swap_with_row_below_for_nonzero_coefficient_if_able(self, given_row, col):
        for row in range(given_row + 1, len(self)):
            c2 = MyDecimal(self[row][col])
            if not c2.is_near_zero():
                self.swap_rows(given_row, row)
                return True
        return False

    def clear_coefficients_below(self, given_row, col):
        c1 = MyDecimal(self[given_row][col])
        for row in range(given_row + 1, len(self)):
            c2 = MyDecimal(self[row][col])
            c = -c2/c1
            self.add_multiple_times_row_to_row(c, given_row, row)

    def clear_coefficients_above(self, given_row, col):
        for row in range(given_row)[::-1]:
            c2 = MyDecimal(-self[row][col])
            self.add_multiple_times_row_to_row(c2, given_row, row)

    def compute_rref(self):
        tf = self.compute_triangular_form()
        pivot_indices = tf.indices_of_first_nonzero_terms_in_each_row()
        for row in range(len(tf))[::-1]:
            pivot_var = pivot_indices[row]
            if pivot_var < 0:
                continue
            tf.scale_row_to_make_coefficient_equal_one(row, pivot_var)
            tf.clear_coefficients_above(row, pivot_var)
        return tf

    def scale_row_to_make_coefficient_equal_one(self, row, col):
        n = self[row].normal_vector
        beta = Decimal('1.0') / n[col]
        self.multiply_coefficient_and_row(beta, row)

    def compute_triangular_form_original(self):
        system = deepcopy(self)
        still_working = True
        print("compute_triangular_form_original_BEGIN")
        print(self)

        while len(self) > 3:
            # Let's reduce this to 3 equations
            # But we don't want to change those
            # planes with some zero coefficients 
            # because they are already candidates
            # for the final triangle.
            # So let's get the planes with the most
            # non-zero coefficients. 
            # First we'll get a sorted list of 
            # indices, from least to most, whose
            # corresponding plane has so many
            # nonzero coefficients
            L = self.sorted_indices_num_nonzero_terms_each_row_min_to_max()
            # Now we get the list of indices we wish to condense
            # The indices we want are at the tail of this list
            # because the list is least to most.
            idxs_to_condense = L[(2 - len(self)):]
            # Now we get the list of planes we wish to condense
            planes_to_condense = [self[i] for i in idxs_to_condense]

            # now we condense em
            for i,p in enumerate(planes_to_condense):
                the_condensed_plane = p if i == 0 else p.add(the_condensed_plane)

            for idx in idxs_to_condense:
                self.planes.pop(idx)

            # And finally, let's add the_condensed_plane to the list
            self.planes.append(the_condensed_plane)

        while still_working:
            n0x = self.indices_of_first_nonzero_terms_in_each_row()
            print("n0x ({}) ({})".format(self.wc, n0x))
            self.wc += 1
            if self.wc >= 5: return system, False
            if 0 in n0x and sum(n0x) == 3:
                # so we have some combination of [0,1,2]
                # first let's make it exactly [0,1,2]
                if n0x[0] != 0:
                    if n0x[1] == 0:
                        self.swap_rows(0,1)
                        n0x[0], n0x[1] = n0x[1], n0x[0]
                    else:
                        self.swap_rows(0,2)
                        n0x[0], n0x[2] = n0x[2], n0x[0]
                if n0x[1] != 1:
                    self.swap_rows(1,2)
                    n0x[1], n0x[2] = n0x[2], n0x[1]
                # and now we have the desired triangular form
                w0x = self.indices_of_first_nonzero_terms_in_each_row()
                print("www000xxx ({})".format(w0x))
                print("nnn000xxx ({})".format(n0x))
                print(self)
                return self, True
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
                    print("s21a.{}".format(self.indices_of_first_nonzero_terms_in_each_row()))
                    print(self)
                    # now let's make [0, 1, 2]
                    p2 = self[1]
                    p3 = self[2]
                    y2,y3 = p2[1], p3[1]
                    coefficient = -y3/y2
                    self.add_multiple_times_row_to_row(coefficient, 1, 2)
                    # and the next loop will see [0, 1, 2]
                    print("s21c.{}".format(self.indices_of_first_nonzero_terms_in_each_row()))
                    print(self)
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
                    print("s22.{}".format(self.indices_of_first_nonzero_terms_in_each_row()))
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
                print("07")
                print(self)
            else:
                # all zeros [0, 0, 0]
                p1 = self[0]
                p3 = self[2]
                x1,x3 = p1[0], p3[0]
                coefficient = -x3/x1
                self.add_multiple_times_row_to_row(coefficient, 0, 2)
                # and the next loop will see [0, 0, 1]

        return system, False

    def sorted_indices_num_nonzero_terms_each_row_min_to_max(self):
        L = self.num_nonzero_terms_each_row()
        print("LLLLLLLLLLLLLLLLLLLLL ({})".format(L))
        return sorted(range(len(L)), key=lambda i:L[i])

    def num_nonzero_terms_each_row(self):
        lst = []
        for p in self.planes:
            lst.append(sum([1 if c != 0 else 0 for c in p]))
        return lst

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
    print("TEST ROW OPS")
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
    print("TEST TRIANGULAR FORM")
    p1 = Plane(normal_vector=Vector(['1','1','1']), constant_term='1')
    p2 = Plane(normal_vector=Vector(['0','1','1']), constant_term='2')
    s = LinearSystem([p1,p2])
    t = s.compute_triangular_form()
    if not (t[0] == p1 and
            t[1] == p2):
        print('test case 1 failed')
    else: print('test case 1 PASSED')
    
    p1 = Plane(normal_vector=Vector(['1','1','1']), constant_term='1')
    p2 = Plane(normal_vector=Vector(['1','1','1']), constant_term='2')
    s = LinearSystem([p1,p2])
    t = s.compute_triangular_form()
    if not (t[0] == p1 and
            t[1] == Plane(constant_term='1')):
        print('test case 2 failed')
    else: print('test case 2 PASSED')
    
    p1 = Plane(normal_vector=Vector(['1','1','1']), constant_term='1')
    p2 = Plane(normal_vector=Vector(['0','1','0']), constant_term='2')
    p3 = Plane(normal_vector=Vector(['1','1','-1']), constant_term='3')
    p4 = Plane(normal_vector=Vector(['1','0','-2']), constant_term='2')
    s = LinearSystem([p1,p2,p3,p4])
    t = s.compute_triangular_form()
    if not (t[0] == p1 and
            t[1] == p2 and
            t[2] == Plane(normal_vector=Vector(['0','0','-2']), constant_term='2') and
            t[3] == Plane()):
        print('test case 3 failed')
    else: print('test case 3 PASSED')
    
    p1 = Plane(normal_vector=Vector(['0','1','1']), constant_term='1')
    p2 = Plane(normal_vector=Vector(['1','-1','1']), constant_term='2')
    p3 = Plane(normal_vector=Vector(['1','2','-5']), constant_term='3')
    s = LinearSystem([p1,p2,p3])
    t = s.compute_triangular_form()
    if not (t[0] == Plane(normal_vector=Vector(['1','-1','1']), constant_term='2') and
            t[1] == Plane(normal_vector=Vector(['0','1','1']), constant_term='1') and
            t[2] == Plane(normal_vector=Vector(['0','0','-9']), constant_term='-2')):
        print('test case 4 failed')
    else: print('test case 4 PASSED')

def test_triangular_form_original():
    p1 = Plane(normal_vector=Vector(['1','1','1']), constant_term='1')
    p2 = Plane(normal_vector=Vector(['0','1','0']), constant_term='2')
    p3 = Plane(normal_vector=Vector(['1','1','-1']), constant_term='3')
    p4 = Plane(normal_vector=Vector(['1','0','-2']), constant_term='2')
    s = LinearSystem([p1,p2,p3,p4])
    t,w = s.compute_triangular_form()
    if not (t[0] == p1 and
            t[1] == p2 and
            t[2] == Plane(normal_vector=Vector(['0','0','-2']), constant_term='2') and
            t[3] == Plane()):
        print('test case 3 failed')

    p1 = Plane(normal_vector=Vector(['0','1','1']), constant_term='1')
    p2 = Plane(normal_vector=Vector(['1','-1','1']), constant_term='2')
    p3 = Plane(normal_vector=Vector(['1','2','-5']), constant_term='3')
    s = LinearSystem([p1,p2,p3])
    t,w = s.compute_triangular_form()
    print("WWWWWWWWWWWWWWWW t({}) w({})".format(t, w))
    if not (t[0] == Plane(normal_vector=Vector(['1','-1','1']), constant_term='2') and
            t[1] == Plane(normal_vector=Vector(['0','1','1']), constant_term='1') and
            t[2] == Plane(normal_vector=Vector(['0','0','-9']), constant_term='-2')):
        print('test case 4 failed')

def test_rref():
    print("TEST REDUCED ROW ECHELON FORM")
    p1 = Plane(normal_vector=Vector(['1', '1', '1']), constant_term='1')
    p2 = Plane(normal_vector=Vector(['0', '1', '1']), constant_term='2')
    s = LinearSystem([p1, p2])
    r = s.compute_rref()
    if not (r[0] == Plane(normal_vector=Vector(['1', '0', '0']),
                          constant_term='-1') and
            r[1] == p2):
        print('test case 1 failed')
    else: print('test case 2 PASSED')
    
    p1 = Plane(normal_vector=Vector(['1', '1', '1']), constant_term='1')
    p2 = Plane(normal_vector=Vector(['1', '1', '1']), constant_term='2')
    s = LinearSystem([p1, p2])
    r = s.compute_rref()
    if not (r[0] == p1 and
            r[1] == Plane(constant_term='1')):
        print('test case 2 failed')
    else: print('test case 2 PASSED')
    
    p1 = Plane(normal_vector=Vector(['1', '1', '1']), constant_term='1')
    p2 = Plane(normal_vector=Vector(['0', '1', '0']), constant_term='2')
    p3 = Plane(normal_vector=Vector(['1', '1', '-1']), constant_term='3')
    p4 = Plane(normal_vector=Vector(['1', '0', '-2']), constant_term='2')
    s = LinearSystem([p1, p2, p3, p4])
    r = s.compute_rref()
    if not (r[0] == Plane(normal_vector=Vector(['1', '0', '0']),
                          constant_term='0') and
            r[1] == p2 and
            r[2] == Plane(normal_vector=Vector(['0', '0', '-2']),
                          constant_term='2') and
            r[3] == Plane()):
        print('test case 3 failed')
    else: print('test case 3 PASSED')
    
    p1 = Plane(normal_vector=Vector(['0', '1', '1']), constant_term='1')
    p2 = Plane(normal_vector=Vector(['1', '-1', '1']), constant_term='2')
    p3 = Plane(normal_vector=Vector(['1', '2', '-5']), constant_term='3')
    s = LinearSystem([p1, p2, p3])
    r = s.compute_rref()
    if not (r[0] == Plane(normal_vector=Vector(['1', '0', '0']),
                          constant_term=Decimal('23') / Decimal('9')) and
            r[1] == Plane(normal_vector=Vector(['0', '1', '0']),
                          constant_term=Decimal('7') / Decimal('9')) and
            r[2] == Plane(normal_vector=Vector(['0', '0', '1']),
                          constant_term=Decimal('2') / Decimal('9'))):
        print('test case 4 failed')
    else: print('test case 4 PASSED')

def test_rref_xtra():
    p1 = Plane(normal_vector=Vector(['17', '14', '-13']), constant_term='13')
    p2 = Plane(normal_vector=Vector(['-3', '1', '-19']), constant_term='-1')
    p3 = Plane(normal_vector=Vector(['8', '7', '21']), constant_term='11')
    s = LinearSystem([p1, p2, p3])
    r = s.compute_rref()
    print("test_rref_xtra")
    print(r)


def test():
    test_row_ops()
#    test_triangular_form_original()
    test_triangular_form()
    test_rref()
    test_rref_xtra()

if __name__ == '__main__':
    test()
