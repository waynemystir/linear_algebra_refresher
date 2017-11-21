from decimal import Decimal, getcontext
from copy import deepcopy

from vector import Vector
from hyperplane import Hyperplane

getcontext().prec = 30

def _get_new_plane(coefficient, plane):
    new_normal_vector = plane.normal_vector.scalar_multiply(coefficient)
    return Hyperplane(normal_vector=new_normal_vector,
                      constant_term=coefficient * plane.constant_term)

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
        self[row] = _get_new_plane(coefficient, self[row])

    def add_multiple_times_row_to_row(self, coefficient, row_to_add,
                                      row_to_be_added_to):
        recipient_plane = self[row_to_be_added_to]
        new_plane = _get_new_plane(coefficient, self[row_to_add])
        new_normal_vector = \
            recipient_plane.normal_vector.plus(new_plane.normal_vector)
        constant_term = new_plane.constant_term + recipient_plane.constant_term
        self[row_to_be_added_to] = Hyperplane(normal_vector=new_normal_vector,
                                              constant_term=constant_term)

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

    def do_gaussian_elimination(self):
        rref = self.compute_rref()

        try:
            rref.raise_excepion_if_contradictory_equation()
            rref.raise_excepion_if_too_few_pivots()
        except Exception as e:
            return str(e)

        num_variables = rref.dimension
        solution_coordinates = [rref.planes[i].constant_term
                                for i in range(num_variables)]

        return Vector(solution_coordinates)

    def raise_excepion_if_contradictory_equation(self):
        for plane in self.planes:
            try:
                plane.first_nonzero_index(plane.normal_vector)

            except Exception as e:
                if str(e) == 'No nonzero elements found':
                    constant_term = MyDecimal(plane.constant_term)
                    if not constant_term.is_near_zero():
                        raise Exception(self.NO_SOLUTIONS_MSG)

                else:
                    raise e

    def raise_excepion_if_too_few_pivots(self):
        pivot_indices = self.indices_of_first_nonzero_terms_in_each_row()
        num_pivots = sum([1 if index >= 0 else 0 for index in pivot_indices])
        num_variables = self.dimension

        if num_pivots < num_variables:
            raise Exception(self.INF_SOLUTIONS_MSG)

    def compute_solution(self):
        try:
            return self.do_gaussian_elimination_and_parametrization()

        except Exception as e:
            if str(e) == self.NO_SOLUTIONS_MSG:
                return str(e)
            else:
                raise e

    def do_gaussian_elimination_and_parametrization(self):
        rref = self.compute_rref()
        rref.raise_excepion_if_contradictory_equation()

        direction_vectors = rref.extract_direction_vectors_for_parametrization()  # NOQA
        basepoint = rref.extract_basepoint_for_parametrization()

        return Parametrization(basepoint, direction_vectors)

    def extract_direction_vectors_for_parametrization(self):
        num_variables = self.dimension
        pivot_indices = self.indices_of_first_nonzero_terms_in_each_row()
        free_variable_indices = set(range(num_variables)) - set(pivot_indices)

        direction_vectors = []

        for free_var in free_variable_indices:
            vector_coords = [0] * num_variables
            vector_coords[free_var] = 1
            for index, plane in enumerate(self.planes):
                pivot_var = pivot_indices[index]
                if pivot_var < 0:
                    break
                vector_coords[pivot_var] = -plane.normal_vector[free_var]

            direction_vectors.append(Vector(vector_coords))

        return direction_vectors

    def extract_basepoint_for_parametrization(self):
        num_variables = self.dimension
        pivot_indices = self.indices_of_first_nonzero_terms_in_each_row()

        basepoint_coords = [0] * num_variables

        for index, plane in enumerate(self.planes):
            pivot_var = pivot_indices[index]
            if pivot_var < 0:
                break
            basepoint_coords[pivot_var] = plane.constant_term

        return Vector(basepoint_coords)

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
                if str(e) == Hyperplane.NO_NONZERO_ELTS_FOUND_MSG:
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
            print("WRONG DIMENSIONS (%d) (%d)" % (x.dimension, self.dimension))
            print("XXX\n{}".format(x))
            print("YYY\n{}".format(self))
            raise Exception(self.ALL_PLANES_MUST_BE_IN_SAME_DIM_MSG)


    def __str__(self):
        ret = 'Linear System:\n'
        temp = ['Equation {}: {}'.format(i+1,p) for i,p in enumerate(self.planes)]
        ret += '\n'.join(temp)
        return ret


class MyDecimal(Decimal):
    def is_near_zero(self, eps=1e-10):
        return abs(self) < eps


class Parametrization(object):

    BASEPT_AND_DIR_VECTORS_MUST_BE_IN_SAME_DIM = (
        'The basepoint and direction vectors should all live in the same '
        'dimension')

    def __init__(self, basepoint, direction_vectors):

        self.basepoint = basepoint
        self.direction_vectors = direction_vectors
        self.dimension = self.basepoint.dimension

        try:
            for v in direction_vectors:
                assert v.dimension == self.dimension

        except AssertionError:
            raise Exception(self.BASEPT_AND_DIR_VECTORS_MUST_BE_IN_SAME_DIM)

    def __str__(self):

        output = ''
        for coord in range(self.dimension):
            output += 'x_{} = {} '.format(coord + 1,
                                          round(self.basepoint[coord], 3))
            for free_var, vector in enumerate(self.direction_vectors):
                output += '+ {} t_{}'.format(round(vector[coord], 3),
                                             free_var + 1)
            output += '\n'
        return output


def test_hyperplanes():
    p1 = Hyperplane(normal_vector=Vector([0.786, 0.786]), constant_term=0.786)
    p2 = Hyperplane(normal_vector=Vector([-0.131, -0.131]), constant_term=-0.131)
    
    system = LinearSystem([p1, p2])
    print("Hyperplane #1 solution")
    print(system.compute_solution())
    
    
    p1 = Hyperplane(normal_vector=Vector([2.102, 7.489, -0.786]),
                    constant_term=-5.113)
    p2 = Hyperplane(normal_vector=Vector([-1.131, 8.318, -1.209]),
                    constant_term=-6.775)
    p3 = Hyperplane(normal_vector=Vector([9.015, 5.873, -1.105]),
                    constant_term=-0.831)
    
    system = LinearSystem([p1, p2, p3])
    print("Hyperplane #2 solution")
    print(system.compute_solution())
    
    p1 = Hyperplane(normal_vector=Vector([0.786, 0.786, 8.123, 1.111, -8.363]),
                    constant_term=-9.955)
    p2 = Hyperplane(normal_vector=Vector([0.131, -0.131, 7.05, -2.813, 1.19]),
                    constant_term=-1.991)
    p3 = Hyperplane(normal_vector=Vector([9.015, -5.873, -1.105, 2.013, -2.802]),
                    constant_term=-3.982)
    
    system = LinearSystem([p1, p2, p3])
    print("Hyperplane #3 solution")
    print(system.compute_solution())


def test():
    test_hyperplanes()

if __name__ == '__main__':
    test()
