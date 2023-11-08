#!/usr/bin/env python3
from fractions import Fraction

import gmpy2
import sympy

import test
from symmetric import SHE
from util import split,recover,get_number_raw,get_number,get_random_unimat
from numpy.linalg import inv,pinv
from line_profiler import LineProfiler
import numpy as np
import functions,time
from phe import paillier


if __name__ == '__main__':
    test.test_divide()
    test.test_distance()
    test.test_matrix_mul()
    test.test_integer_map()
    test.test_matrix_inverse()
    test.test_compare()
    test.test_square()
    test.test_absolute_difference()
    test.test_integer_root()
    test.test_oblivious_shuffle()
    test.test_quicksort()



    



