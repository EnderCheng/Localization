#!/usr/bin/env python3
import math
import time

import gmpy2
import numpy as np

import util
from symmetric import SHE
from line_profiler import LineProfiler
import functions, random

def test_matrix_mul():
    she = SHE(500,160,80)
    bit = 70
    module = pow(2,bit)
    mat_0 = util.generate_random_matrix(3,2,-3000,3000)
    print(mat_0)
    mat_0_inv = mat_0.T
    print(mat_0_inv)


    [a_0, a_1] = util.split_matrix(mat_0_inv,module)
    [b_0, b_1] = util.split_matrix(mat_0,module)

    # lp = LineProfiler()
    # lp_wrapper = lp(functions.matrix_mul)
    # lp_wrapper(a_0, b_0, a_1, b_1, she, bit, module)
    # lp.print_stats()

    # start = time.process_time()
    # for i in range(0,100):
    #     [c_0,c_1] = functions.matrix_mul(a_0,b_0,a_1,b_1,she,bit,module)
    # end = time.process_time()
    # print("matric_mul processing time:"+ str((end-start)/100))

    [c_0, c_1] = functions.matrix_mul(a_0, b_0, a_1, b_1, she, bit, module)
    print(np.matmul(mat_0_inv,mat_0))
    print(util.recover_matrix(c_0,c_1,module))

def test_integer_map():
    she = SHE(500,160,80)
    bit = 70
    module = pow(2,bit)
    test_value = random.SystemRandom().randrange(-10000000,10000000)
    [value_0, value_1] = util.split(test_value,module)
    c_0 = she.encrypt(value_0)
    c_1 = she.encrypt(value_1)

    c = c_0+c_1
    half_module = module // 2
    new_c = functions.integer_mapping(c,she,bit,module, half_module)
    print(test_value)
    print(she.decrypt(new_c))

def test_matrix_inverse():
    she = SHE(500,200,100)
    bit = 95
    module = pow(2,bit)
    scale = 10**29
    value_max = 1500
    half_module = module // 2
    mat = util.generate_random_matrix(2,2,-value_max,value_max)
    [mat_0, mat_1] = util.split_matrix(mat,module)
    [c_0,c_1] = functions.matrix_inverse(mat_0,mat_1,she,bit,module,half_module,1500,scale)
    print(util.recover_matrix(c_0,c_1,module))

    # print(mat)
    mat_sp = util.npmat_to_spmat(mat)
    mat_sp_inv = mat_sp.inv()
    mat_np = util.spmat_to_npmat(mat_sp_inv,float)
    print(mat_np*scale)

def test_compare():
    she = SHE(500,160,80)
    bit = 70
    module = pow(2,bit)
    test_a= random.SystemRandom().randrange(-1000,1000)
    [a_0, a_1] = util.split(test_a,module)
    test_b= random.SystemRandom().randrange(-1000,1000)
    [b_0, b_1] = util.split(test_b,module)
    [c_0,c_1] = functions.compare(a_0,b_0,a_1,b_1,she,bit,'SGET')
    print(c_0^c_1)
    print(test_a < test_b)

def test_divide():
    she = SHE(500, 160, 80)
    bit = 70
    module = pow(2, bit)
    divisor = 1
    test_a = random.SystemRandom().randrange(-1000, 1000)
    [a_0, a_1] = util.split(test_a,module)
    [c_0,c_1] = functions.divide(a_0,a_1,divisor,she,bit,module)
    print(util.recover(c_0,c_1,module))
    print(test_a/divisor)


def test_square():
    she = SHE(500,160,80)
    bit = 70
    module = pow(2,bit)
    test_a= random.SystemRandom().randrange(-10000000,10000000)
    [a_0, a_1] = util.split(test_a,module)
    [c_0,c_1] = functions.square(a_0,a_1,she,bit,module)
    print(util.recover(c_0,c_1,module))
    print(test_a**2)

def test_absolute_difference():
    she_0 = SHE(500,160,80)
    she_1 = SHE(500, 160, 80)
    bit = 70
    module = pow(2,bit)
    test_a= random.SystemRandom().randrange(-1000,1000)
    [a_0, a_1] = util.split(test_a,module)
    test_b= random.SystemRandom().randrange(-1000,1000)
    [b_0, b_1] = util.split(test_b,module)

    # lp = LineProfiler()
    # lp_wrapper = lp(functions.absolute_difference)
    # lp_wrapper(a_0,b_0,a_1,b_1,she_0,she_1,bit,module)
    # lp.print_stats()

    [c_0,c_1] = functions.absolute_difference(a_0,b_0,a_1,b_1,she_0,she_1,bit,module)
    print(util.recover(c_0, c_1, module))
    print(abs(test_a-test_b))

def test_integer_root():
    she = SHE(500,160,80)
    bit = 70
    module = pow(2,bit)
    test_a= random.SystemRandom().randrange(0,10000)
    print(functions.integer_root_square_plain(test_a,bit))
    [a_0, a_1] = util.split(test_a,module)

    # lp = LineProfiler()
    # lp_wrapper = lp(functions.integer_root_square)
    # lp_wrapper(a_0,a_1,she,bit,module)
    # lp.print_stats()

    [c_0,c_1] = functions.integer_root_square(a_0,a_1,she,bit,module)
    print(util.recover(c_0, c_1, module))
    print(math.sqrt(test_a))

def test_oblivious_shuffle():
    she_0 = SHE(500,160,80)
    she_1 = SHE(500,160,80)
    bit = 70
    module = pow(2,bit)
    num = 10

    vec = np.empty(num, dtype=object)
    for i in range(0,num):
        vec[i] = random.SystemRandom().randrange(0,10000)
        print(str(i)+","+str(vec[i]))

    print("-------------")
    seed_0 = 12389
    seed_1 = 45609

    [vec_0,vec_1] = util.split_vector(vec,module)

    lp = LineProfiler()
    lp_wrapper = lp( functions.oblivious_shuffle)
    lp_wrapper(vec_0,vec_1,num,she_0,she_1,bit,module,seed_0,seed_1)
    lp.print_stats()

    [sf_ind_0,sf_vec_0,sf_ind_1,sf_vec_1] = functions.oblivious_shuffle\
        (vec_0,vec_1,num,she_0,she_1,bit,module,seed_0,seed_1)

    for i in range(0,num):
        print(str(util.recover(sf_ind_0[i],sf_ind_1[i],module))+","+str(util.recover(sf_vec_0[i],sf_vec_1[i],module)))

def test_quicksort():
    she_0 = SHE(500,160,80)
    she_1 = SHE(500,160,80)
    bit = 70
    module = pow(2,bit)
    num = 5

    vec = np.empty(num, dtype=object)
    for i in range(0,num):
        vec[i] = random.SystemRandom().randrange(0,10000)
        print(str(i)+","+str(vec[i]))

    print("-------------")
    seed_0 = 12389
    seed_1 = 45609

    [vec_0,vec_1] = util.split_vector(vec,module)
    [sf_ind_0,sf_vec_0,sf_ind_1,sf_vec_1] = functions.oblivious_shuffle\
        (vec_0,vec_1,num,she_0,she_1,bit,module,seed_0,seed_1)

    [theta_0,theta_1] = functions.quicksort(sf_ind_0,sf_vec_0,sf_ind_1,sf_vec_1,num,she_0,bit,module)

    for i in range(0,num):
        print(util.recover(theta_0[i],theta_1[i],module))

def test_distance():
    she = SHE(500,160,80)
    bit = 70
    module = pow(2,bit)
    test_x= random.SystemRandom().randrange(0,3000)
    test_y= random.SystemRandom().randrange(0,3000)

    test_xx= random.SystemRandom().randrange(0,3000)
    test_yy= random.SystemRandom().randrange(0,3000)
    [x_0, x_1] = util.split(test_x, module)
    [xx_0, xx_1] = util.split(test_xx, module)

    [y_0, y_1] = util.split(test_y, module)
    [yy_0, yy_1] = util.split(test_yy, module)

    dist = util.distance([test_x,test_y],[test_xx,test_yy])
    [dist_0,dist_1] = functions.distance(x_0,y_0,xx_0,yy_0,x_1,y_1,xx_1,yy_1,she,bit,module)
    print(dist)
    print(util.recover(dist_0,dist_1,module))

