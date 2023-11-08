#!/usr/bin/env python3

import random, gmpy2
from fractions import Fraction

import numpy as np
import sympy as sp
from sympy.combinatorics import Permutation
from sympy import Matrix, zeros, randMatrix
import math
from itertools import combinations

import functions


def spmat_to_npmat(mat, type):
    return np.array(mat).astype(type)

def npmat_to_spmat(mat):
    return Matrix(mat)

def get_random_mat(n,value_max):
    randfunc = random.SystemRandom()
    while True:
        mat_r = np.empty([n,n],dtype='float')
        for i in range(n):
            for j in range(n):
                mat_r[i,j] = randfunc.randrange(-value_max, value_max)
        if np.linalg.det(mat_r) != 0:
            break
    return mat_r

def get_random_unimat(n,value_max):
    randfunc = random.SystemRandom()

    p_1 = sp.eye(n)
    p_2 = sp.eye(n)
    # print(npmat_to_spmat(spmat_to_npmat(p_1,'int')))

    r = sp.ones(n,n)
    for i in range(n):
        for j in range(n):
            if i!=j:
                r[i,j] = randfunc.randrange(-value_max, value_max)

    r_1 = r.lower_triangular()
    r_2 = r.upper_triangular()
    # print(r_1)
    # print(r_2)
    return p_1*r_1*r_2*p_2

# def get_random_unimat(n,value_max):
#     randfunc = random.SystemRandom()
#
#     p_1 = np.identity(n, dtype='int')
#     np.random.shuffle(p_1)
#     p_2 = np.identity(n, dtype='int')
#     np.random.shuffle(p_2)
#
#     r = np.ones([n,n], dtype='int')
#     for i in range(n):
#         for j in range(n):
#             if i!=j:
#                 r[i,j] = randfunc.randrange(0, value_max)
#
#     r_1 = np.tril(r)
#     r_2 = np.triu(r)
#     print(np.matrix(r_1))
#     print(np.matrix(r_2))
#     return np.matmul(np.matmul(np.matmul(p_1,r_1),r_2),p_2)




def get_number(n):
    randfunc = random.SystemRandom()
    r = gmpy2.mpz(randfunc.getrandbits(n))
    r = gmpy2.bit_set(r, n - 1)
    return int(r)

def generate_matrix_pps(row,col,num,min,max):
    matrix_zero = zeros(row,col)
    matrices = []
    matrix_all = zeros(row,col)
    for i in range(0,num-1):
        matrix_r = randMatrix(row,col,min,max)
        matrix_all = matrix_all + matrix_r
        matrices.append(matrix_r)
    matrices.append(matrix_zero-matrix_all)
    return matrices

def get_number_raw(n):
    randfunc = random.SystemRandom()
    r = randfunc.getrandbits(n)
    return r

def split(value, module):
    half_module = module // 2
    abs_value = abs(value)
    if abs_value > half_module:
        raise Exception("util.split: plaintext overflow")

    if value < 0:
        value = module - abs_value

    s_1 = int(random.SystemRandom().randrange(1, module))
    s_2 = (value - s_1) % module

    return [s_1,s_2]


def recover(s_1,s_2,module):
    half_module = module // 2
    value = (s_1 + s_2) % module
    if value > half_module:
        value = value - module
    return value

def split_vector(vec, module):
    num = len(vec)
    vec_0 = np.empty(num, dtype=object)
    vec_1 = np.empty(num, dtype=object)
    for i in range(num):
        s_1,s_2 = split(vec[i],module)
        vec_0[i] = s_1
        vec_1[i] = s_2

    return [vec_0, vec_1]

def split_coordinate(cds, module):
    num = len(cds)
    coord_0 = np.empty(num, dtype=object)
    coord_1 = np.empty(num, dtype=object)
    for i in range(num):
        s_1,s_2 = split(cds[i][0],module)
        s_3,s_4 = split(cds[i][1],module)
        coord_0[i] = [s_1,s_3]
        coord_1[i] = [s_2,s_4]

    return [coord_0,coord_1]


def split_matrix(mat, module):
    if not isinstance(mat, np.ndarray):
        raise Exception("util.split_matrix: input is not a matrix")

    [rows,cols] =  mat.shape
    a = np.zeros(mat.shape,dtype='object')
    b = np.zeros(mat.shape,dtype='object')
    for i in range(rows):
        for j in range(cols):
            s_1,s_2 = split(mat[i,j],module)
            a[i,j] = s_1
            b[i,j] = s_2

    return [a, b]

def recover_matrix(a, b, module):
    if not isinstance(a, np.ndarray) or not isinstance(b, np.ndarray) :
        raise Exception("util.split_matrix: inputs are not a matrix")

    if a.shape[0] != b.shape[0] or a.shape[1] != b.shape[1]:
        raise Exception("util.recover_matrix: inputs have different matrix shape")

    [rows,cols] =  a.shape
    mat = np.zeros(a.shape, dtype='object')
    for i in range(rows):
        for j in range(cols):
            s = recover(a[i,j],b[i,j],module)
            mat[i,j] = s
    return mat

def generate_random_matrix(rows, cols, value_min, value_max):
    mat = np.zeros([rows, cols],dtype='object')
    randfunc = random.SystemRandom()
    for i in range(rows):
        for j in range(cols):
            mat[i,j] = int(randfunc.randrange(value_min, value_max))
    return mat

def generate_random_matrix_bit(rows, cols, bit):
    mat = np.zeros([rows, cols],dtype='object')
    randfunc = random.SystemRandom()
    for i in range(rows):
        for j in range(cols):
            mat[i, j] = int(get_number(bit))
            # mat[i,j] = int(randfunc.getrandbits(bit))
    return mat

def subset_all(num, min):
    all_lists = []
    for r in range(min,num):
        arr = range(1,num)
        lists =  list(combinations(arr, r))
        for i in range(0,len(lists)):
            tmp = list(lists[i])
            all_lists.append(tmp)
    # print(all_lists)
    return all_lists


def distance(p1,p2):
    dist = math.sqrt((p1[0] - p2[0]) ** 2 + (p1[1] - p2[1]) ** 2)
    return dist

def generate_toa_measure(ue,acs,num, los_err, nlos_err):
    radius = np.zeros(num)
    for i in range(0,num):
        radius[i] = distance(ue,acs[i]) + los_err
        if i == num-1:
            radius[i] =  radius[i] + nlos_err
    return radius

def toa_localize(acs,toas,num):
    matrix_a = zeros(num-1,2)
    for i in range(0,num-1):
        for j in range(0,2):
            matrix_a[i,j] = 2*(acs[i+1][j] - acs[0][j])

    vector_b = zeros(num-1,1)
    for i in range(0,num-1):
        vector_b[i,0] = acs[i+1][0]**2+acs[i+1][1]**2-acs[0][0]**2-acs[0][1]**2-toas[i+1]**2+toas[0]**2

    matrix_a_trans = matrix_a.transpose()
    z = (matrix_a_trans*matrix_a).inv()*matrix_a_trans*vector_b
    # print(z)
    return z

def secure_toa_localize(acs, toas, num, she, bit, module):
    acs_sq = [[i[0] ** 2, i[1] ** 2] for i in acs]
    toas_sq =[round(i) ** 2 for i in toas]
    [acs_0,acs_1] = split_coordinate(acs,module)
    [acs_0_sq,acs_1_sq] = split_coordinate(acs_sq,module)
    [toas_0_sq,toas_1_sq] = split_vector(toas_sq, module)

    matrix_a_0 = np.zeros([num-1,2],dtype='object')
    matrix_a_1 = np.zeros([num-1,2],dtype='object')
    for i in range(0,num-1):
        for j in range(0,2):
            matrix_a_0[i,j] = 2*(acs_0[i+1][j] - acs_0[0][j]) % module
            matrix_a_1[i,j] = 2*(acs_1[i+1][j] - acs_1[0][j]) % module


    vector_b_0 = np.zeros([num-1,1],dtype='object')
    vector_b_1 = np.zeros([num-1,1],dtype='object')
    for i in range(0,num-1):
        vector_b_0[i,0] = (acs_0_sq[i+1][0]+acs_0_sq[i+1][1]-acs_0_sq[0][0]-acs_0_sq[0][1]-toas_0_sq[i+1]+toas_0_sq[0]) % module
        vector_b_1[i,0] = (acs_1_sq[i+1][0]+acs_1_sq[i+1][1]-acs_1_sq[0][0]-acs_1_sq[0][1]-toas_1_sq[i+1]+toas_1_sq[0]) % module

    matrix_a_trans_0 = matrix_a_0.T
    matrix_a_trans_1 = matrix_a_1.T
    scale = 10**20
    half_module = module // 2
    [tmp_0,tmp_1] = functions.matrix_mul(matrix_a_trans_0,matrix_a_0,matrix_a_trans_1,matrix_a_1,she,bit,module)
    [tmp_0_inv,tmp_1_inv] = functions.matrix_inverse(tmp_0,tmp_1,she,bit,module,half_module,500,scale)
    [tmp_2,tmp_3] = functions.matrix_mul(matrix_a_trans_0,vector_b_0,matrix_a_trans_1,vector_b_1, she, bit, module)
    [z_0,z_1] = functions.matrix_mul(tmp_0_inv,tmp_2,tmp_1_inv,tmp_3,she,bit,module)
    # print(recover_matrix(z_0,z_1,module))
    return [z_0,z_1]

def toa_identify(acs,toas,num):
    ref_loc = list(toa_localize(acs,toas,num))
    # ref_loc = [1500,1500]
    # print(ref_loc)
    all_subsets = subset_all(num,2)
    set_n = len(all_subsets)
    residuals = np.zeros(set_n, dtype=object)
    diff_ac = np.zeros(num-1, dtype=object)

    # for i in range(0,num):
    #     print("real4:"+"("+str(i)+"):"+str(toas[i]))

    for j in range(1,num):
        dist = distance(acs[j], ref_loc)
        # print("real:"+"("+str(j)+"):"+str(dist))
        diff_ac[j-1] = (toas[j] - dist) ** 2
        # print("real3:" + "(" + str(j) + "):" + str((toas[j] - dist)))
        # print("real2:" + "(" + str(j) + "):" + str(diff_ac[j-1]))
        # print(diff_ac[j-1])

    for i in range(0,set_n):
        residuals[i] = 0
        for j in all_subsets[i]:
            residuals[i] = residuals[i] +  diff_ac[j-1]
        # print(residuals[i])
    # print(residuals)
    acs_res = np.zeros(num, dtype=float)
    # print(acs_res)
    set_num = 0
    for subset in all_subsets:
        for i in range(0,num):
            if i in subset:
                acs_res[i] = acs_res[i] + residuals[set_num]
        set_num = set_num + 1
    # for i in range(0, num):
    #     print(acs_res[i])
    return acs_res

def secure_toa_identify(acs,toas,num,she_0,she_1,bit,module):
    [z_0,z_1] = list(secure_toa_localize(acs,toas,num,she_0,bit,module))
    scale = 10**20
    [x_0,x_1] = functions.divide(z_0[0],z_1[0],scale,she_0,bit,module)
    [y_0,y_1] = functions.divide(z_0[1],z_1[1],scale,she_0,bit,module)

    all_subsets = subset_all(num,2)
    set_n = len(all_subsets)
    residuals_0 = np.zeros(set_n, dtype=object)
    residuals_1 = np.zeros(set_n, dtype=object)

    diff_ac_0 = np.zeros(num-1, dtype=object)
    diff_ac_1 = np.zeros(num-1, dtype=object)

    toas =[round(i) for i in toas]
    [acs_0,acs_1] = split_coordinate(acs,module)
    [toas_0,toas_1] = split_vector(toas, module)
    for i in range(1, num):
        [xx_0,xx_1] = [acs_0[i][0], acs_1[i][0]]
        [yy_0,yy_1] = [acs_0[i][1], acs_1[i][1]]
        [dist_0, dist_1] = functions.distance(x_0, y_0, xx_0, yy_0, x_1, y_1, xx_1, yy_1, she_0, bit, module)
        diff_0 = (toas_0[i] - dist_0) % module
        diff_1 = (toas_1[i] - dist_1) % module
        [diff_ac_tmp_0,diff_ac_tmp_1] = functions.square(diff_0,diff_1,she_0,bit,module)
        diff_ac_0[i - 1] = diff_ac_tmp_0
        diff_ac_1[i - 1] = diff_ac_tmp_1
        # print(recover(diff_ac_0[i-1],diff_ac_1[i-1],module))

    for i in range(0,set_n):
        residuals_0[i] = 0
        residuals_1[i] = 0
        for j in all_subsets[i]:
            residuals_0[i] = (residuals_0[i] +  diff_ac_0[j-1]) % module
            residuals_1[i] = (residuals_1[i] +  diff_ac_1[j-1]) % module
        # print(recover(residuals_0[i], residuals_1[i], module))
    acs_res_0 = np.zeros(num, dtype=object)
    acs_res_1 = np.zeros(num, dtype=object)
    # print(acs_res)
    set_num = 0
    for subset in all_subsets:
        for i in range(0,num):
            if i in subset:
                acs_res_0[i] = (acs_res_0[i] + residuals_0[set_num]) % module
                acs_res_1[i] = (acs_res_1[i] + residuals_1[set_num]) % module
        set_num = set_num + 1

    seed_0 = 12389
    seed_1 = 45609
    [sf_ind_0,sf_vec_0,sf_ind_1,sf_vec_1] = functions.oblivious_shuffle(acs_res_0,acs_res_1,num,she_0,she_1,bit,module,seed_0,seed_1)
    [theta_0, theta_1] = functions.quicksort(sf_ind_0, sf_vec_0, sf_ind_1, sf_vec_1, num, she_0, bit, module)
    # for i in range(0, num):
    #     print(recover(theta_0[i], theta_1[i], module))
    # return [theta_0,theta_1]
    acs_res = []
    for i in range(0, num):
        acs_res.append(recover(acs_res_0[i],acs_res_1[i],module))
    return acs_res

def generate_tdoa_measure(ue,acs,num, los_err, nlos_err):
    radius = np.zeros(num)
    for i in range(0,num):
        radius[i] = distance(ue,acs[i])
    radius_diff = np.zeros(num-1)
    for i in range(0, num-1):
        radius_diff[i] = abs(radius[i+1]-radius[0]) + los_err
        if i == num-2:
            radius_diff[i] = radius_diff[i] + nlos_err
    return radius_diff


def tdoa_localize(acs,tdoas,num):
    matrix_a = zeros(num-1,3)
    for i in range(0,num-1):
        for j in range(0,3):
            if j == 0 or j == 1:
                matrix_a[i,j] = acs[0][j] - acs[i+1][j]
            else:
                matrix_a[i,j] = tdoas[i]

    vector_b = zeros(num-1,1)
    for i in range(0,num-1):
        vector_b[i,0] = tdoas[i]**2 - acs[i+1][0]**2 - acs[i+1][1]**2 + acs[0][0]**2 + acs[0][1]**2

    matrix_a_trans = matrix_a.transpose()
    z = (matrix_a_trans*matrix_a).inv()*matrix_a_trans*vector_b
    # print(matrix_a_trans*matrix_a)
    # print(matrix_a_trans*vector_b)
    # print(z/2)
    return z/2

def secure_tdoa_localize(acs,tdoas,num, she, bit, module):
    acs_sq = [[i[0] ** 2, i[1] ** 2] for i in acs]
    tdoas =[round(i) for i in tdoas]
    tdoas_sq =[round(i) ** 2 for i in tdoas]
    [acs_0,acs_1] = split_coordinate(acs,module)
    [tdoas_0,tdoas_1] = split_vector(tdoas, module)
    [acs_0_sq,acs_1_sq] = split_coordinate(acs_sq,module)
    [tdoas_0_sq,tdoas_1_sq] = split_vector(tdoas_sq, module)

    matrix_a_0 = np.zeros([num-1,3], dtype='object')
    matrix_a_1 = np.zeros([num-1,3], dtype='object')
    for i in range(0,num-1):
        for j in range(0,3):
            if j == 0 or j == 1:
                matrix_a_0[i,j] = (acs_0[0][j] - acs_0[i+1][j]) % module
                matrix_a_1[i,j] = (acs_1[0][j] - acs_1[i+1][j]) % module
            else:
                matrix_a_0[i,j] = tdoas_0[i] % module
                matrix_a_1[i,j] = tdoas_1[i] % module

    vector_b_0 = np.zeros([num-1,1], dtype='object')
    vector_b_1 = np.zeros([num-1,1], dtype='object')
    for i in range(0,num-1):
        vector_b_0[i,0] = (tdoas_0_sq[i] - acs_0_sq[i+1][0] - acs_0_sq[i+1][1] + acs_0_sq[0][0] + acs_0_sq[0][1]) % module
        vector_b_1[i,0] = (tdoas_1_sq[i] - acs_1_sq[i+1][0] - acs_1_sq[i+1][1] + acs_1_sq[0][0] + acs_1_sq[0][1]) % module

    matrix_a_trans_0 = matrix_a_0.T
    matrix_a_trans_1 = matrix_a_1.T
    scale = 10**20
    half_module = module // 2
    [tmp_0,tmp_1] = functions.matrix_mul(matrix_a_trans_0,matrix_a_0,matrix_a_trans_1,matrix_a_1,she,bit,module)
    [tmp_0_inv,tmp_1_inv] = functions.matrix_inverse(tmp_0,tmp_1,she,bit,module,half_module,500,scale)
    [tmp_2,tmp_3] = functions.matrix_mul(matrix_a_trans_0,vector_b_0,matrix_a_trans_1,vector_b_1, she, bit, module)
    [z_0,z_1] = functions.matrix_mul(tmp_0_inv,tmp_2,tmp_1_inv,tmp_3,she,bit,module)
    # print(recover_matrix(z_0,z_1,module))
    return [z_0,z_1]

def tdoa_identify(acs,tdoas,num):
    ref_loc = list(tdoa_localize(acs,tdoas,num))
    ref_loc= ref_loc[0:2]
    all_subsets = subset_all(num,3)
    set_n = len(all_subsets)
    residuals = np.zeros(set_n, dtype=object)
    diff_ac = np.zeros(num-1, dtype=object)

    dist_base = distance(acs[0], ref_loc)
    for j in range(1,num):
        dist = abs(distance(acs[j], ref_loc)-dist_base)
        diff_ac[j-1] = (tdoas[j-1] - dist) ** 2

    for i in range(0,set_n):
        residuals[i] = 0
        for j in all_subsets[i]:
            residuals[i] = residuals[i] +  diff_ac[j-1]
    acs_res = np.zeros(num, dtype=float)
    set_num = 0
    for subset in all_subsets:
        for i in range(0,num):
            if i in subset:
                acs_res[i] = acs_res[i] + residuals[set_num]
        set_num = set_num + 1
    # for i in range(0,num):
    #     print(acs_res[i])
    return acs_res

def secure_tdoa_identify(acs,tdoas,num,she_0,she_1,bit,module):
    [z_0,z_1] = list(secure_tdoa_localize(acs,tdoas,num,she_0,bit,module))
    scale = 2*10**20
    [x_0,x_1] = functions.divide(z_0[0],z_1[0],scale,she_0,bit,module)
    [y_0,y_1] = functions.divide(z_0[1],z_1[1],scale,she_0,bit,module)

    all_subsets = subset_all(num,3)
    set_n = len(all_subsets)
    residuals_0 = np.zeros(set_n, dtype=object)
    residuals_1 = np.zeros(set_n, dtype=object)

    diff_ac_0 = np.zeros(num-1, dtype=object)
    diff_ac_1 = np.zeros(num-1, dtype=object)

    tdoas =[round(i) for i in tdoas]
    [acs_0,acs_1] = split_coordinate(acs,module)
    [tdoas_0,tdoas_1] = split_vector(tdoas, module)
    [dist_base_0,dist_base_1] = functions.distance(x_0, y_0, acs_0[0][0], acs_0[0][1], x_1, y_1, acs_1[0][0], acs_1[0][1], she_0, bit, module)
    for i in range(1, num):
        [xx_0,xx_1] = [acs_0[i][0], acs_1[i][0]]
        [yy_0,yy_1] = [acs_0[i][1], acs_1[i][1]]
        [dist_00, dist_11] = functions.distance(x_0, y_0, xx_0, yy_0, x_1, y_1, xx_1, yy_1, she_0, bit, module)
        [dist_0,dist_1] = functions.absolute_difference(dist_00,dist_base_0,dist_11,dist_base_1,she_0,she_1,bit,module)
        diff_0 = (tdoas_0[i-1] - dist_0) % module
        diff_1 = (tdoas_1[i-1] - dist_1) % module
        [diff_ac_tmp_0,diff_ac_tmp_1] = functions.square(diff_0,diff_1,she_0,bit,module)
        diff_ac_0[i - 1] = diff_ac_tmp_0
        diff_ac_1[i - 1] = diff_ac_tmp_1

    for i in range(0,set_n):
        residuals_0[i] = 0
        residuals_1[i] = 0
        for j in all_subsets[i]:
            residuals_0[i] = (residuals_0[i] +  diff_ac_0[j-1]) % module
            residuals_1[i] = (residuals_1[i] +  diff_ac_1[j-1]) % module
        # print(recover(residuals_0[i], residuals_1[i], module))
    acs_res_0 = np.zeros(num, dtype=object)
    acs_res_1 = np.zeros(num, dtype=object)
    # print(acs_res)
    set_num = 0
    for subset in all_subsets:
        for i in range(0,num):
            if i in subset:
                acs_res_0[i] = (acs_res_0[i] + residuals_0[set_num]) % module
                acs_res_1[i] = (acs_res_1[i] + residuals_1[set_num]) % module
        set_num = set_num + 1
    # for i in range(0, num):
    #     print(recover(acs_res_0[i],acs_res_1[i],module))
    seed_0 = 12389
    seed_1 = 45609
    [sf_ind_0,sf_vec_0,sf_ind_1,sf_vec_1] = functions.oblivious_shuffle(acs_res_0,acs_res_1,num,she_0,she_1,bit,module,seed_0,seed_1)
    [theta_0, theta_1] = functions.quicksort(sf_ind_0, sf_vec_0, sf_ind_1, sf_vec_1, num, she_0, bit, module)
    # for i in range(0, num):
    #     print(recover(theta_0[i], theta_1[i], module))
    ret = []
    for i in range(0, num):
        ret.append(recover(acs_res_0[i],acs_res_1[i],module))
    # return [theta_0,theta_1]
    return ret

def outside_circle_to_vectors(circle_x, circle_y, radius, num_pts):
    b = [0] * num_pts
    A = [0] * num_pts
    for i in range(num_pts):
        A[i] = [0] * 2

    for i in range(num_pts):
        x = round(circle_x + radius * math.cos(2 * math.pi * i / num_pts),5)
        y = round(circle_y + radius * math.sin(2 * math.pi * i / num_pts),5)
        [w,v,z] = tangent_line(circle_x,circle_y,x,y)
        A[i][0] = w
        A[i][1] = v
        b[i] = z
    return [A, b]

def tangent_line(circle_x, circle_y, x, y):
    if x == circle_x:
        if y >= circle_y:
            return [0, 1, y]
        else:
            return [0, -1, -y]
    elif y == circle_y:
        if x >=circle_x:
            return [1, 0, x]
        else:
            return [-1, 0, -x]
    else:
        k = -((x - circle_x)+0.0)/(y-circle_y)
    b = y-k*x # y = kx+b
    if x< circle_x and y > circle_y:
        return [-k, 1, b]
    if x> circle_x and y > circle_y:
        return [-k, 1, b]
    if x< circle_x and y < circle_y:
        return [k, -1, -b]
    if x> circle_x and y < circle_y:
        return [k, -1, -b]
