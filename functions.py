#!/usr/bin/env python3
import copy
import math
import time

import gmpy2
import numpy as np
import random

import functions
import util
from symmetric import SHE
from line_profiler import LineProfiler
from numpy.linalg import inv

def matrix_mul(a_0, b_0, a_1, b_1, she, bit, module):
    if a_0.shape[0] != a_1.shape[0] or a_0.shape[1] != a_1.shape[1] \
            or b_0.shape[0] != b_1.shape[0] or b_0.shape[1] != b_1.shape[1]:
        raise Exception("functions.matrix_mul: matrices are not correct ")

    enc_a_0 = she.encrypt_matrix(a_0)
    enc_b_0 = she.encrypt_matrix(b_0)
    w = (np.matmul(enc_a_0, b_1) + np.matmul(a_1, enc_b_0)) % she.N
    mat_r = util.generate_random_matrix_bit(a_0.shape[0],b_0.shape[1],bit)
    c_1 = (np.matmul(a_1,b_1) - mat_r) %  module
    share_w = w + mat_r
    c_0 = (she.decrypt_matrix(share_w) + np.matmul(a_0,b_0)) % module

    return [c_0, c_1]

def integer_mapping(cipher, she, bit, module, module_half):
    coins = np.empty(2, dtype=object)
    ct =np.empty(2, dtype=object)

    for b in range(0, 2):
        coins[b] = random.randint(0, 1)

        r = np.empty(2, dtype=object)
        r[0] = util.get_number(bit)
        r[1] = util.get_number(bit-1)
        if coins[b] == 1:
            ct[b] = (r[0]*(cipher-(1+2*b)*module_half)+r[1]) % she.N
        else:
            ct[b] = (r[0] * (she.negate(cipher) + (1+2*b)*module_half) - r[1])  % she.N

    w = np.empty(2, dtype=object)
    for b in range(0, 2):
        tmp = she.decrypt(ct[b])
        if tmp > 0:
            w[b] = she.encrypt(-1)
        else:
            w[b] = she.encrypt(0)

    h = np.empty(2, dtype=object)
    for b in range(0, 2):
        if coins[b] == 1:
            h[b] = coins[b] - 1 + w[b]
        else:
            h[b] = coins[b] - 1 + she.negate(w[b])
    # print(she.decrypt(h[0]))
    # print(she.decrypt(h[1]))
    return (cipher + (h[0] + h[1]) * module) % she.N

def matrix_inverse(mat_0, mat_1, she, bit, module, module_half,value_max,scale):
    enc_mat_0 = she.encrypt_matrix(mat_0)
    enc_mat = enc_mat_0 + mat_1
    [rows, cols] = enc_mat.shape
    new_enc_mat = np.zeros(enc_mat.shape,dtype='object')
    for i in range(rows):
        for j in range(cols):
            new_enc_mat[i,j] = integer_mapping(enc_mat[i,j],she, bit, module, module_half)

    mat_r_1 = util.get_random_unimat(rows,value_max)
    mat_r_2 = util.get_random_unimat(rows,value_max)

    # print("randomness")
    # print(mat_r_1)
    # print(mat_r_2)
    # mat_r_1_inv = mat_r_1.inv()
    # mat_r_2_inv = mat_r_2.inv()

    new_enc_mat = np.matmul(np.matmul(mat_r_1,new_enc_mat),mat_r_2)

    for i in range(rows):
        for j in range(cols):
            new_enc_mat[i,j] = she.negate_end(new_enc_mat[i,j])
    dec_mat = she.decrypt_matrix(new_enc_mat)
    dec_mat = util.npmat_to_spmat(dec_mat)
    dec_mat_inv = dec_mat.inv()
    dec_mat = util.spmat_to_npmat(dec_mat_inv*scale,'object')

    # tmp = np.matmul(np.matmul(mat_r_2,dec_mat),mat_r_1)
    # print(int(tmp[0][0]*scale).bit_length())

    enc_mat_ret = she.encrypt_matrix(dec_mat)

    enc_mat_inv = np.matmul(np.matmul(mat_r_2,enc_mat_ret),mat_r_1)

    for i in range(rows):
        for j in range(cols):
            enc_mat_inv[i,j] = she.negate_end(enc_mat_inv[i,j])
    # print(she.decrypt_matrix(enc_mat_inv))
    mat_r = util.generate_random_matrix_bit(enc_mat_inv.shape[0],enc_mat_inv.shape[0],bit)
    enc_mat_inv = enc_mat_inv - mat_r
    c_1 = mat_r % module
    c_0 = (she.decrypt_matrix(enc_mat_inv)) % module

    return [c_0,c_1]

def compare(a_0, b_0, a_1, b_1, she, bit, protocol):
    if protocol == 'SLT':
        d_0 = a_0 - b_0
    else:
        d_0 = b_0 - a_0
    enc_d_0 = she.encrypt(d_0)

    coin = random.randint(0, 1)
    r = np.empty(2, dtype=object)
    r[0] = util.get_number(bit)
    r[1] = util.get_number(bit - 1)
    if protocol == 'SLT':
        d_1 = a_1-b_1
    else:
        d_1 = b_1-a_1
    if coin == 1:
        ct = (r[0] * (enc_d_0+d_1) + r[1]) % she.N
        c_1 = 1
    else:
        ct = (r[0] * (she.negate(enc_d_0) -d_1) - r[1]) % she.N
        c_1 = 0

    dec_d = she.decrypt(ct)
    if dec_d < 0:
        c_0 = 0
    else:
        c_0 = 1
    return [c_0,c_1]

def square(a_0,a_1,she,bit,module):
    enc_a_0 = she.encrypt(a_0)
    r = util.get_number(bit-1)
    c_1 = (a_1**2-r) % module

    enc_ret = (2*a_1*enc_a_0+r) % she.N
    dec_ret = she.decrypt(enc_ret)
    c_0 = (a_0**2 + dec_ret) % module
    return [c_0,c_1]

def divide(a_0,a_1,divisor,she,bit,module):
    module_half = module // 2
    enc_a_0 = she.encrypt(a_0)
    # print(util.recover(a_0,a_1,module))
    r = util.get_number(bit)
    c_0 = -r // divisor
    # share = enc_a_0+a_1 - r
    share = integer_mapping(enc_a_0+a_1,she, bit, module, module_half) + r
    c_1 = (she.decrypt(share) // divisor) % module
    return [c_0,c_1]

def absolute_difference(a_0,b_0,a_1,b_1,she_0,she_1,bit,module):
    [phi_0,phi_1] = compare(a_0,b_0,a_1,b_1,she_0,bit,'SGET')
    enc_phi_0 = she_0.encrypt(phi_0)
    enc_diff = she_0.encrypt(a_0-b_0)

    if phi_1 == 0:
        gamma = enc_phi_0
    else:
        gamma = 1+she_0.negate(enc_phi_0)

    coin = random.randint(0, 1)
    if coin == 0:
        coin = -1

    r_1 = util.get_number(bit)

    if coin == 1:
        w = enc_diff+a_1-b_1+r_1
        h = 2*gamma-1
    else:
        w = she_0.negate(enc_diff) -a_1+b_1-r_1
        h = she_0.negate_end(2*gamma-1)

    enc_ret = she_1.encrypt(-coin*r_1)

    r_0 = util.get_number(bit)
    dec_h = she_0.decrypt(h)
    c_0 = (she_0.decrypt(w)*dec_h-r_0) % module
    if dec_h == 1:
        z = enc_ret + r_0
    else:
        z = she_1.negate(enc_ret)+r_0

    c_1 = she_1.decrypt(z)
    return [c_0,c_1]


def integer_root_square_plain(test,bit):
    w = 0
    expo = math.floor((bit-1)/2-1)
    zeta = pow(2,expo)
    for i in range(0,expo+1):
        # print("test:" + str(test))
        # print("w:"+ str(w))
        c = (2*w+zeta)*pow(2,expo-i)
        # print("c:"+str(c))
        # print(c)
        if test >= c:
            # print(i)
            w = w + zeta
            test = test - c
        zeta = zeta // 2
    return w

def integer_root_square(in_0,in_1,she,bit,module):
    a_0 = copy.deepcopy(in_0)
    a_1 = copy.deepcopy(in_1)
    w_0 = 0
    w_1 = 0
    expo = math.floor((bit-1)/2-1)
    zeta = pow(2,expo)
    for i in range(0,expo+1):
        # print("test:" + str((a_0+a_1) % module))
        enc_w_0 = she.encrypt(w_0)
        enc_w = w_1+enc_w_0
        c = (2*enc_w+zeta)*pow(2,expo-i)
        b_1 = util.get_number(bit) % module
        b_0 = she.decrypt(c-b_1) % module
        [phi_0,phi_1] = functions.compare(a_0,b_0,a_1,b_1,she,bit,'SGET')
        enc_phi_0 = she.encrypt(phi_0)
        if phi_1 == 0:
            rho = enc_phi_0
        else:
            rho = 1+she.negate(enc_phi_0)

        # if she.decrypt(rho) == 1:
        #     print(she.decrypt(c))
        #     print(i)
        rr_1 = util.get_number(bit//2)
        w_1 = w_1 - rr_1
        z_1 = util.get_number(bit) % module
        r_1 = util.get_number(bit)
        chi = (c-b_1)*z_1+(rho-z_1)*b_1-r_1
        gamma_1 = (z_1*b_1+r_1) % module
        a_1 = (a_1 - gamma_1) % module
        z_0 = she.decrypt(rho-z_1) % module
        gamma_0 = (she.decrypt(chi)+z_0*b_0) % module
        # print("c:"+str(she.decrypt(c)))
        # print("beta:"+str((b_0+b_1) % module))
        # print("gamma:"+str((gamma_0+gamma_1) % module))
        rr_0 = she.decrypt(rho*zeta+rr_1)
        w_0 = w_0 + rr_0
        # print("w:"+str((w_0+w_1) % module))
        a_0 = (a_0-gamma_0) % module
        zeta = zeta // 2
    return [w_0,w_1]

def oblivious_shuffle(vec_0,vec_1, num, she_0,she_1,bit,module,seed_0,seed_1):
    enc_vec = np.empty(num, dtype=object)
    for i in range(0,num):
        # print("test:" + str((a_0+a_1) % module))
        enc_vec[i] = she_0.encrypt(vec_0[i])

    enc_ind = np.empty(num, dtype=object)
    enc_r_0 = np.empty(num, dtype=object)
    enc_r_1 = np.empty(num, dtype=object)
    for i in range(0,num):
        r_0 = util.get_number(bit)
        r_1 = util.get_number(bit)
        enc_ind[i] = she_0.pub_encrypt(i+r_0)
        enc_vec[i] = enc_vec[i] + vec_1[i] + r_1
        enc_r_0[i] = she_1.encrypt(r_0)
        enc_r_1[i] = she_1.encrypt(r_1)

    # random.seed(seed_0)

    random.Random(seed_0).shuffle(enc_ind)
    random.Random(seed_0).shuffle(enc_vec)
    random.Random(seed_0).shuffle(enc_r_0)
    random.Random(seed_0).shuffle(enc_r_1)

    # random.seed(seed_1)

    random.Random(seed_1).shuffle(enc_ind)
    random.Random(seed_1).shuffle(enc_vec)
    random.Random(seed_1).shuffle(enc_r_0)
    random.Random(seed_1).shuffle(enc_r_1)

    for i in range(0,num):
        r_0 = util.get_number(bit)
        r_1 = util.get_number(bit)
        enc_ind[i] = enc_ind[i] + r_0
        enc_vec[i] = enc_vec[i] + r_1
        enc_r_0[i] = enc_r_0[i] + r_0
        enc_r_1[i] = enc_r_1[i] + r_1

    sf_ind_1 = np.empty(num, dtype=object)
    sf_vec_1 = np.empty(num, dtype=object)
    for i in range(0,num):
        sf_ind_1[i] = util.get_number(bit)
        sf_vec_1[i] = util.get_number(bit)
        dec_r_0 = she_1.decrypt(enc_r_0[i])
        dec_r_1 = she_1.decrypt(enc_r_1[i])
        enc_ind[i] = enc_ind[i] - dec_r_0 - sf_ind_1[i]
        enc_vec[i] = enc_vec[i] - dec_r_1 - sf_vec_1[i]

    sf_ind_0 = np.empty(num, dtype=object)
    sf_vec_0 = np.empty(num, dtype=object)
    for i in range(0,num):
        sf_ind_0[i] = she_0.decrypt(enc_ind[i]) % module
        sf_vec_0[i] = she_0.decrypt(enc_vec[i]) % module

    return [sf_ind_0,sf_vec_0,sf_ind_1,sf_vec_1]


def quicksort(ind_0,vec_0,ind_1,vec_1,num,she,bit,module):
    enc_ind = np.empty(num, dtype=object)
    for i in range(0,num):
        # print("test:" + str((a_0+a_1) % module))
        enc_ind[i] = she.encrypt(ind_0[i])
    for i in range(0,num):
        enc_ind[i] = enc_ind[i] + ind_1[i]

    [enc_ind, vec_0, vec_1] = sort(enc_ind,vec_0,vec_1,0,num-1,she,bit)

    theta_1 = np.empty(num, dtype=object)
    enc_theta = np.empty(num, dtype=object)
    for i in range(0, num):
        theta_1[i] = util.get_number(bit)
        enc_theta[i] = enc_ind[i] - theta_1[i]

    theta_0 = np.empty(num, dtype=object)
    for i in range(0, num):
        theta_0[i] = she.decrypt(enc_theta[i]) % module

    return [theta_0,theta_1]

def sort(enc_ind, vec_0, vec_1,left,right,she,bit):
    if left < right:
        [mid, enc_ind, vec_0, vec_1] = partition(enc_ind,vec_0,vec_1,left,right,she,bit)
        [enc_ind, vec_0, vec_1] = sort(enc_ind,vec_0,vec_1,left,mid-1,she,bit)
        [enc_ind, vec_0, vec_1] = sort(enc_ind, vec_0, vec_1, mid+1, right,she,bit)
    return [enc_ind, vec_0, vec_1]

def partition(enc_ind, vec_0, vec_1,begin,end,she,bit):
    p_value_1 = vec_1[end]
    p_value_0 = vec_0[end]
    i = begin-1
    for j in range(begin,end):
        [c_0, c_1] = functions.compare(vec_0[j], p_value_0, vec_1[j], p_value_1, she, bit, 'SLT')
        rho = c_0^c_1
        if rho == 1:
            i = i + 1
            enc_ind[i], enc_ind[j] = enc_ind[j], enc_ind[i]
            vec_0[i], vec_0[j] = vec_0[j], vec_0[i]
            vec_1[i], vec_1[j] = vec_1[j], vec_1[i]

    enc_ind[i+1], enc_ind[end] = enc_ind[end], enc_ind[i+1]
    vec_0[i+1], vec_0[end] = vec_0[end], vec_0[i+1]
    vec_1[i+1], vec_1[end] = vec_1[end], vec_1[i+1]
    mid = i+1
    return [mid, enc_ind, vec_0, vec_1]

def distance(x_0,y_0,xx_0,yy_0,x_1,y_1,xx_1,yy_1,she,bit,module):
    diff_x_0 = (x_0-xx_0) % module
    diff_y_0 = (y_0-yy_0) % module

    diff_x_1 = (x_1-xx_1) % module
    diff_y_1 = (y_1-yy_1) % module

    [dist_x_sq_0, dist_x_sq_1] = square(diff_x_0,diff_x_1,she,bit,module)
    [dist_y_sq_0, dist_y_sq_1] = square(diff_y_0, diff_y_1, she, bit, module)

    in_0 = (dist_x_sq_0+dist_y_sq_0) % module
    in_1 = (dist_x_sq_1 + dist_y_sq_1) % module
    [dist_0, dist_1] = integer_root_square(in_0,in_1,she,bit,module)

    return [dist_0,dist_1]



















