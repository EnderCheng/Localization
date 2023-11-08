#!/usr/bin/env python3
import math
import random

import numpy as np
from line_profiler import LineProfiler

import util
from pyfhe.ckks.ckks_decryptor import CKKSDecryptor
from pyfhe.ckks.ckks_encoder import CKKSEncoder
from pyfhe.ckks.ckks_encryptor import CKKSEncryptor
from pyfhe.ckks.ckks_evaluator import CKKSEvaluator
from pyfhe.ckks.ckks_key_generator import CKKSKeyGenerator
from pyfhe.ckks.ckks_parameters import CKKSParameters
import time

# def encode(value, num_slot, encoder,scaling_factor):
#     message = [0] * num_slot
#     message[0] = value
#     return encoder.encode(message,scaling_factor)
#
# def decode(value, encoder):
#     return encoder.decode(value)

def flatten(l):
    return [item for sublist in l for item in sublist]

def secure_abolute(enc_c, iterations, num_slot,encoder,encryptor,decryptor,evaluator,relin_key):
    enc_c_sq = evaluator.multiply(enc_c,enc_c,relin_key)
    enc_c_sq = evaluator.rescale(enc_c_sq,evaluator.scaling_factor)
    return ckks_sqrt_abs(enc_c_sq,iterations,num_slot,encoder,encryptor,decryptor,evaluator,relin_key)

def secure_distance(enc_p1, enc_p2, iterations, num_slot,encoder,encryptor,decryptor,evaluator,relin_key):
    [p1_x,p1_y] = enc_p1
    [p2_x,p2_y] = enc_p2

    left = evaluator.subtract(p1_x, p2_x)
    right = evaluator.subtract(p1_y, p2_y)

    left_sq = evaluator.multiply(left,left,relin_key)
    left_sq = evaluator.rescale(left_sq,evaluator.scaling_factor)
    right_sq = evaluator.multiply(right,right,relin_key)
    right_sq = evaluator.rescale(right_sq,evaluator.scaling_factor)

    tmp =  evaluator.add(left_sq,right_sq)
    return  ckks_sqrt(tmp, iterations, num_slot,encoder,encryptor,decryptor,evaluator,relin_key)

def sqrt(value, iterations):
    x = 0.5*(1.0/math.sqrt(3000**2+3000**2))
    for i in range(0,iterations):
        x = 1.5*x - 0.5*value*(x**3)
        # print(x)
    return value*x

def ckks_sqrt_abs(enc_c, iterations, num_slot,encoder,encryptor,decryptor,evaluator,relin_key):
    x = 0.5 * (1.0 / math.sqrt(3000 ** 2 + 3000 ** 2))
    message = [0] * num_slot
    message[0] = x
    enc_x = encryptor.encrypt(encoder.encode(message, evaluator.scaling_factor))
    plain_half = evaluator.create_constant_plain(0.5)
    value_half = evaluator.multiply_plain(enc_c, plain_half)
    value_half = evaluator.rescale(value_half, evaluator.scaling_factor)
    # value_half = evaluator.lower_modulus(value_half, evaluator.scaling_factor)
    for i in range(0, iterations):
        enc_x_double = evaluator.multiply(enc_x, enc_x, relin_key)
        enc_x_double = evaluator.rescale(enc_x_double, evaluator.scaling_factor)
        enc_x = evaluator.lower_modulus(enc_x, evaluator.scaling_factor)
        enc_x_triple = evaluator.multiply(enc_x_double, enc_x, relin_key)
        enc_x_triple = evaluator.rescale(enc_x_triple, evaluator.scaling_factor)
        enc_x_triple = evaluator.lower_modulus(enc_x_triple, evaluator.scaling_factor ** 5)
        tmp_2 = evaluator.multiply(value_half, enc_x_triple, relin_key)
        tmp_2 = evaluator.rescale(tmp_2, evaluator.scaling_factor)
        plain_onehalf = evaluator.create_constant_plain(1.5)
        tmp_1 = evaluator.multiply_plain(enc_x, plain_onehalf)
        tmp_1 = evaluator.rescale(tmp_1, evaluator.scaling_factor)
        tmp_1 = evaluator.lower_modulus(tmp_1, evaluator.scaling_factor ** 6)
        enc_x = evaluator.subtract(tmp_1, tmp_2)
        # value_half = evaluator.lower_modulus(value_half, evaluator.scaling_factor)
        # value_half = evaluator.lower_modulus(value_half, evaluator.scaling_factor)
        # print(math.log(value_half.modulus,2))
        r = random.random()
        plain_r = evaluator.create_constant_plain(-r)
        enc_x_r = evaluator.add_plain(enc_x, plain_r)
        x_r = encoder.decode(decryptor.decrypt(enc_x_r))[0].real
        message[0] = x_r
        re_enc_x_r = encryptor.encrypt(encoder.encode(message, evaluator.scaling_factor))
        message[0] = r
        plain_r = encoder.encode(message, evaluator.scaling_factor)
        enc_x = evaluator.add_plain(re_enc_x_r, plain_r)
        # print(encoder.decode(decryptor.decrypt(enc_x))[0].real)
    enc_x = evaluator.lower_modulus(enc_x, evaluator.scaling_factor ** 6)
    ret = evaluator.multiply(enc_c, enc_x, relin_key)
    ret = evaluator.rescale(ret, evaluator.scaling_factor)
    return ret

def ckks_sqrt_tdoa(enc_c, iterations, num_slot,encoder,encryptor,decryptor,evaluator,relin_key):
    x = 0.5*(1.0/math.sqrt(3000**2+3000**2))
    message = [0]*num_slot
    message[0] = x
    enc_x = encryptor.encrypt(encoder.encode(message,evaluator.scaling_factor))
    plain_half = evaluator.create_constant_plain(0.5)
    value_half = evaluator.multiply_plain(enc_c, plain_half)
    value_half = evaluator.rescale(value_half,evaluator.scaling_factor)
    # value_half = evaluator.lower_modulus(value_half, evaluator.scaling_factor)
    for i in range(0,iterations):
        enc_x_double = evaluator.multiply(enc_x,enc_x,relin_key)
        enc_x_double = evaluator.rescale(enc_x_double,evaluator.scaling_factor)
        enc_x = evaluator.lower_modulus(enc_x,evaluator.scaling_factor)
        enc_x_triple = evaluator.multiply(enc_x_double,enc_x,relin_key)
        enc_x_triple = evaluator.rescale(enc_x_triple,evaluator.scaling_factor)
        enc_x_triple = evaluator.lower_modulus(enc_x_triple,evaluator.scaling_factor**3)
        tmp_2 = evaluator.multiply(value_half,enc_x_triple,relin_key)
        tmp_2 = evaluator.rescale(tmp_2,evaluator.scaling_factor)
        plain_onehalf = evaluator.create_constant_plain(1.5)
        tmp_1 = evaluator.multiply_plain(enc_x,plain_onehalf)
        tmp_1 = evaluator.rescale(tmp_1,evaluator.scaling_factor)
        tmp_1 = evaluator.lower_modulus(tmp_1,evaluator.scaling_factor**4)
        enc_x = evaluator.subtract(tmp_1,tmp_2)
        # value_half = evaluator.lower_modulus(value_half, evaluator.scaling_factor)
        # value_half = evaluator.lower_modulus(value_half, evaluator.scaling_factor)
        # print(math.log(value_half.modulus,2))
        r = random.random()
        plain_r = evaluator.create_constant_plain(-r)
        enc_x_r = evaluator.add_plain(enc_x,plain_r)
        x_r =encoder.decode(decryptor.decrypt(enc_x_r))[0].real
        message[0] = x_r
        re_enc_x_r = encryptor.encrypt(encoder.encode(message,scaling_factor))
        message[0] = r
        plain_r = encoder.encode(message,scaling_factor)
        enc_x = evaluator.add_plain(re_enc_x_r,plain_r)
        # print(encoder.decode(decryptor.decrypt(enc_x))[0].real)
    enc_x = evaluator.lower_modulus(enc_x,evaluator.scaling_factor**4)
    ret = evaluator.multiply(enc_c,enc_x,relin_key)
    ret = evaluator.rescale(ret,evaluator.scaling_factor)
    return ret
def ckks_sqrt(enc_c, iterations, num_slot,encoder,encryptor,decryptor,evaluator,relin_key):
    x = 0.5*(1.0/math.sqrt(3000**2+3000**2))
    message = [0]*num_slot
    message[0] = x
    enc_x = encryptor.encrypt(encoder.encode(message,evaluator.scaling_factor))
    plain_half = evaluator.create_constant_plain(0.5)
    value_half = evaluator.multiply_plain(enc_c, plain_half)
    value_half = evaluator.rescale(value_half,evaluator.scaling_factor)
    # value_half = evaluator.lower_modulus(value_half, evaluator.scaling_factor)
    for i in range(0,iterations):
        enc_x_double = evaluator.multiply(enc_x,enc_x,relin_key)
        enc_x_double = evaluator.rescale(enc_x_double,evaluator.scaling_factor)
        enc_x = evaluator.lower_modulus(enc_x,evaluator.scaling_factor)
        enc_x_triple = evaluator.multiply(enc_x_double,enc_x,relin_key)
        enc_x_triple = evaluator.rescale(enc_x_triple,evaluator.scaling_factor)
        enc_x_triple = evaluator.lower_modulus(enc_x_triple,evaluator.scaling_factor**3)
        tmp_2 = evaluator.multiply(value_half,enc_x_triple,relin_key)
        tmp_2 = evaluator.rescale(tmp_2,evaluator.scaling_factor)
        plain_onehalf = evaluator.create_constant_plain(1.5)
        tmp_1 = evaluator.multiply_plain(enc_x,plain_onehalf)
        tmp_1 = evaluator.rescale(tmp_1,evaluator.scaling_factor)
        tmp_1 = evaluator.lower_modulus(tmp_1,evaluator.scaling_factor**4)
        enc_x = evaluator.subtract(tmp_1,tmp_2)
        # value_half = evaluator.lower_modulus(value_half, evaluator.scaling_factor)
        # value_half = evaluator.lower_modulus(value_half, evaluator.scaling_factor)
        # print(math.log(value_half.modulus,2))
        r = random.random()
        plain_r = evaluator.create_constant_plain(-r)
        enc_x_r = evaluator.add_plain(enc_x,plain_r)
        x_r =encoder.decode(decryptor.decrypt(enc_x_r))[0].real
        message[0] = x_r
        re_enc_x_r = encryptor.encrypt(encoder.encode(message,evaluator.scaling_factor))
        message[0] = r
        plain_r = encoder.encode(message,evaluator.scaling_factor)
        enc_x = evaluator.add_plain(re_enc_x_r,plain_r)
        # print(encoder.decode(decryptor.decrypt(enc_x))[0].real)
    enc_x = evaluator.lower_modulus(enc_x,evaluator.scaling_factor**4)
    ret = evaluator.multiply(enc_c,enc_x,relin_key)
    ret = evaluator.rescale(ret,evaluator.scaling_factor)
    return ret

def divide_two(enc_c,evaluator):
    plain_half = evaluator.create_constant_plain(0.5)
    tmp = evaluator.multiply_plain(enc_c,plain_half)
    return evaluator.rescale(tmp,evaluator.scaling_factor)
def homo_matrix_mul(X,Y,evaluator,relin_key):
    row = len(X)
    col = len(Y[0])
    result = np.empty([row,col],dtype='object')
    # for i in range(row):
    #     for j in range(col):
    #         message = [0] * num_slot
    #         result[i][j] = encryptor.encrypt(encoder.encode(message, 160))

    for i in range(row):
        for j in range(col):
            for k in range(len(Y)):
                tmp = evaluator.multiply(X[i][k], Y[k][j], relin_key)
                tmp = evaluator.rescale(tmp,evaluator.scaling_factor)
                if k == 0:
                    result[i][j] = tmp
                else:
                    result[i][j] = evaluator.add(result[i][j],tmp)

    return  result
def homo_plain_matrix_mul(X,Y,evaluator,type):
    row = len(X)
    col = len(Y[0])
    result = np.empty([row,col],dtype='object')

    for i in range(row):
        for j in range(col):
            for k in range(len(Y)):
                if type == 'left':
                    tmp = evaluator.multiply_plain(Y[k][j],X[i][k])
                else:
                    tmp = evaluator.multiply_plain(X[i][k],Y[k][j])
                tmp = evaluator.rescale(tmp,evaluator.scaling_factor)
                if k == 0:
                    result[i][j] = tmp
                else:
                    result[i][j] = evaluator.add(result[i][j],tmp)

    return  result
def ckks_tdoa_localize(acs,tdoas,num,num_slot,encoder,encryptor,decryptor,evaluator,relin_key):
    acs_sq = [[i[0] ** 2, i[1] ** 2] for i in acs]
    tdoas_sq = [i ** 2 for i in tdoas]
    enc_tdoas = []
    enc_tdoas_sq = [0] * num
    enc_acs = []
    enc_acs_sq = []
    for i in range(0,num-1):
        message = [0] * num_slot
        message[0] = tdoas[i]
        enc_tdoas.append(encryptor.encrypt(encoder.encode(message,evaluator.scaling_factor)))
        message[0] = tdoas_sq[i]
        enc_tdoas_sq[i] = encryptor.encrypt(encoder.encode(message,evaluator.scaling_factor))

    for i in range(0,num):
        message = [0] * num_slot
        message[0] = acs[i][0]
        enc_acs_x = encryptor.encrypt(encoder.encode(message, evaluator.scaling_factor))
        message[0] = acs[i][1]
        enc_acs_y = encryptor.encrypt(encoder.encode(message, evaluator.scaling_factor))
        enc_acs.append([enc_acs_x,enc_acs_y])
        message[0] = acs_sq[i][0]
        enc_acs_x_sq = encryptor.encrypt(encoder.encode(message, evaluator.scaling_factor))
        message[0] = acs_sq[i][1]
        enc_acs_y_sq = encryptor.encrypt(encoder.encode(message, evaluator.scaling_factor))
        enc_acs_sq.append([enc_acs_x_sq,enc_acs_y_sq])

    matrix_a = np.empty([num-1,3],dtype='object')
    for i in range(0,num-1):
        for j in range(0,3):
            if j == 0 or j == 1:
                matrix_a[i,j] = evaluator.subtract(enc_acs[0][j], enc_acs[i+1][j])
            else:
                matrix_a[i,j] = enc_tdoas[i]

    vector_b = np.empty([num-1,1],dtype='object')
    for i in range(0,num-1):
        vector_b[i][0] = evaluator.add(
            evaluator.subtract(
                enc_tdoas_sq[i]
                ,evaluator.add(
                    enc_acs_sq[i+1][0],enc_acs_sq[i+1][1])
            ),enc_acs_sq[0][0])
        vector_b[i][0] = evaluator.add(vector_b[i][0],enc_acs_sq[0][1])

    matrix_a_trans = matrix_a.T
    left = homo_matrix_mul(matrix_a_trans, matrix_a, evaluator, relin_key)
    size = len(left)
    random_matrix = util.get_random_mat(size, 100)
    random_matrix_plain = np.empty([size, size], dtype='object')
    for i in range(0, size):
        for j in range(0, size):
            message = [0] * num_slot
            message[0] = random_matrix[i, j]
            random_matrix_plain[i, j] = encoder.encode(message, evaluator.scaling_factor)

    random_left = homo_plain_matrix_mul(left, random_matrix_plain, evaluator, 'right')

    matrix_random = np.empty([size, size], dtype='float64')
    for i in range(0, size):
        for j in range(0, size):
            matrix_random[i, j] = encoder.decode(decryptor.decrypt(random_left[i, j]))[0].real

    matrix_random_inverse = util.spmat_to_npmat(util.npmat_to_spmat(matrix_random).inv(), 'float')

    enc_mat_r_inv = np.empty([size, size], dtype='object')
    for i in range(0, size):
        for j in range(0, size):
            message = [0] * num_slot
            message[0] = matrix_random_inverse[i, j]
            enc_mat_r_inv[i, j] = encryptor.encrypt(encoder.encode(message, evaluator.scaling_factor))

    left_inverse = homo_plain_matrix_mul(random_matrix_plain, enc_mat_r_inv, evaluator, 'left')
    right = homo_matrix_mul(matrix_a_trans, vector_b, evaluator, relin_key)
    enc_ret = homo_matrix_mul(left_inverse, right, evaluator, relin_key)
    # for i in range(0, len(enc_ret)):
        # tmp = divide_two(enc_ret[i, 0],num_slot,encoder,evaluator,scaling_factor)
        # print(encoder.decode(decryptor.decrypt(enc_ret[i, 0]))[0].real)
    return enc_ret

def ckks_toa_identify(acs,toas,num,num_slot,encoder,encryptor,decryptor,evaluator,relin_key):
    enc_ref_loc = flatten(ckks_toa_localize(acs,toas,num,num_slot,encoder,encryptor,decryptor,evaluator,relin_key).tolist())
    all_subsets = util.subset_all(num,2)
    set_n = len(all_subsets)
    residuals = np.empty(set_n, dtype=object)
    enc_diff_ac = np.empty(num-1, dtype=object)

    enc_acs = []
    enc_toas = []
    for i in range(0,num):
        message = [0] * num_slot
        message[0] = acs[i][0]
        enc_acs_x = encryptor.encrypt(encoder.encode(message, enc_ref_loc[0].scaling_factor))
        enc_acs_x = evaluator.lower_modulus(enc_acs_x,evaluator.scaling_factor**3)
        message[0] = acs[i][1]
        enc_acs_y = encryptor.encrypt(encoder.encode(message, enc_ref_loc[1].scaling_factor))
        enc_acs_y = evaluator.lower_modulus(enc_acs_y,evaluator.scaling_factor**3)
        enc_acs.append([enc_acs_x,enc_acs_y])
        message = [0] * num_slot
        message[0] = toas[i]
        tmp_toas = encryptor.encrypt(encoder.encode(message,evaluator.scaling_factor))
        tmp_toas = evaluator.lower_modulus(tmp_toas,evaluator.scaling_factor**5)
        enc_toas.append(tmp_toas)
        # print("test4:"+"("+str(i)+"):"+str(encoder.decode(decryptor.decrypt(tmp_toas))[0].real))

    for j in range(1,num):
        enc_dist = secure_distance(enc_acs[j],enc_ref_loc,15,num_slot,encoder,encryptor,decryptor,evaluator,relin_key)
        # print("test:"+"("+str(j)+"):"+str(encoder.decode(decryptor.decrypt(enc_dist))[0].real))
        tmp = evaluator.subtract(enc_toas[j],enc_dist)
        # print("test3:"+"("+str(j)+"):"+str(encoder.decode(decryptor.decrypt(tmp))[0].real))
        enc_diff_ac[j - 1] = evaluator.multiply(tmp,tmp,relin_key)
        # print("test2:"+"("+str(j)+"):"+str(encoder.decode(decryptor.decrypt(enc_diff_ac[j - 1]))[0].real))
        enc_diff_ac[j - 1] = evaluator.rescale(enc_diff_ac[j - 1],evaluator.scaling_factor)


    for i in range(0,set_n):
        residuals[i] = encryptor.encrypt(evaluator.create_constant_plain(0))
        residuals[i] = evaluator.lower_modulus(residuals[i],evaluator.scaling_factor**6)
        for j in all_subsets[i]:
            residuals[i] = evaluator.add(residuals[i],enc_diff_ac[j-1])

    acs_res = np.empty(num, dtype=object)
    for i in range(0, num):
        acs_res[i] = encryptor.encrypt(evaluator.create_constant_plain(0))
        acs_res[i] = evaluator.lower_modulus(acs_res[i], evaluator.scaling_factor ** 6)
    # print(acs_res)
    set_num = 0
    for subset in all_subsets:
        for i in range(0,num):
            if i in subset:
                acs_res[i] = evaluator.add(acs_res[i],residuals[set_num])
        set_num = set_num + 1

    ret = []

    for i in range(0, num):
        ret.append(encoder.decode(decryptor.decrypt(acs_res[i]))[0].real)
    return ret

def ckks_toa_localize(acs,toas,num,num_slot,encoder,encryptor,decryptor,evaluator,relin_key):
    acs_sq = [[i[0] ** 2, i[1] ** 2] for i in acs]
    toas_sq =[i ** 2 for i in toas]

    enc_toas_sq = []
    enc_acs = []
    enc_acs_sq = []
    for i in range(0,num):
        message = [0] * num_slot
        message[0] = acs[i][0]
        enc_acs_x = encryptor.encrypt(encoder.encode(message, evaluator.scaling_factor))
        message[0] = acs[i][1]
        enc_acs_y = encryptor.encrypt(encoder.encode(message, evaluator.scaling_factor))
        enc_acs.append([enc_acs_x,enc_acs_y])
        message[0] = acs_sq[i][0]
        enc_acs_x_sq = encryptor.encrypt(encoder.encode(message, evaluator.scaling_factor))
        message[0] = acs_sq[i][1]
        enc_acs_y_sq = encryptor.encrypt(encoder.encode(message, evaluator.scaling_factor))
        enc_acs_sq.append([enc_acs_x_sq,enc_acs_y_sq])
        message[0] = toas_sq[i]
        enc_toas_sq.append(encryptor.encrypt(encoder.encode(message,evaluator.scaling_factor)))

    matrix_a = np.empty([num-1,2],dtype='object')
    message = [0] * num_slot
    message[0] = 2
    plain_2 = encoder.encode(message, evaluator.scaling_factor)

    for i in range(0, num - 1):
        for j in range(0, 2):
            matrix_a[i, j] = evaluator.multiply_plain(evaluator.subtract(enc_acs[i + 1][j], enc_acs[0][j]),plain_2)
            matrix_a[i, j] = evaluator.rescale(matrix_a[i, j],evaluator.scaling_factor)
            # print(encoder.decode(decryptor.decrypt(matrix_a[i, j]))[0].real)
    vector_b = np.empty([num-1,1],dtype='object')
    for i in range(0,num-1):
        vector_b[i][0] = evaluator.subtract(
            evaluator.subtract(
                evaluator.add(
                    enc_acs_sq[i+1][0],enc_acs_sq[i+1][1]
                ),evaluator.add(
                    enc_acs_sq[0][0],enc_acs_sq[0][1])
            ),enc_toas_sq[i+1])
        vector_b[i][0] = evaluator.add(vector_b[i][0],enc_toas_sq[0])
        vector_b[i][0] = evaluator.lower_modulus(vector_b[i][0],evaluator.scaling_factor)

    matrix_a_trans = matrix_a.T
    left = homo_matrix_mul(matrix_a_trans, matrix_a, evaluator, relin_key)
    size = len(left)
    random_matrix = util.get_random_mat(size, 100)
    random_matrix_plain = np.empty([size, size], dtype='object')
    for i in range(0, size):
        for j in range(0, size):
            message = [0] * num_slot
            message[0] = random_matrix[i, j]
            random_matrix_plain[i, j] = encoder.encode(message, evaluator.scaling_factor)

    random_left = homo_plain_matrix_mul(left, random_matrix_plain, evaluator, 'right')

    matrix_random = np.empty([size, size], dtype='float64')
    for i in range(0, size):
        for j in range(0, size):
            matrix_random[i, j] = encoder.decode(decryptor.decrypt(random_left[i, j]))[0].real

    matrix_random_inverse = util.spmat_to_npmat(util.npmat_to_spmat(matrix_random).inv(), 'float')

    enc_mat_r_inv = np.empty([size, size], dtype='object')
    for i in range(0, size):
        for j in range(0, size):
            message = [0] * num_slot
            message[0] = matrix_random_inverse[i, j]
            enc_mat_r_inv[i, j] = encryptor.encrypt(encoder.encode(message, evaluator.scaling_factor))

    left_inverse = homo_plain_matrix_mul(random_matrix_plain, enc_mat_r_inv, evaluator, 'left')
    for i in range(0, size):
        for j in range(0, size):
            left_inverse[i][j] = evaluator.lower_modulus(left_inverse[i][j], evaluator.scaling_factor)
    right = homo_matrix_mul(matrix_a_trans, vector_b, evaluator, relin_key)
    enc_ret = homo_matrix_mul(left_inverse, right, evaluator, relin_key)
    # for i in range(0, len(enc_ret)):
    #     print(encoder.decode(decryptor.decrypt(enc_ret[i,0]))[0].real)
    return enc_ret
def ckks_tdoa_identify(acs,tdoas,num,num_slot,encoder,encryptor,decryptor,evaluator,relin_key):
    enc_ref_loc = flatten(ckks_tdoa_localize(acs,tdoas,num,num_slot,encoder,encryptor,decryptor,evaluator,relin_key).tolist())
    enc_ref_loc = enc_ref_loc[0:2]
    enc_ref_loc[0] = divide_two(enc_ref_loc[0],evaluator)
    enc_ref_loc[1] = divide_two(enc_ref_loc[1],evaluator)
    all_subsets = util.subset_all(num,3)
    set_n = len(all_subsets)
    residuals = np.empty(set_n, dtype=object)
    enc_diff_ac = np.empty(num-1, dtype=object)

    enc_acs = []
    enc_tdoas = []
    for i in range(0,num):
        message = [0] * num_slot
        message[0] = acs[i][0]
        enc_acs_x = encryptor.encrypt(encoder.encode(message, enc_ref_loc[0].scaling_factor))
        enc_acs_x = evaluator.lower_modulus(enc_acs_x,evaluator.scaling_factor**3)
        message[0] = acs[i][1]
        enc_acs_y = encryptor.encrypt(encoder.encode(message, enc_ref_loc[1].scaling_factor))
        enc_acs_y = evaluator.lower_modulus(enc_acs_y,evaluator.scaling_factor**3)
        enc_acs.append([enc_acs_x,enc_acs_y])

    for i in range(0,num-1):
        message = [0] * num_slot
        message[0] = tdoas[i]
        tmp_tdoas = encryptor.encrypt(encoder.encode(message,evaluator.scaling_factor))
        tmp_tdoas = evaluator.lower_modulus(tmp_tdoas,evaluator.scaling_factor**7)
        enc_tdoas.append(tmp_tdoas)
        # print("test4:"+"("+str(i)+"):"+str(encoder.decode(decryptor.decrypt(tmp_toas))[0].real))

    enc_dist_base = secure_distance(enc_acs[0],enc_ref_loc,15,num_slot,encoder,encryptor,decryptor,evaluator,relin_key)

    for j in range(1,num):
        enc_dist = secure_distance(enc_acs[j],enc_ref_loc,15,num_slot,encoder,encryptor,decryptor,evaluator,relin_key)
        enc_diff = secure_abolute(evaluator.subtract(enc_dist,enc_dist_base)
                                  ,15,num_slot,encoder,encryptor,decryptor,evaluator,relin_key)
        # print("test:"+"("+str(j)+"):"+str(encoder.decode(decryptor.decrypt(enc_dist))[0].real))
        tmp = evaluator.subtract(enc_tdoas[j-1],enc_diff)
        # print("test3:"+"("+str(j)+"):"+str(encoder.decode(decryptor.decrypt(tmp))[0].real))
        enc_diff_ac[j - 1] = evaluator.multiply(tmp,tmp,relin_key)
        # print("test2:"+"("+str(j)+"):"+str(encoder.decode(decryptor.decrypt(enc_diff_ac[j - 1]))[0].real))
        enc_diff_ac[j - 1] = evaluator.rescale(enc_diff_ac[j - 1],evaluator.scaling_factor)


    for i in range(0,set_n):
        residuals[i] = encryptor.encrypt(evaluator.create_constant_plain(0))
        residuals[i] = evaluator.lower_modulus(residuals[i],evaluator.scaling_factor**8)
        for j in all_subsets[i]:
            residuals[i] = evaluator.add(residuals[i],enc_diff_ac[j-1])

    acs_res = np.empty(num, dtype=object)
    for i in range(0, num):
        acs_res[i] = encryptor.encrypt(evaluator.create_constant_plain(0))
        acs_res[i] = evaluator.lower_modulus(acs_res[i], evaluator.scaling_factor ** 8)
    # print(acs_res)
    set_num = 0
    for subset in all_subsets:
        for i in range(0,num):
            if i in subset:
                acs_res[i] = evaluator.add(acs_res[i],residuals[set_num])
        set_num = set_num + 1

    # for i in range(0, num):
    #     print(encoder.decode(decryptor.decrypt(acs_res[i]))[0].real)
    ret = []

    for i in range(0, num):
        ret.append(encoder.decode(decryptor.decrypt(acs_res[i]))[0].real)
    return ret

if __name__ == '__main__':
    poly_degree = 16
    ciph_modulus = 1 << 600
    big_modulus = 1 << 1200
    scaling_factor = 1 << 45

    params = CKKSParameters(poly_degree=poly_degree,
                            ciph_modulus=ciph_modulus,
                            big_modulus=big_modulus,
                            scaling_factor=scaling_factor)
    params.print_parameters()
    key_generator = CKKSKeyGenerator(params)
    public_key = key_generator.public_key
    secret_key = key_generator.secret_key
    relin_key = key_generator.relin_key
    encoder = CKKSEncoder(params)
    encryptor = CKKSEncryptor(params, public_key, secret_key)
    decryptor = CKKSDecryptor(params, secret_key)
    evaluator = CKKSEvaluator(params)

    rot_keys = {}
    for i in range(poly_degree//2):
        rot_keys[i] = key_generator.generate_rot_key(i)

    conj_key = key_generator.generate_conj_key()

    ue = [1500,1500]
    acs = np.empty(5, dtype=object)
    toas = np.empty(5, dtype=object)
    tdoas = np.empty(5,dtype=object)

    acs[0] = [[2221, 2273],[701, 1644],[996, 829],[2322, 913]]
    acs[1] = [[1911, 2551],[1171, 2328],[500, 1466],[1436, 3],[2479, 1384]]
    acs[2] = [[2232, 2169],[976, 2877],[244, 1961],[608, 351],[1514, 643],[2269, 1363]]
    acs[3] = [[2531, 2135],[1439, 2589],[1049, 1809],[561, 1352],[849, 474],[1947, 859],[2132, 1423]]
    acs[4] = [[2181, 2026],[1695, 2796],[1273, 2037],[348, 1511],[701, 1080],[1240, 777],[1894, 768],[2301, 1223]]

    # sigma = 30
    # los_error = random.gauss(0, sigma)
    los_err = 5
    nlos_err = 300

    for i in range(1,5):
        toas[i] = util.generate_toa_measure(ue,acs[i],i+4,los_err,nlos_err)
        tdoas[i] = util.generate_tdoa_measure(ue,acs[i],i+4,los_err,nlos_err)

    for i in range(1,5):
        start = time.process_time()
        # ckks_toa_localize(acs[i],toas[i],i+4,poly_degree // 2, encoder,encryptor,decryptor,evaluator,relin_key)
        # util.toa_localize(acs[i],toas[i],i+4)
        # ckks_tdoa_localize(acs[i],tdoas[i],i+4,poly_degree // 2, encoder,encryptor,decryptor,evaluator,relin_key)
        # util.tdoa_localize(acs[i],tdoas[i],i+4)
        end = time.process_time()
        print("processing time:"+ str((end-start)*1000))

    for i in range(1,5):
        # ckks_toa_identify(acs[i],toas[i],i+4,poly_degree // 2,encoder,encryptor,decryptor,evaluator,relin_key)
        # util.toa_identify(acs[i],toas[i],i+4)
        ckks_tdoa_identify(acs[i],tdoas[i],i+4,poly_degree // 2,encoder,encryptor,decryptor,evaluator,relin_key)
        util.tdoa_identify(acs[i],tdoas[i],i+4)



