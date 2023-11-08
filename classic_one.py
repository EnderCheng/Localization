#!/usr/bin/env python3
import time

import numpy as np
from sympy import Matrix

import util
import sympy as sp
from operator import add


def up_tdoa_identify(acs, tdoas, num):
    min = -2 ** 32
    max = 2 ** 32
    # data collect step-1
    value_1_agg = Matrix(2, 1, [0, 0])
    value_2_agg = Matrix([0])
    value_3_agg = Matrix([0])
    value_4_agg = Matrix([0])
    value_7_agg = Matrix([0])
    value_8_agg = Matrix([0])

    value_1 = [0] * (num - 1)
    value_2 = [0] * (num - 1)
    value_3 = [0] * (num - 1)
    value_4 = [0] * (num - 1)
    value_7 = [0] * (num - 1)
    value_8 = [0] * (num - 1)
    value_1_r = [0] * (num - 1)
    value_2_r = [0] * (num - 1)
    value_3_r = [0] * (num - 1)
    value_4_r = [0] * (num - 1)
    value_7_r = [0] * (num - 1)
    value_8_r = [0] * (num - 1)

    for i in range(0, num - 1):
        value_1[i] = Matrix(2, 1, acs[i + 1])
        value_2[i] = value_1[i].transpose() * value_1[i]
        value_3[i] = Matrix([tdoas[i]])
        value_4[i] = Matrix([tdoas[i] ** 2])
        value_7[i] = Matrix([tdoas[i] ** 2])
        value_8[i] = Matrix([-tdoas[i] ** 3])
        value_1_r[i] = util.generate_matrix_pps(2, 1, num - 1, min, max)
        value_2_r[i] = util.generate_matrix_pps(1, 1, num - 1, min, max)
        value_3_r[i] = util.generate_matrix_pps(1, 1, num - 1, min, max)
        value_4_r[i] = util.generate_matrix_pps(1, 1, num - 1, min, max)
        value_7_r[i] = util.generate_matrix_pps(1, 1, num - 1, min, max)
        value_8_r[i] = util.generate_matrix_pps(1, 1, num - 1, min, max)

    # data aggregation step-1
    for i in range(0, num - 1):
        value_1_agg = value_1_agg + value_1[i]
        value_2_agg = value_2_agg + value_2[i]
        value_3_agg = value_3_agg + value_3[i]
        value_4_agg = value_4_agg + value_4[i]
        value_7_agg = value_7_agg + value_7[i]
        value_8_agg = value_8_agg + value_8[i]
        for j in range(0, num - 1):
            value_1_agg = value_1_agg + value_1_r[i][j]
            value_2_agg = value_2_agg + value_2_r[i][j]
            value_3_agg = value_3_agg + value_3_r[i][j]
            value_4_agg = value_4_agg + value_4_r[i][j]
            value_7_agg = value_7_agg + value_7_r[i][j]
            value_8_agg = value_8_agg + value_8_r[i][j]
    A22 = value_7_agg
    B22 = value_8_agg

    # data collect step-2
    value_5 = [0] * num
    value_6 = [0] * num
    value_9 = [0] * num
    value_10 = [0] * num
    value_11 = [0] * num
    value_5_r = [0] * num
    value_6_r = [0] * num
    value_9_r = [0] * num
    value_10_r = [0] * num
    value_11_r = [0] * num
    value_5_agg = Matrix([[0, 0], [0, 0]])
    value_6_agg = Matrix(2, 1, [0, 0])
    value_9_agg = Matrix(2, 1, [0, 0])
    value_10_agg = Matrix(2, 1, [0, 0])
    value_11_agg = Matrix([0])

    for i in range(0, num - 1):
        value_5[i] = value_1[i] * value_1[i].transpose()
        value_6[i] = value_1[i] * (value_3[i][0, 0])
        value_9[i] = value_1[i] * value_1[i].transpose() * value_1[i]
        value_10[i] = value_1[i] * (-value_4[i][0, 0])
        value_11[i] = value_1[i].transpose() * value_1[i] * (value_3[i][0, 0])
    tmp = Matrix(2, 1, acs[0])
    tmp_t = tmp.transpose()
    value_5[num - 1] = (num - 1) * tmp * tmp_t - value_1_agg * tmp_t - tmp * value_1_agg.transpose()
    value_6[num - 1] = tmp * (-value_3_agg[0, 0])
    value_9[num - 1] = (num - 1) * tmp * tmp_t * tmp - value_1_agg * tmp_t * tmp - tmp * value_2_agg
    value_10[num - 1] = tmp * value_7_agg[0, 0]
    value_11[num - 1] = tmp_t * tmp * (-value_3_agg[0, 0])

    for i in range(0, num):
        value_5_r[i] = util.generate_matrix_pps(2, 2, num, min, max)
        value_6_r[i] = util.generate_matrix_pps(2, 1, num, min, max)
        value_9_r[i] = util.generate_matrix_pps(2, 1, num, min, max)
        value_10_r[i] = util.generate_matrix_pps(2, 1, num, min, max)
        value_11_r[i] = util.generate_matrix_pps(1, 1, num, min, max)

    for i in range(0, num):
        value_5_agg = value_5_agg + value_5[i]
        value_6_agg = value_6_agg + value_6[i]
        value_9_agg = value_9_agg + value_9[i]
        value_10_agg = value_10_agg + value_10[i]
        value_11_agg = value_11_agg + value_11[i]
        for j in range(0, num):
            value_5_agg = value_5_agg + value_5_r[i][j]
            value_6_agg = value_6_agg + value_6_r[i][j]
            value_9_agg = value_9_agg + value_9_r[i][j]
            value_10_agg = value_10_agg + value_10_r[i][j]
            value_11_agg = value_11_agg + value_11_r[i][j]

    A11 = value_5_agg
    A12 = value_6_agg
    A21 = value_6_agg.transpose()
    B11 = value_9_agg
    B12 = value_10_agg
    B21 = value_11_agg

    left_top = A11.row_join(A12)
    left_btm = A21.row_join(A22)
    left = left_top.col_join(left_btm)

    right_top = B11 + B12
    right_btm = B21 + B22
    right = right_top.col_join(right_btm)

    # print(left)
    # print(right)

    ret = left.inv() * right
    ret = ret / 2
    return ret


def down_tdoa_identity(acs, tdoas, num):
    min = -2 ** 10
    max = 2 ** 10
    value_1 = [0] * (num - 1)
    value_2 = [0] * (num - 1)
    value_3 = [0] * (num - 1)
    value_5 = [0] * num
    value_9 = [0] * num
    value_7 = [0] * (num - 1)
    value_8 = [0] * (num - 1)
    value_1_r = [0] * (num - 1)
    value_2_r = [0] * (num - 1)
    value_5_r = [0] * num
    value_9_r = [0] * num
    A22 = Matrix([0])
    B22 = Matrix([0])
    value_1_agg = Matrix(2, 1, [0, 0])
    value_2_agg = Matrix([0])
    value_3_agg = Matrix([0])
    value_5_agg = Matrix([[0, 0], [0, 0]])
    value_9_agg = Matrix(2, 1, [0, 0])

    for i in range(0, num - 1):
        value_3[i] = Matrix([tdoas[i]])
        value_7[i] = Matrix([tdoas[i] ** 2])
        value_8[i] = Matrix([-tdoas[i] ** 3])
        A22 = A22 + value_7[i]
        B22 = B22 + value_8[i]
        value_3_agg = value_3_agg + value_3[i]

    for i in range(0, num - 1):
        value_1[i] = Matrix(2, 1, acs[i + 1])
        value_2[i] = value_1[i].transpose() * value_1[i]
        value_5[i] = value_1[i] * value_1[i].transpose()
        value_9[i] = value_1[i] * value_1[i].transpose() * value_1[i]
        value_1_r[i] = util.generate_matrix_pps(2, 1, num - 1, min, max)
        value_2_r[i] = util.generate_matrix_pps(1, 1, num - 1, min, max)

    for i in range(0, num - 1):
        value_1_agg = value_1_agg + value_1[i]
        value_2_agg = value_2_agg + value_2[i]
        for j in range(0, num - 1):
            value_1_agg = value_1_agg + value_1_r[i][j]
            value_2_agg = value_2_agg + value_2_r[i][j]

    tmp = Matrix(2, 1, acs[0])
    tmp_t = tmp.transpose()
    value_5[num - 1] = (num - 1) * tmp * tmp_t - value_1_agg * tmp_t - tmp * value_1_agg.transpose()
    value_9[num - 1] = (num - 1) * tmp * tmp_t * tmp - value_1_agg * tmp_t * tmp - tmp * value_2_agg

    for i in range(0, num):
        value_5_r[i] = util.generate_matrix_pps(2, 2, num, min, max)
        value_9_r[i] = util.generate_matrix_pps(2, 1, num, min, max)

    for i in range(0, num):
        value_5_agg = value_5_agg + value_5[i]
        value_9_agg = value_9_agg + value_9[i]
        for j in range(0, num):
            value_5_agg = value_5_agg + value_5_r[i][j]
            value_9_agg = value_9_agg + value_9_r[i][j]

    A11 = value_5_agg
    B11 = value_9_agg

    radius1 = [0] * (num - 1)
    radius2 = [0] * (num - 1)
    for i in range(0, num - 1):
        radius1[i] = value_3[i] / value_3_agg[0, 0]
        radius2[i] = value_7[i] / A22[0, 0]

    sigma1_r = [0] * num
    sigma2_r = [0] * num
    rho_r = [0] * num
    for i in range(0, num):
        sigma1_r[i] = util.generate_matrix_pps(2, 1, num, min, max)
        sigma2_r[i] = util.generate_matrix_pps(2, 1, num, min, max)
        rho_r[i] = util.generate_matrix_pps(1, 1, num, min, max)

    sigma1 = [0] * num
    sigma2 = [0] * num
    rho = [0] * num
    for i in range(0, num):
        tmp_1 = Matrix(2, 1, [0, 0])
        tmp_2 = Matrix(2, 1, [0, 0])
        tmp_3 = Matrix(1, 1, [0])
        for j in range(0, num):
            tmp_1 = tmp_1 + sigma1_r[j][i]
            tmp_2 = tmp_2 + sigma2_r[j][i]
            tmp_3 = tmp_3 + rho_r[j][i]
        sigma1[i] = tmp_1
        sigma2[i] = tmp_2
        rho[i] = tmp_3

    chi1 = [0] * num
    chi2 = [0] * num
    eta = [0] * num
    for i in range(0, num - 1):
        chi1[i] = value_1[i] + (sigma1[i + 1] / radius1[i][0, 0])
        chi2[i] = -1 * value_1[i] + (sigma2[i + 1] / radius2[i][0, 0])
        eta[i] = (value_1[i].transpose() * value_1[i]) + (rho[i + 1] / radius1[i][0, 0])

    chi_0_1 = -1 * tmp + sigma1[0]
    chi_0_2 = tmp + sigma2[0]
    eta_0 = -1 * tmp_t * tmp + rho[0]

    A12_left = Matrix(2, 1, [0, 0])
    A12_right = Matrix(2, 1, [0, 0])
    B12_left = Matrix(2, 1, [0, 0])
    B12_right = Matrix(2, 1, [0, 0])
    B21_left = Matrix(1, 1, [0])
    B21_right = Matrix(1, 1, [0])

    for i in range(0, num - 1):
        A12_left = A12_left + chi1[i] * value_3[i][0, 0]
        B12_left = B12_left + chi2[i] * value_7[i][0, 0]
        B21_left = B21_left + eta[i] * value_3[i][0, 0]

    A12_right = A12_right + chi_0_1 * value_3_agg[0, 0]
    B12_right = B12_right + chi_0_2 * A22[0, 0]
    B21_right = B21_right + eta_0 * value_3_agg[0, 0]

    A12 = A12_left + A12_right
    A21 = A12.transpose()
    B12 = B12_left + B12_right
    B21 = B21_left + B21_right

    left_top = A11.row_join(A12)
    left_btm = A21.row_join(A22)
    left = left_top.col_join(left_btm)

    right_top = B11 + B12
    right_btm = B21 + B22
    right = right_top.col_join(right_btm)

    # print(left)
    # print(right)

    ret = left.inv() * right
    ret = ret / 2
    # print(ret)
    return ret


if __name__ == '__main__':
    ue = [1500, 1500]
    acs = np.empty(5, dtype=object)
    # toas = np.empty(5, dtype=object)
    tdoas = np.empty(5, dtype=object)

    acs[0] = [[2221, 2273], [701, 1644], [996, 829], [2322, 913]]
    acs[1] = [[1911, 2551], [1171, 2328], [500, 1466], [1436, 3], [2479, 1384]]
    acs[2] = [[2232, 2169], [976, 2877], [244, 1961], [608, 351], [1514, 643], [2269, 1363]]
    acs[3] = [[2531, 2135], [1439, 2589], [1049, 1809], [561, 1352], [849, 474], [1947, 859], [2132, 1423]]
    acs[4] = [[2181, 2026], [1695, 2796], [1273, 2037], [348, 1511], [701, 1080], [1240, 777], [1894, 768],
              [2301, 1223]]

    # sigma = 30
    # los_error = random.gauss(0, sigma)
    los_err = 5
    nlos_err = 300

    for i in range(4, 5):
        tdoas[i] = util.generate_tdoa_measure(ue, acs[i], i + 4, los_err, nlos_err)

    for i in range(4, 5):
        # start = time.process_time()
        up_tdoa_identify(acs[i], tdoas[i], i + 4)
        # end = time.process_time()
        # print("up tdoa time:" + str((end - start) * 1000))
        # util.tdoa_localize(acs[i],tdoas[i],i+4)
        start = time.process_time()
        down_tdoa_identity(acs[i], tdoas[i], i + 4)
        end = time.process_time()
        print("down tdoa time:" + str((end - start) * 1000))
