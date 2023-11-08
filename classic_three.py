#!/usr/bin/env python3
import numpy as np
from sympy import Matrix, zeros

import util
from phe import paillier


def ton_toa_localize(acs, toas, num, public_key, private_key):
    min = -2**32
    max = 2**32

    p_r = [0]* num
    v_r = [0]* num
    w_r = [0]* (num-1)
    t_r = [0]* (num-1)
    z_r = [0]* num

    for i in range(0,num):
        p_r[i] = util.generate_matrix_pps(2,2,num,min,max)
        v_r[i] = util.generate_matrix_pps(2,1,num,min,max)
        z_r[i] = util.generate_matrix_pps(2,1,num,min,max)
    for i in range(0, num-1):
        w_r[i] = util.generate_matrix_pps(2,1,num-1,min,max)
        t_r[i] = util.generate_matrix_pps(1,1,num-1,min,max)

    P = [0]* num
    V = [0]* num
    W = [0]* (num-1)
    T = [0]* (num-1)
    Z = [0]* num
    for i in range(0,num):
        tmp_1 = Matrix([[0,0],[0,0]])
        tmp_2 = Matrix(2,1,[0,0])
        tmp_3 = Matrix(2,1,[0,0])
        for j in range(0,num):
            tmp_1 = tmp_1 + p_r[j][i]
            tmp_2 = tmp_2 + v_r[j][i]
            tmp_3 = tmp_3 + z_r[j][i]
        P[i] = tmp_1
        V[i] = tmp_2
        Z[i] = tmp_3
    for i in range(0, num-1):
        tmp_3 = Matrix(2, 1, [0, 0])
        tmp_4 = Matrix(1, 1, [0])
        for j in range(0, num-1):
            tmp_3 = tmp_3 + w_r[j][i]
            tmp_4 = tmp_4 + t_r[j][i]
        W[i] = tmp_3
        T[i] = tmp_4

    enc_toas = [0] * (num-1)
    enc_one =public_key.encrypt(1)
    MDZ = [0]* num
    d_sum = 0
    for i in range(0,num-1):
        d_tmp = toas[i+1]**2-toas[0]**2
        d_sum = d_sum + d_tmp
        enc_toas[i] = public_key.encrypt(-d_tmp)

        values_tmp = [0] * 2
        for j in range(0,2):
            values_tmp[j] = enc_toas[i]*acs[i+1][j]+enc_one*int(Z[i+1][j])

        for j in range(0,2):
            if j == 0:
                MDZ[i+1] = Matrix(1,1,[private_key.decrypt(values_tmp[j])])
            else:
                MDZ[i+1] = MDZ[i+1].col_join(Matrix(1,1,[private_key.decrypt(values_tmp[j])]))

    enc_d_sum = public_key.encrypt(d_sum)
    values_sum = [0] * 2
    for j in range(0, 2):
        values_sum[j] =  enc_d_sum*acs[0][j]+enc_one*int(Z[0][j])

    for j in range(0, 2):
        if j == 0:
            MDZ[0] = Matrix(1, 1, [private_key.decrypt(values_sum[j])])
        else:
            MDZ[0] = MDZ[0].col_join(Matrix(1, 1, [private_key.decrypt(values_sum[j])]))

    MDZ_sum = Matrix(2, 1, [0,0])
    for i in range(num):
        MDZ_sum = MDZ_sum + MDZ[i]

    ####

    ####

    value_1 = [0] * num
    h = [0] * num
    Omega = [0] * num
    Phi = [0] * num
    alpha = [0] * (num-1)
    beta = [0] * (num-1)
    for i in range(0,num-1):
        value_1[i] = Matrix(1,2,acs[i+1])
        h[i] = acs[i+1][0]**2+acs[i+1][1]**2
        # h[i] = acs[i+1][0]**2+acs[i+1][1]**2 - toas[i+1]**2
        Omega[i] = (value_1[i].transpose()*value_1[i])+P[i]
        Phi[i] = h[i]*value_1[i].transpose()+V[i]
        alpha[i] = value_1[i].transpose()+W[i]
        beta[i] = h[i] + T[i][0,0]

    alpha_agg = Matrix(2, 1, [0, 0])
    beta_agg = 0
    for i in range(0,num-1):
        alpha_agg = alpha_agg + alpha[i]
        beta_agg = beta_agg + beta[i]

    tmp = value_1[num-1] = Matrix(1,2,acs[0])
    tmp_t = tmp.transpose()
    Omega[num-1] = (num-1)*tmp_t*tmp - alpha_agg*tmp - tmp_t*alpha_agg.transpose() + P[num-1]
    # h[num-1] = h_m = acs[0][0]**2+acs[0][1]**2 - toas[0]**2
    h[num - 1] = h_m = acs[0][0]**2+acs[0][1]**2
    Phi[num-1] = (num-1)*h_m*tmp_t-h_m*alpha_agg-beta_agg*tmp_t+V[num-1]

    ######
    # value_1 = [0] * (num-1)
    # h = [0] *(num-1)
    # tmp = Matrix(2,1,acs[0])
    # tmp_t = tmp.transpose()
    # for i in range(0,num-1):
    #     value_1[i] = Matrix(2,1,acs[i+1])
    #     h[i] = acs[i + 1][0] ** 2 + acs[i + 1][1] ** 2 - toas[i + 1]
    # value_5_agg = Matrix([[0,0],[0,0]])
    # value_1_agg =  Matrix(2,1,[0,0])
    # value_5 = [0] * num
    # for i in range(0,num-1):
    #     value_1_agg = value_1_agg + value_1[i]
    #     value_5[i] = value_1[i]*value_1[i].transpose()
    # value_5[num - 1] = (num - 1) * tmp*tmp_t   - value_1_agg*tmp_t - tmp*value_1_agg.transpose()
    # h_m = acs[0][0] ** 2 + acs[0][1] ** 2 - toas[0]
    # for i in range(0,num):
    #     value_5_agg = value_5_agg + value_5[i]
    # print(4*value_5_agg)
    # h_all = 0
    # phi_all = Matrix(2,1,[0,0])
    # for i in range(0,num-1):
    #     h_all = h_all + h[i]
    #     phi_all = phi_all + h[i]*value_1[i].transpose()
    # test = (num-1)*h_m*tmp_t+phi_all-h_m*alpha_agg-h_all*tmp_t
    # print("xxxxxx")
    # print(test*2)
    # print("xxxxxx")
    ######

    Omega_agg = Matrix([[0,0],[0,0]])
    Phi_agg = Matrix(2, 1, [0, 0])
    for i in range(0, num):
        Omega_agg = Omega_agg + Omega[i]
        Phi_agg = Phi_agg + Phi[i]

    ret = (Omega_agg*4).inv()*((Phi_agg+MDZ_sum)*2)
    # print(ret)
    return ret

if __name__ == '__main__':
    ue = [1500,1500]
    acs = np.empty(5, dtype=object)
    toas = np.empty(5, dtype=object)

    acs[0] = [[2221, 2273],[701, 1644],[996, 829],[2322, 913]]
    acs[1] = [[1911, 2551],[1171, 2328],[500, 1466],[1436, 3],[2479, 1384]]
    acs[2] = [[2232, 2169],[976, 2877],[244, 1961],[608, 351],[1514, 643],[2269, 1363]]
    acs[3] = [[2531, 2135],[1439, 2589],[1049, 1809],[561, 1352],[849, 474],[1947, 859],[2132, 1423]]
    acs[4] = [[2181, 2026],[1695, 2796],[1273, 2037],[348, 1511],[701, 1080],[1240, 777],[1894, 768],[2301, 1223]]

    los_err = 5
    nlos_err = 300

    public_key, private_key = paillier.generate_paillier_keypair(n_length=2048)

    for i in range(0,5):
        toas[i] = util.generate_toa_measure(ue,acs[i],i+4,los_err,nlos_err)

    for i in range(0,5):
        ret_1 = util.toa_localize(acs[i],toas[i],i+4)
        ret_2 = ton_toa_localize(acs[i],toas[i],i+4,public_key,private_key)
        print(ret_1)
        print(ret_2)