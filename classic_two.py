#!/usr/bin/env python3
import numpy as np
from sympy import Matrix

import util
from phe import paillier

def homo_plain_matrix_mul(X,Y):
    row = X.shape[0]
    col = len(Y[0])
    result = np.empty([row,col],dtype='object')

    for i in range(row):
        for j in range(col):
            for k in range(len(Y)):
                tmp = float(X[i,k])*Y[k,j]
                if k == 0:
                    result[i][j] = tmp
                else:
                    result[i][j] = result[i][j]+tmp

    return  result

def plain_toa_localize(acs, toas, num, num_pts):
    [A, b] = util.outside_circle_to_vectors(acs[0][0], acs[0][1], toas[0], num_pts)
    A = Matrix(A)
    b = Matrix(b)
    # step = A.copy()
    for i in range(1,num):
        [A_add, b_add] = util.outside_circle_to_vectors(acs[i][0], acs[i][1], toas[i], num_pts)
        A = A.col_join(Matrix(A_add))
        b = b.col_join(Matrix(b_add))
    ret = ((A.transpose()*A).inv())*A.transpose()*b
    print(ret)
    return ret

def encrypt_toa_localize(acs, toas, num, num_pts, public_key, private_key):

    [A, b] = util.outside_circle_to_vectors(acs[0][0], acs[0][1], toas[0], num_pts)
    A = Matrix(A)
    b = Matrix(b)
    for i in range(1,num):
        [A_add, b_add] = util.outside_circle_to_vectors(acs[i][0], acs[i][1], toas[i], num_pts)
        A = A.col_join(Matrix(A_add))
        b = b.col_join(Matrix(b_add))

    enc_b = util.spmat_to_npmat(b, 'object')

    for i in range(len(b)):
        enc_b[i] = public_key.encrypt(round(float(b[i]),5))

    A_t = A.transpose()
    left = (A_t*A).inv()*A_t
    tmp = homo_plain_matrix_mul(left,enc_b)
    # print(private_key.decrypt(tmp[0][0]))
    # print(private_key.decrypt(tmp[1][0]))
    ret =  [private_key.decrypt(tmp[0][0]),private_key.decrypt(tmp[1][0])]
    return ret


if __name__ == '__main__':
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

    # [A,b] = util.outside_circle_to_vectors(0,0,5,10)
    # for i in range(len(A)):
    #     print("("+str(A[i][0])+")*x+("+str(A[i][1])+")*y<="+str(b[i]))
    public_key, private_key = paillier.generate_paillier_keypair(n_length=2048)

    for i in range(0,1):
        toas[i] = util.generate_toa_measure(ue,acs[i],i+4,los_err,nlos_err)

    for i in range(0,1):
        plain_toa_localize(acs[i],toas[i],i+4,100)
        # util.toa_localize(acs[i],toas[i],i+4)
        encrypt_toa_localize(acs[i],toas[i],i+4,100,public_key,private_key)