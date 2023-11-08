#!/usr/bin/env python3

import util
from phe.util import getprimeover
import math, gmpy2, numpy as np

class SHE(object):

    def __init__(self, k_0, k_1, k_2):
        self.k_0 = k_0
        self.k_3 = k_0 * k_0 - k_0
        self.p = getprimeover(k_0)
        self.q = 1
        for i in range(k_0):
            self.q = self.q * getprimeover(k_0)
        self.L = util.get_number(k_1)
        self.N = self.p * self.q
        self.k_1 = k_1
        self.k_2 = k_2
        self.L_half = self.L // 2
        self.mul_depth = math.floor(k_0 / (k_1 + k_2)) - 1
        self.E01 = self.encrypt(0)
        self.E02 = self.encrypt(0)
        self.coeff = pow(2,k_0-k_1-k_2-5)
        print("bit length input:" + str(self.L_half.bit_length()-1))
        print("mul_depth:" + str(self.mul_depth))

    def encrypt(self, value):
        gamma = gmpy2.mpz(util.get_number(self.k_2))
        r = gmpy2.mpz(util.get_number(self.k_3))
        tmp = gmpy2.add(gmpy2.mul(gamma, gmpy2.mpz(self.L)), gmpy2.mpz(value))
        return int(gmpy2.mod(gmpy2.add(tmp, gmpy2.mul(r, gmpy2.mpz(self.p))), gmpy2.mpz(self.N)))

    def pub_encrypt(self, value):
        gamma_1 = gmpy2.mpz(util.get_number(self.k_2))
        gamma_2 = gmpy2.mpz(util.get_number(self.k_2))
        tmp = gmpy2.add(gmpy2.mul(gamma_1, self.E01), gmpy2.mpz(value))
        tmp = gmpy2.add(gmpy2.mul(gamma_2, self.E02),tmp)
        return int(gmpy2.mod(tmp, gmpy2.mpz(self.N)))

    def decrypt(self, ciphertext):
        tmp = gmpy2.mod(gmpy2.mod(gmpy2.mpz(ciphertext), gmpy2.mpz(self.p)), gmpy2.mpz(self.L))
        L_half = gmpy2.mpz(self.L_half)
        L = gmpy2.mpz(self.L)
        return int(gmpy2.sub(gmpy2.mod(gmpy2.add(tmp, L_half), L), L_half))

    def negate(self, ciphertext):
        return int(gmpy2.mod(gmpy2.mpz(-ciphertext+4*self.E01), gmpy2.mpz(self.N)))

    def negate_end(self, ciphertext):
        return int(gmpy2.mod(gmpy2.mpz(-ciphertext+self.coeff*self.E01), gmpy2.mpz(self.N)))

    def encrypt_matrix(self, mat):
        [row, col] = mat.shape
        enc_mat = np.zeros(mat.shape, dtype='object')
        for i in range(row):
            for j in range(col):
                enc_mat[i, j] = self.encrypt(mat[i, j])
        return enc_mat

    def decrypt_matrix(self, enc_mat):
        [row, col] = enc_mat.shape
        dec_mat = np.zeros(enc_mat.shape, dtype='object')
        for i in range(row):
            for j in range(col):
                dec_mat[i, j] = self.decrypt(enc_mat[i, j])
        return dec_mat



