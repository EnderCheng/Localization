import numpy as np
from sympy import Matrix

import classic_one
import util
from classic_three import ton_toa_localize
from phe import paillier
from pyfhe.ckks.ckks_decryptor import CKKSDecryptor
from pyfhe.ckks.ckks_encoder import CKKSEncoder
from pyfhe.ckks.ckks_encryptor import CKKSEncryptor
from pyfhe.ckks.ckks_evaluator import CKKSEvaluator
from pyfhe.ckks.ckks_key_generator import CKKSKeyGenerator
from pyfhe.ckks.ckks_parameters import CKKSParameters
from symmetric import SHE
from util import split_coordinate, split_vector
import time

def our_ue_operation(acs, toas, tdoas, num,module):

    start = time.process_time()
    acs_sq = [[i[0] ** 2, i[1] ** 2] for i in acs]
    toas_sq =[round(i) ** 2 for i in toas]
    [acs_0,acs_1] = split_coordinate(acs,module)
    [acs_0_sq,acs_1_sq] = split_coordinate(acs_sq,module)
    [toas_0_sq,toas_1_sq] = split_vector(toas_sq, module)
    end = time.process_time()
    print("Scheme DL/UL RM (UE/ACs)"+",AC number("+str(num)+"):"+str((end-start)/ (num+1)))

    start = time.process_time()
    acs_sq = [[i[0] ** 2, i[1] ** 2] for i in acs]
    tdoas =[round(i) for i in tdoas]
    tdoas_sq =[round(i) ** 2 for i in tdoas]
    [acs_0,acs_1] = split_coordinate(acs,module)
    [tdoas_0,tdoas_1] = split_vector(tdoas, module)
    [acs_0_sq,acs_1_sq] = split_coordinate(acs_sq,module)
    [tdoas_0_sq,tdoas_1_sq] = split_vector(tdoas_sq, module)
    end = time.process_time()
    print("Scheme DL/UL RDM (UE/ACs)"+",AC number("+str(num)+"):"+str((end-start)/(num+1)))

def fhe_ue_operation(acs,toas, tdoas,num,num_slot,encryptor,encoder,evaluator):
    start = time.process_time()
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
    end = time.process_time()
    print("CKKS DL/UL RM (UE/ACs)"+",AC number("+str(num)+"):"+str((end-start)/(num+1)))

    start = time.process_time()
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
    end = time.process_time()
    print("CKKS DL/UL RDM (UE/ACs)"+",AC number("+str(num)+"):"+str((end-start)/(num+1)))

def pps_operation(acs,tdoas,num):
    start = time.process_time()
    classic_one.down_tdoa_identity(acs, tdoas, num)
    end = time.process_time()
    print("PPSPP DL RDM (UE/ACs)" + ",AC number(" + str(num) + "):" + str((end - start) / (num+1)))


    start = time.process_time()
    classic_one.up_tdoa_identify(acs, tdoas, num)
    end = time.process_time()
    print("PPSPP UL RDM (UE/ACs)" + ",AC number(" + str(num) + "):" + str((end - start) / (num+1)))

def phe_operation(acs,toas,num,num_pts,public_key):
    start = time.process_time()
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
    end = time.process_time()
    print("PHE UL RM (UE/ACs)" + ",AC number(" + str(num) + "):" + str((end - start) / (num+1)))

def ppsphe_operation(acs,toas,num,public_key,private_key):
    start = time.process_time()
    ton_toa_localize(acs,toas,i+4,public_key,private_key)
    end = time.process_time()
    print("PHE+PPS DL RM (UE/ACs)" + ",AC number(" + str(num) + "):" + str((end - start) / (num+1)))

if __name__ == '__main__':
    ue = [1500, 1500]
    acs = np.empty(5, dtype=object)
    toas = np.empty(5, dtype=object)
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

    k_0 = 500
    k_1 = 200
    k_2 = 100
    bit = 95

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
    for i in range(poly_degree // 2):
        rot_keys[i] = key_generator.generate_rot_key(i)

    conj_key = key_generator.generate_conj_key()

    she_0 = SHE(k_0, k_1, k_2)
    she_1 = SHE(k_0, k_1, k_2)

    public_key, private_key = paillier.generate_paillier_keypair(n_length=2048)

    module = pow(2, bit)

    for i in range(0, 5):
        toas[i] = util.generate_toa_measure(ue, acs[i], i + 4, los_err, nlos_err)
        tdoas[i] = util.generate_tdoa_measure(ue, acs[i], i + 4, los_err, nlos_err)

    for i in range(0,5):
        our_ue_operation(acs[i], toas[i], tdoas[i], i+4, module)
        fhe_ue_operation(acs[i], toas[i], tdoas[i], i+4, poly_degree//2, encryptor, encoder, evaluator)
        pps_operation(acs[i], tdoas[i], i+4)
        phe_operation(acs[i], toas[i], i+4, 100, public_key)
        ppsphe_operation(acs[i],toas[i],i+4,public_key,private_key)
