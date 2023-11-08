
import numpy as np
import util
from ckks_localization import ckks_toa_localize, ckks_tdoa_localize, ckks_toa_identify, ckks_tdoa_identify
from classic_one import up_tdoa_identify, down_tdoa_identity
from classic_three import ton_toa_localize
from classic_two import encrypt_toa_localize
from phe import paillier
from pyfhe.ckks.ckks_decryptor import CKKSDecryptor
from pyfhe.ckks.ckks_encoder import CKKSEncoder
from pyfhe.ckks.ckks_encryptor import CKKSEncryptor
from pyfhe.ckks.ckks_evaluator import CKKSEvaluator
from pyfhe.ckks.ckks_key_generator import CKKSKeyGenerator
from pyfhe.ckks.ckks_parameters import CKKSParameters
from symmetric import SHE
import time

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

    start = time.process_time()
    for i in range(0,5):
        toas[i] = util.generate_toa_measure(ue, acs[i], i + 4, los_err, nlos_err)
        tdoas[i] = util.generate_tdoa_measure(ue, acs[i], i + 4, los_err, nlos_err)
    end = time.process_time()

    # for i in range(0,5):
        # start = time.process_time()
        # for j in range(100):
        #     util.secure_toa_localize(acs[i],toas[i],i+4,she_0,bit,module)
        # end = time.process_time()
        # print((end-start)*10,end=' ')

        # start = time.process_time()
        # for j in range(1):
        #     ckks_toa_localize(acs[i],toas[i],i+4,poly_degree // 2, encoder,encryptor,decryptor,evaluator,relin_key)
        # end = time.process_time()
        # print((end-start)*1000,end=' ')

        # start = time.process_time()
        # for j in range(100):
        #     ton_toa_localize(acs[i], toas[i], i + 4, public_key, private_key)
        # end = time.process_time()
        # print((end-start)*10,end=' ')

        # start = time.process_time()
        # for j in range(1):
        #     encrypt_toa_localize(acs[i], toas[i], i + 4, 100, public_key, private_key)
        # end = time.process_time()
        # print((end-start)*1000,end=' ')

        # start = time.process_time()
        # for j in range(100):
        #     util.secure_tdoa_localize(acs[i],tdoas[i],i+4,she_0,bit,module)
        # end = time.process_time()
        # print((end-start)*10,end=' ')

        # start = time.process_time()
        # for j in range(1):
        #     ckks_tdoa_localize(acs[i],tdoas[i],i+4,poly_degree // 2, encoder,encryptor,decryptor,evaluator,relin_key)
        # end = time.process_time()
        # print((end-start)*1000,end=' ')

        # start = time.process_time()
        # for j in range(1000):
        #     up_tdoa_identify(acs[i],tdoas[i],i+4)
        # end = time.process_time()
        # print((end-start),end=' ')

        # start = time.process_time()
        # for j in range(100):
        #     down_tdoa_identity(acs[i], tdoas[i], i + 4)
        # end = time.process_time()
        # print((end-start)*10,end=' ')


    for i in range(1,5):
        # util.toa_identify(acs[i], toas[i], i+4)

        # start = time.process_time()
        # for j in range(1):
        #     util.secure_toa_identify(acs[i], toas[i], i + 4, she_0, she_1, bit, module)
        # end = time.process_time()
        # print((end-start)*1000,end=' ')

        # start = time.process_time()
        # for j in range(1):
        #     ckks_toa_identify(acs[i],toas[i],i+4,poly_degree // 2,encoder,encryptor,decryptor,evaluator,relin_key)
        # end = time.process_time()
        # print((end-start)*1000,end=' ')

        # start = time.process_time()
        # for j in range(1):
        #     util.secure_tdoa_identify(acs[i], tdoas[i], i + 4,she_0,she_1,bit,module)
        # end = time.process_time()
        # print((end-start)*1000,end=' ')

        # util.tdoa_identify(acs[i],tdoas[i], i+4)

        start = time.process_time()
        for j in range(1):
            ckks_tdoa_identify(acs[i],tdoas[i],i+4,poly_degree // 2,encoder,encryptor,decryptor,evaluator,relin_key)
        end = time.process_time()
        print((end-start)*1000,end=' ')