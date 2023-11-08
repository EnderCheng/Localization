import math

import classic_one
import classic_three
import classic_two
import util
import numpy as np

from ckks_localization import ckks_toa_localize, ckks_tdoa_localize, ckks_toa_identify, ckks_tdoa_identify
from phe import paillier
from symmetric import SHE
from pyfhe.ckks.ckks_decryptor import CKKSDecryptor
from pyfhe.ckks.ckks_encoder import CKKSEncoder
from pyfhe.ckks.ckks_encryptor import CKKSEncryptor
from pyfhe.ckks.ckks_evaluator import CKKSEvaluator
from pyfhe.ckks.ckks_key_generator import CKKSKeyGenerator
from pyfhe.ckks.ckks_parameters import CKKSParameters

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

    # for i in range(0, 5):
    #     value_rm = 0
    #     value_rm_ckks = 0
    #     value_rm_phe = 0
    #     value_rm_phepps = 0
    #
    #     [x_0, y_0] = util.toa_localize(acs[i], toas[i], i + 4)
    #     # print("Real location (RM):"+str([x_0,y_0]))
    #     [z_0, z_1] = util.secure_toa_localize(acs[i], toas[i], i + 4, she_0, bit, module)
    #     [xx_0, yy_0] = util.recover_matrix(z_0, z_1, module)
    #     xx_0 = (xx_0 + 0.0) / (10 ** 20)
    #     yy_0 = (yy_0 + 0.0) / (10 ** 20)
    #     # print("our Scheme location (RM):"+str([xx_0,yy_0]))
    #     value_rm = util.distance([x_0, y_0], list([xx_0[0], yy_0[0]]))
    #
    #     enc_ret = ckks_toa_localize(acs[i], toas[i], i + 4, poly_degree // 2, encoder, encryptor, decryptor, evaluator,
    #                                 relin_key)
    #     x_ckks = encoder.decode(decryptor.decrypt(enc_ret[0, 0]))[0].real
    #     y_ckks = encoder.decode(decryptor.decrypt(enc_ret[1, 0]))[0].real
    #     value_rm_ckks = util.distance([x_0, y_0], [x_ckks, y_ckks])
    #
    #     [x_phe, y_phe] = classic_two.encrypt_toa_localize(acs[i],toas[i],i+4,100,public_key,private_key)
    #     value_rm_phe = util.distance([x_0, y_0], [x_phe, y_phe])
    #     ton_ret = classic_three.ton_toa_localize(acs[i],toas[i],i+4,public_key,private_key)
    #     value_rm_phepps = util.distance([x_0, y_0], [ton_ret[0], ton_ret[1]])
    #
    #     print("Scheme Accuracy RM(" + str(i) + "): " + str(math.sqrt(value_rm)))
    #     print("CKKS Accuracy RM(" + str(i) + "): " + str(math.sqrt(value_rm_ckks)))
    #     print("PHE Accuracy RM(" + str(i) + "): " + str(math.sqrt(value_rm_phe)))
    #     print("PHEPPS Accuracy RM(" + str(i) + "): " + str(math.sqrt(value_rm_phepps)))
    #
    #
    # for i in range(0, 5):
    #     value_rdm = 0
    #     value_rdm_ckks = 0
    #     value_rdm_pps = 0
    #
    #     [x_0, y_0, d_0] = util.tdoa_localize(acs[i], tdoas[i], i + 4)
    #     [z_0, z_1] = util.secure_tdoa_localize(acs[i], tdoas[i], i + 4, she_0, bit, module)
    #     [xx_0, yy_0, dd_0] = util.recover_matrix(z_0, z_1, module)
    #     xx_0 = (xx_0 + 0.0) / (2 * 10 ** 20)
    #     yy_0 = (yy_0 + 0.0) / (2 * 10 ** 20)
    #     value_rdm = util.distance([x_0, y_0], list([xx_0[0], yy_0[0]]))
    #
    #     enc_ret = ckks_tdoa_localize(acs[i], tdoas[i], i + 4, poly_degree // 2, encoder, encryptor, decryptor,
    #                                  evaluator,
    #                                  relin_key)
    #     x_ckks = encoder.decode(decryptor.decrypt(enc_ret[0, 0]))[0].real / 2
    #     y_ckks = encoder.decode(decryptor.decrypt(enc_ret[1, 0]))[0].real / 2
    #     value_rdm_ckks = util.distance([x_0, y_0], [x_ckks, y_ckks])
    #     pps_ret = classic_one.down_tdoa_identity(acs[i],tdoas[i],i+4)
    #     value_rdm_pps = util.distance([x_0, y_0], [pps_ret[0], pps_ret[1]])
    #
    #     print("Scheme Accuracy RDM(" + str(i) + "): " + str(math.sqrt(value_rdm)))
    #     print("CKKS Accuracy RDM(" + str(i) + "): " + str(math.sqrt(value_rdm_ckks)))
    #     print("PPS Accuracy RDM(" + str(i) + "): " + str(math.sqrt(value_rdm_pps)))

    # for i in range(1, 5):
    #     value_rm = 0
    #     value_rm_ckks = 0
    #
    #     plain_res_rm = util.toa_identify(acs[i], toas[i], i + 4)
    #     our_res_rm = util.secure_toa_identify(acs[i], toas[i], i + 4,she_0,she_1,bit,module)
    #     for j in range(1,i+4):
    #         value_rm = value_rm + abs(our_res_rm[j]-plain_res_rm[j])/plain_res_rm[j]
    #     value_rm = value_rm/(i+3)
    #     ckks_res_rm = ckks_toa_identify(acs[i],toas[i],i+4,poly_degree // 2,encoder,encryptor,decryptor,evaluator,relin_key)
    #     for j in range(1,i+4):
    #         value_rm_ckks = value_rm_ckks + abs(ckks_res_rm[j]-plain_res_rm[j])/plain_res_rm[j]
    #     value_rm_ckks = value_rm_ckks/(i+3)
    #
    #     print("Scheme RE RM(" + str(i) + "): " + str(value_rm))
    #     print("CKKS RE RM(" + str(i) + "): " + str(value_rm_ckks))

    for i in range(1, 5):
        value_rdm = 0
        value_rdm_ckks = 0

        plain_res_rdm = util.tdoa_identify(acs[i], tdoas[i], i + 4)
        our_res_rdm = util.secure_tdoa_identify(acs[i], tdoas[i], i + 4,she_0,she_1,bit,module)
        for j in range(1,i+4):
            value_rdm = value_rdm + abs(our_res_rdm[j]-plain_res_rdm[j])/plain_res_rdm[j]
        value_rdm = value_rdm/(i+3)
        print("Scheme RE RDM(" + str(i) + "): " + str(value_rdm))

    for i in range(1, 5):
        value_rdm_ckks = 0

        plain_res_rdm = util.tdoa_identify(acs[i], tdoas[i], i + 4)
        ckks_res_rdm = ckks_tdoa_identify(acs[i],tdoas[i],i+4,poly_degree // 2,encoder,encryptor,decryptor,evaluator,relin_key)
        for j in range(1,i+4):
            value_rdm_ckks = value_rdm_ckks + abs(ckks_res_rdm[j]-plain_res_rdm[j])/plain_res_rdm[j]
        value_rdm_ckks = value_rdm_ckks/(i+3)
        print("CKKS RE RDM(" + str(i) + "): " + str(value_rdm_ckks))