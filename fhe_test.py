#!/usr/bin/env python3

from line_profiler import LineProfiler

from pyfhe.ckks.ckks_decryptor import CKKSDecryptor
from pyfhe.ckks.ckks_encoder import CKKSEncoder
from pyfhe.ckks.ckks_encryptor import CKKSEncryptor
from pyfhe.ckks.ckks_evaluator import CKKSEvaluator
from pyfhe.ckks.ckks_key_generator import CKKSKeyGenerator
from pyfhe.ckks.ckks_parameters import CKKSParameters
import time

if __name__ == '__main__':
    poly_degree = 16
    ciph_modulus = 1 << 1200
    big_modulus = 1 << 2400
    scaling_factor = 1 << 80

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

    message1 = [0] * (poly_degree // 2)
    message2 = [0] * (poly_degree // 2)
    message3 = [0] * (poly_degree // 2)

    message1[0] = -13
    message2[0] = -10
    message3[0] = 10

    start = time.process_time()
    plain1 = encoder.encode(message1, scaling_factor)
    ciph1 = encryptor.encrypt(plain1)
    end = time.process_time()
    print("encryption time:" + str((end - start)*1000))

    plain2 = encoder.encode(message2, scaling_factor)
    ciph2 = encryptor.encrypt(plain2)

    plain3 = encoder.encode(message3, scaling_factor)
    ciph3 = encryptor.encrypt(plain3)

    start = time.process_time()
    ciph_add = evaluator.add(ciph1,ciph2)
    end = time.process_time()
    print("addition time:" + str((end - start)*1000))

    start = time.process_time()
    ciph_prod = evaluator.multiply(ciph1, ciph2, relin_key)
    end = time.process_time()
    print("multiplication time:" + str((end - start)*1000))


    for i in range(0,12):
        ciph3 = evaluator.multiply(ciph3,ciph2,relin_key)

    dec_ciph3 = decryptor.decrypt(ciph3)
    print("--------")
    print(encoder.decode(dec_ciph3))
    print("--------")

    start = time.process_time()
    decrypted_prod = decryptor.decrypt(ciph_prod)
    decoded_prod = encoder.decode(decrypted_prod)
    end = time.process_time()
    print("multi_decryption time:" + str((end - start)*1000))

    print(decoded_prod[0].real)

    start = time.process_time()
    decrypted_add = decryptor.decrypt(ciph_add)
    decoded_add = encoder.decode(decrypted_add)
    end = time.process_time()
    print("multi_decryption time:" + str((end - start)*1000))

    print(decoded_add)