import math

from phe import paillier
import time, util

from symmetric import SHE

value = util.get_number(70)

public_key, private_key = paillier.generate_paillier_keypair(n_length=2048)

start = time.process_time()
for i in range(0, 1000):
    enc_value = public_key.encrypt(value + i)
end = time.process_time()
print("paillier encryption processing time:" + str((end - start)))

start = time.process_time()
for i in range(0, 1000):
    dec_value = private_key.decrypt(enc_value)
end = time.process_time()
print("paillier decryption processing time:" + str((end - start)))

enc_value_1 = public_key.encrypt(12435)
enc_value_2 = public_key.encrypt(678910)

start = time.process_time()
for i in range(0, 1000):
    add_value = enc_value_1 + enc_value_2
end = time.process_time()
print("paillier add-I time:" + str((end - start)))

start = time.process_time()
for i in range(0, 1000):
    add_value_2 = enc_value_1 + 12345
end = time.process_time()
print("paillier add-II time:" + str((end - start)))

start = time.process_time()
for i in range(0, 1000):
    mul_value = enc_value_1 *243564
end = time.process_time()
print("paillier mul-II time:" + str((end - start)))

she = SHE(750, 250, 100)

tau = math.floor((she.k_3-she.k_1-she.k_2)/(she.k_0-she.k_1-she.k_2))-1
print(tau)
start = time.process_time()
for i in range(0, 1000):
    enc_value = she.encrypt(value + i)
end = time.process_time()
print("SHE symmetric encryption processing time:" + str((end - start)))

start = time.process_time()
for i in range(0, 1000):
    enc_value = she.pub_encrypt(value + i)
end = time.process_time()
print("SHE public key encryption processing time:" + str((end - start)))

start = time.process_time()
for i in range(0, 1000):
    dec_value = she.decrypt(enc_value)
end = time.process_time()
print("SHE decryption processing time:" + str((end - start)))

enc_value_1 = she.encrypt(12435)
enc_value_2 = she.encrypt(678910)

start = time.process_time()
for i in range(0, 1000):
    add_value = enc_value_1 + enc_value_2
end = time.process_time()
print("SHE add-I time:" + str((end - start)))

start = time.process_time()
for i in range(0, 1000):
    add_value_2 = enc_value_1 + 12345
end = time.process_time()
print("SHE add-II time:" + str((end - start)))

start = time.process_time()
for i in range(0, 1000):
    mul_value_1 = enc_value_1 *243564
end = time.process_time()
print("SHE mul-II time:" + str((end - start)))

start = time.process_time()
for i in range(0, 1000):
    mul_value_2 = enc_value_1 *enc_value_2
end = time.process_time()
print("paillier mul-I time:" + str((end - start)))

