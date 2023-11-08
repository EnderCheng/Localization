#!/usr/bin/env python3
import random

from line_profiler import LineProfiler

import util
import numpy as np

from symmetric import SHE

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

    she_0 = SHE(500,200,100)
    she_1 = SHE(500,200,100)
    bit = 95
    module = pow(2,bit)

    for i in range(4,5):
        toas[i] = util.generate_toa_measure(ue,acs[i],i+4,los_err,nlos_err)
        tdoas[i] = util.generate_tdoa_measure(ue,acs[i],i+4,los_err,nlos_err)

    for i in range(4,5):
        util.toa_localize(acs[i], toas[i], i+4)
        util.secure_toa_localize(acs[i],toas[i],i+4,she_0,bit,module)
        util.tdoa_localize(acs[i], tdoas[i], i+4)
        util.secure_tdoa_localize(acs[i],tdoas[i],i+4,she_0,bit,module)

    for i in range(4,5):
        util.toa_identify(acs[i], toas[i], i+4)
        util.secure_toa_identify(acs[i], toas[i], i + 4,she_0,she_1,bit,module)
        util.tdoa_identify(acs[i],tdoas[i], i+4)
        util.secure_tdoa_identify(acs[i], tdoas[i], i + 4,she_0,she_1,bit,module)


































