# PPRP: Preserving Location Privacy for Range-based Positioning in Mobile Networks
This repository contains the source code accompanying the paper "PPRP: Preserving Location Privacy for Range-based Positioning in Mobile Networks," submitted to IEEE Transactions on Mobile Computing.

## Primary Entry Files:
- **main.py**: Main file for testing the proposed basic modules.
- **plain.py**: Main file for testing both the proposed secure positioning scheme and the non-secure positioning scheme.
- **classic_one.py**: Main file for testing PPSPP.
- **classic_two.py**: Main file for testing PHEPP.
- **classic_three.py**: Main file for testing  PHE+PPS.
- **ckks_localization.py**: Main file for testing the FHEBASE algorithm.
- **accuracy_benchmark.py**: Main file for evaluating the accuracy of both the proposed and classical schemes.
- **exec_benchmark.py**: Main file for assessing the computational efficiency of both the proposed and classical schemes.

## Getting Started
The prototype is developed using Python 3.9 in PyCharm. You can import the project into PyCharm or another IDE of your choice that supports Python development.

### Package Requirements
To run the prototype, the following packages must be installed:

- `gmpy2==2.1.5`
- `mpmath==1.3.0`
- `numpy==1.24.3`
- `sympy==1.11.1`

You can install these packages using pip with the following command:

```bash
pip install gmpy2==2.1.5 mpmath==1.3.0 numpy==1.24.3 sympy==1.11.1
