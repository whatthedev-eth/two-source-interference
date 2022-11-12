import pytest
import os
from starkware.starknet.testing.starknet import Starknet
import asyncio
import numpy as np

#################################################
# Equivalent to parts of constants.cairo contract
#
# Fixed point math constants
SCALE_FP = 10**20
SCALE_FP_SQRT = 10**10
RANGE_CHECK_BOUND = 2**120

# Constants for felts
PRIME = 3618502788666131213697322783095070105623107215331596699973092056135872020481
HALF_PRIME = (
    1809251394333065606848661391547535052811553607665798349986546028067936010240
)

#
# Math constants
#
# Fixed point values
TWO_PI_fp = 628318530 * SCALE_FP / 100000000
PI_fp = 314159265 * SCALE_FP / 100000000
# Non fixed point values
TWO_PI = TWO_PI_fp / SCALE_FP
PI = PI_fp / SCALE_FP

# n = Number of terms in cosine approximation
n = 5


#################################################
# Equivalent to parts of math.cairo contract
#
def theta_shifter(theta):
    # Shifts theta so it is in range -pi <= theta <= +pi
    # using PI instead of np.pi to match cairo file
    theta_abs = abs(theta)
    if theta_abs >= PI:
        cycles_to_shift = 1 + ((theta_abs - PI) / (TWO_PI)) // 1
        theta_abs_shifted = theta_abs - (TWO_PI * cycles_to_shift)
        if theta >= 0:
            theta_shifted = theta_abs_shifted
        else:
            theta_shifted = -theta_abs_shifted
    else:
        theta_shifted = theta
    return theta_shifted


def cosine_n_terms(theta, n):
    # n = number of terms (not order)
    # 2n = order
    # cosine(theta) ~= ((-1)^n)*(theta^(2n))/(2n)!
    #               ~= 1 - theta^2/2! + theta^4/4! - theta^6/6! + ...
    cos_nth = 0
    # Must have -pi <= theta <= +pi for cosine approx., so shift theta as needed
    theta_shifted = theta_shifter(theta)
    for i in range(n):
        power_neg_one = (-1) ** i
        power_theta = theta_shifted ** (2 * i)
        fact = np.math.factorial(2 * i)
        cos_nth += power_neg_one * power_theta / fact
    return cos_nth


#
# Pytest
#

# The path to the contract source code.
CONTRACT_FILE = os.path.join("contracts", "math.cairo")


@pytest.mark.asyncio
async def test():

    starknet = await Starknet.empty()
    contract = await starknet.deploy(
        source=CONTRACT_FILE,
    )
    print()  # print a blank line

    #
    # Test theta_shifter and cosine_8th
    #
    # Pick an odd number of tests to run
    num_tests = 51
    for i in range(num_tests):
        print()

        theta_fp = int((i - ((num_tests - 1) / 2)) * PI_fp / 4)
        theta = theta_fp / SCALE_FP

        print(f"> theta_fp = {theta_fp}")
        print(f"> theta    = {theta}")
        # Add '>' before our print messages to indicate they are our messages

        # Call theta_shifter_fp(theta_fp) and print out return value
        ret = await contract.theta_shifter_fp(theta_fp=theta_fp).call()
        if ret.result.value_fp >= HALF_PRIME:
            result_fp = ret.result.value_fp - PRIME
        else:
            result_fp = ret.result.value_fp
        print(f">         theta_shifter_fp(theta_fp) return: {result_fp}")

        theta_shifter_py = theta_shifter(theta)
        # Scale up by SCALE_FP for comparison
        theta_shifter_py_fp = int(theta_shifter_py * SCALE_FP)
        if theta_shifter_py_fp >= HALF_PRIME:
            result_py_fp = theta_shifter_py_fp - PRIME
        else:
            result_py_fp = theta_shifter_py_fp
        print(f"> theta_shifter_py(theta) return * SCALE_FP: {result_py_fp}")

        # Call cosine_8th() with theta_fp and print out return value
        ret = await contract.cosine_8th_fp(theta_fp=theta_fp).call()
        if ret.result.value_fp >= HALF_PRIME:
            result_fp = ret.result.value_fp - PRIME
        else:
            result_fp = ret.result.value_fp
        print(f">             cosine_8th_fp(theta_fp) return: {result_fp}")

        cos_8th_py = cosine_n_terms(theta, n)
        # Scale up by SCALE_FP for comparison
        cos_8th_py_fp = int(cos_8th_py * SCALE_FP)
        print(f"> cosine_n_terms(theta, 5) return * SCALE_FP: {cos_8th_py_fp}")
        print(f">     or (without mult. by SCALE_FP), return: {cos_8th_py}")
