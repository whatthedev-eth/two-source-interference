import pytest
import os
from starkware.starknet.testing.starknet import Starknet
import asyncio
import numpy as np

# Constants from constants.cairo contract
#
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
TWO_PI = 6283185 * SCALE_FP / 1000000
PI = TWO_PI / 2
# Non fixed point values
TWO_PI_no_fp = TWO_PI / SCALE_FP
PI_no_fp = PI / SCALE_FP


# n = number of terms in cosine approximation
n = 5


#
# Python functions
#
def theta_shifter(theta):
    # shifts theta so it is in range -pi <= theta <= +pi
    # using PI_no_fp instead of np.pi to match cairo file
    theta_abs = abs(theta)
    if theta_abs >= PI_no_fp:
        cycles_to_shift = 1 + ((theta_abs - PI_no_fp) / (TWO_PI_no_fp)) // 1
        theta_abs_shifted = theta_abs - (TWO_PI_no_fp * cycles_to_shift)
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
    # must have -pi <= theta <= +pi for cosine approx., so shift theta as needed
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
    print()  # grab a newline here

    #
    # Test theta_shifter and cosine_8th
    #
    # pick an odd number of tests to run
    num_tests = 51
    for i in range(num_tests):
        print()  # grab a newline here

        theta = int((i - ((num_tests - 1) / 2)) * PI / 4)
        theta_no_fp = theta / SCALE_FP

        print(f"> theta       = {theta}")
        print(f"> theta_no_fp = {theta_no_fp}")
        # add '>' before our print messages to indicate they are our messages

        # Call theta_shifter() with theta and print out return value
        ret = await contract.theta_shifter(theta=theta).call()
        if ret.result.value >= HALF_PRIME:
            result = ret.result.value - PRIME
        else:
            result = ret.result.value
        print(f">       theta_shifter(theta) returns: {result}")

        theta_shifter_py_no_fp = theta_shifter(theta_no_fp)
        # Scale up by SCALE_FP for comparison
        theta_shifter_py = int(theta_shifter_py_no_fp * SCALE_FP)
        if theta_shifter_py >= HALF_PRIME:
            result_py = theta_shifter_py - PRIME
        else:
            result_py = theta_shifter_py
        print(f">    theta_shifter_py(theta) returns: {result_py}")

        # Call cosine_8th() with theta and print out return value
        ret = await contract.cosine_8th(theta=theta).call()
        if ret.result.value >= HALF_PRIME:
            result = ret.result.value - PRIME
        else:
            result = ret.result.value
        print(f">          cosine_8th(theta) returns: {result}")

        cos_8th_py_no_fp = cosine_n_terms(theta_no_fp, n)
        # Scale up by SCALE_FP for comparison
        cos_8th_py = int(cos_8th_py_no_fp * SCALE_FP)
        # if cos_8th_py >= HALF_PRIME:
        #    result_py = cos_8th_py - PRIME
        # else:
        #    result_py = cos_8th_py
        # print(f">   cosine_n_terms(theta, 5) returns: {result_py}")
        print(f">   cosine_n_terms(theta, 5) returns: {cos_8th_py}")
        print(f">            or, with no fp, returns: {cos_8th_py_no_fp}")
