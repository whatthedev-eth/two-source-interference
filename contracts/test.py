import pytest
import os
from starkware.starknet.testing.starknet import Starknet
import asyncio
import numpy as np

# Constants from constants.cairo contract
#
# Fixed point math constants
#
SCALE_FP = 10**20
SCALE_FP_SQRT = 10**10
RANGE_CHECK_BOUND = 2**120
# these are not used yet
PRIME = 3618502788666131213697322783095070105623107215331596699973092056135872020481
HALF_PRIME = (
    1809251394333065606848661391547535052811553607665798349986546028067936010240
)

#
# Math constants
#
TWO_PI = 6283185 * SCALE_FP / 1000000
PI = TWO_PI / 2


#
# Python functions
#
def intensity(wave_fn):
    # total intensity ~ (total wave_fn)^2
    intensity = wave_fn**2
    return intensity


def wave_sum(wave_1_fn, wave_2_fn):
    wave_fn = wave_1_fn + wave_2_fn
    return wave_fn


def theta_shifter(theta):
    # shifts theta so it is in range -pi <= theta <= +pi
    theta_abs = abs(theta)
    if theta_abs >= np.pi:
        cycles_to_shift = 1 + ((theta_abs - np.pi) / (2 * np.pi)) // 1
        theta_abs_shifted = theta_abs - (2 * np.pi * cycles_to_shift)
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


def wave_function(t, x, y, params, wave):
    n = params["n"]
    k = params["k"]
    omega = params["omega"]
    decay_exp = params["decay_exp"]

    x0 = wave["x0"]
    y0 = wave["y0"]
    phi = wave["phi"]

    # distance from source
    r = ((x - x0) ** 2 + (y - y0) ** 2) ** 0.5
    # wave function argument
    theta = k * r - omega * t + phi
    # must have -pi <= theta <= +pi for cosine approx., so shift theta as needed
    theta_shifted = theta_shifter(theta)
    # wave function undecayed
    wave_fn = cosine_n_terms(theta_shifted, n)
    # decayed wave function
    # decayed_wave_fn = wave_fn * r**decay_exp
    # return decayed_wave_fn
    return wave_fn


#
# Physical parameters
#
# wave speed
v = 3 * SCALE_FP / 10
# phase shifts of sources
phi_1 = 0 * SCALE_FP
phi_2 = 0 * SCALE_FP
# decay exponent (power of r to show wave dissipation)
# for spherical wave: (ideal decay ->) -1 <= decay_exp <= 0 (<- no decay)
# but -1 is too strong for graphic display
decay_exp = 0
# number of terms in cosine approximation
n = 5
# time; consider only t=0 for now
t = 0.0 * SCALE_FP


#
# Input parameters
#
# frequency
f = 509 * SCALE_FP / 1000
# source separation
d = 2 * SCALE_FP


#
# Wave parameters
#
# angular frequency
omega = 2 * (PI / SCALE_FP) * f
# wave number
k = (omega / v) * SCALE_FP

#
# Plot parameters
#
# number of points to plot along each axis
num_pts = 10
# min and max values for axes
x_min = 0 * SCALE_FP
x_max = 10 * SCALE_FP
y_min = -5 * SCALE_FP
y_max = 5 * SCALE_FP

# np.linspace(min, max, num) returns num evenly spaced numbers over interval min,max
# creates row array
x_s = np.linspace(x_min, x_max, num_pts)
# creates column array if '.reshape(-1,1)'
y_s = np.linspace(y_min, y_max, num_pts).reshape(-1, 1)
# need to create empty arrays to use in nested loop below
wave_1_fn_s = np.empty((num_pts, num_pts))
wave_2_fn_s = np.empty((num_pts, num_pts))


# Source positions
x0_1 = 0.0
y0_1 = d / 2.0
x0_2 = 0.0
y0_2 = -y0_1

# dicts for wave parameters
common = {"omega": omega / SCALE_FP, "k": k / SCALE_FP, "decay_exp": decay_exp, "n": n}
wave_1 = {"x0": x0_1 / SCALE_FP, "y0": y0_1 / SCALE_FP, "phi": phi_1 / SCALE_FP}
wave_2 = {"x0": x0_2 / SCALE_FP, "y0": y0_2 / SCALE_FP, "phi": phi_2 / SCALE_FP}


def intensity_s(f, d):
    # wave function arrays filled w/loops because cannot use ">" w arrays in functions
    for p in range(0, num_pts):
        x = x_s[p] / SCALE_FP
        for q in range(0, num_pts):
            y = y_s[q] / SCALE_FP
            wave_1_fn_s[q, p] = wave_function(t / SCALE_FP, x, y, common, wave_1)
            wave_2_fn_s[q, p] = wave_function(t / SCALE_FP, x, y, common, wave_2)
            # Entire arrays filled by single calls
            # total wave function array due to superposition
            wave_fn_s = wave_sum(wave_1_fn_s, wave_2_fn_s)
            # total intensities array
            intensity_s = intensity(wave_fn_s)
    return intensity_s


#
# Pytest
#

# The path to the contract source code.
CONTRACT_FILE = os.path.join(os.path.dirname(__file__), "math_testing.cairo")


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
    num_tests = 17
    for i in range(num_tests):
        print()  # grab a newline here

        theta = int((i * PI / 4) - (2 * PI))
        theta_no_fp = theta / SCALE_FP

        print(f"> theta       = {theta}")
        print(f"> theta_no_fp = {theta_no_fp}")
        # add '>' before our print messages to indicate they are our messages

        # Call theta_shifter() with theta and print out return value
        ret = await contract.theta_shifter(theta=theta).call()
        if ret.result.theta_shifted >= HALF_PRIME:
            result = ret.result.theta_shifted - PRIME
        else:
            result = ret.result.theta_shifted
        print(f">      theta_shifter(theta) returns: {result}")

        theta_shifter_py_no_fp = theta_shifter(theta_no_fp)
        print(f">   theta_shifter_py(theta) returns: {theta_shifter_py_no_fp}")

        theta_shifter_py = theta_shifter_py_no_fp * SCALE_FP
        if theta_shifter_py >= HALF_PRIME:
            result_py = theta_shifter_py - PRIME
        else:
            result_py = theta_shifter_py
        print(f">   theta_shifter_py(theta) returns: {result_py}")

        # Call cosine_8th() with theta and print out return value
        ret = await contract.cosine_8th(theta=theta).call()
        if ret.result.value >= HALF_PRIME:
            result = ret.result.value - PRIME
        else:
            result = ret.result.value
        print(f">          cosine_8th(theta) returns: {result}")

        cos_8th_py = cosine_n_terms(theta_no_fp, n)
        # if cos_8th_py >= HALF_PRIME:
        #    result_py = cos_8th_py - PRIME
        # else:
        #    result_py = cos_8th_py
        # print(f">   cosine_n_terms(theta, 5) returns: {result_py}")
        print(f">   cosine_n_terms(theta, 5) returns: {cos_8th_py}")

    # print(f">      intensity_plot returns: {result}")
