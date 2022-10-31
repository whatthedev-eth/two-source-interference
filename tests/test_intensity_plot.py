from re import S
import pytest
import os
from starkware.starknet.testing.starknet import Starknet
import asyncio
import numpy as np
import json

#
# Constants from constants.cairo contract
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
# Input parameters
#
# number of points to plot along each axis
num_pts = 25
# wavelength (cannot use `lambda` in python) scaled up by 100 for Cairo input
lambda1 = 150
# source separation scaled up by 100 for Cairo input
d = 500
# scale up lambda1 and d to simulate what happens in Cairo code
lambda_fp = int(lambda1 * SCALE_FP / 100)
d_fp = int(d * SCALE_FP / 100)

#
# Math constants
#
# Fixed point values
TWO_PI_fp = 6283185 * SCALE_FP / 1000000
PI_fp = TWO_PI_fp / 2
# Non fixed point values
TWO_PI = TWO_PI_fp / SCALE_FP
PI = PI_fp / SCALE_FP


#
# Physical parameters
#
# wave speed
v_fp = 3 * SCALE_FP / 10
# phase shifts of sources
phi_1_fp = 0 * SCALE_FP
phi_2_fp = 0 * SCALE_FP
# decay exponent (power of r to show wave dissipation)
# for spherical wave: (ideal decay ->) -1 <= decay_exp <= 0 (<- no decay)
# but -1 is too strong for graphic display
decay_exp = 0
# number of terms in cosine approximation
n = 5
# time; consider only t=0 for now
t_fp = 0.0 * SCALE_FP


#
# Calculated wave parameters
#
# wave number
k_fp = (TWO_PI_fp / lambda_fp) * SCALE_FP
# angular frequency
omega_fp = (v_fp / SCALE_FP) * lambda_fp


#
# Plot parameters
#
# min and max values for axes
x_min_fp = 0 * SCALE_FP
x_max_fp = 10 * SCALE_FP
y_min_fp = -5 * SCALE_FP
y_max_fp = 5 * SCALE_FP

# np.linspace(min, max, num) returns num evenly spaced numbers over interval min,max
# creates row array
x_s_fp = np.linspace(x_min_fp, x_max_fp, num_pts)
# creates column array if '.reshape(-1,1)'
y_s_fp = np.linspace(y_min_fp, y_max_fp, num_pts).reshape(-1, 1)
# need to create empty arrays to use in nested loop below
wave_1_fn_s = np.empty((num_pts, num_pts))
wave_2_fn_s = np.empty((num_pts, num_pts))

# Source positions
x0_1_fp = 0.0
y0_1_fp = d_fp / 2.0
x0_2_fp = 0.0
y0_2_fp = -d_fp / 2.0


# dicts for wave parameters
common = {
    "omega": omega_fp / SCALE_FP,
    "k": k_fp / SCALE_FP,
    "decay_exp": decay_exp,
    "n": n,
}
wave_1 = {
    "x0": x0_1_fp / SCALE_FP,
    "y0": y0_1_fp / SCALE_FP,
    "phi": phi_1_fp / SCALE_FP,
}
wave_2 = {
    "x0": x0_2_fp / SCALE_FP,
    "y0": y0_2_fp / SCALE_FP,
    "phi": phi_2_fp / SCALE_FP,
}


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
    # using PI_no_fp instead of np.pi to match cairo file
    theta_abs = abs(theta)
    if theta_abs >= PI:
        cycles_to_shift = 1 + ((theta_abs - PI) / TWO_PI) // 1
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
    # 2(n-1) = order
    # cosine(theta) ~= ((-1)^n)*(theta^(2n))/(2n)!
    #               ~= 1 - theta^2/2! + theta^4/4! - theta^6/6! + ...
    cos_n_terms = 0
    # must have -pi <= theta <= +pi for cosine approx., so shift theta as needed
    theta_shifted = theta_shifter(theta)
    for i in range(n):
        power_neg_one = (-1) ** i
        power_theta = theta_shifted ** (2 * i)
        fact = np.math.factorial(2 * i)
        cos_n_terms += power_neg_one * power_theta / fact
    return cos_n_terms


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


def intensity_plot_arr(num_pts, lambda1, d):
    # wave function arrays filled via loops because cannot use ">" w arrays in functions
    for p in range(0, num_pts):
        x = x_s_fp[p] / SCALE_FP
        for q in range(0, num_pts):
            y = y_s_fp[q] / SCALE_FP
            wave_1_fn_s[q, p] = wave_function(t_fp / SCALE_FP, x, y, common, wave_1)
            wave_2_fn_s[q, p] = wave_function(t_fp / SCALE_FP, x, y, common, wave_2)
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
CONTRACT_FILE = os.path.join("contracts", "intensity_plot.cairo")


@pytest.mark.asyncio
async def test():

    starknet = await Starknet.empty()
    contract = await starknet.deploy(
        source=CONTRACT_FILE,
    )
    print()  # print blank line

    # Cairo intensity array (input lambda, d that were already scaled up by 100)
    ret = await contract.intensity_plot_arr(num_pts, lambda1, d).call()
    # dump to json file
    with open("tests/test_intensity_plot_cairo.json", "w") as outfile:
        json.dump(ret.result, outfile)

    # Python intensity array (input unscaled lambda, d)
    intensity_s = intensity_plot_arr(num_pts, lambda_fp / SCALE_FP, d_fp / SCALE_FP)
    # dump to a different json file
    with open("tests/test_intensity_plot_python.json", "w") as outfile:
        json.dump(intensity_s.tolist(), outfile)

    print(f"> intensities for num_pts={num_pts}, lambda={lambda1}, d={d}) returns:")
    print()
    print(f"> intensity_plot_arr from cairo     member")
    print(f"> intensity_s_arr from python")

    # Print array members one line at a time:
    # Cairo calculation above Python calculation, with array member # to the right
    for p in range(0, num_pts):
        for q in range(0, num_pts):
            print()

            intensity_cairo = ret.result[0][q + p * num_pts]
            print(f"> {intensity_cairo}                 {q + p * num_pts}")

            intensity_py = int(intensity_s[q, p] * SCALE_FP)
            print(f"> {intensity_py}")

    # print ret.call_info.execution_resources to get n_steps
    print()
    print(ret.call_info.execution_resources)
