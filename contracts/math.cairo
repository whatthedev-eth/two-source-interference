%lang starknet

from starkware.cairo.common.cairo_builtins import HashBuiltin
from starkware.cairo.common.alloc import alloc
from starkware.cairo.common.math import abs_value, signed_div_rem, sqrt
from starkware.cairo.common.math_cmp import is_le

from contracts.constants import SCALE_FP, SCALE_FP_SQRT, RANGE_CHECK_BOUND, TWO_PI_fp, PI_fp

// Takes square root of fixed point quantity "x_fp"
func sqrt_fp{range_check_ptr}(x_fp: felt) -> felt {
    let x_ = sqrt(x_fp);  // notice: sqrt() now returns a single felt, not a tuple anymore (tuple is returned for cairo < 0.10)
    let y_fp = x_ * SCALE_FP_SQRT;  // compensate for the square root
    return y_fp;
}

// Multiplies fixed point quantity "a_fp" by fixed point quantity "b_fp", with range check
func mul_fp{range_check_ptr}(a_fp: felt, b_fp: felt) -> felt {
    // signed_div_rem by SCALE_FP after multiplication
    tempvar product = a_fp * b_fp;
    let (c_fp, _) = signed_div_rem(product, SCALE_FP, RANGE_CHECK_BOUND);
    return c_fp;
}

// Divides fixed point quantity "a_fp" by fixed point quantity "b_fp", with range check
func div_fp{range_check_ptr}(a_fp: felt, b_fp: felt) -> felt {
    // multiply by SCALE_FP before signed_div_rem
    tempvar a_scaled = a_fp * SCALE_FP;
    let (c_fp, _) = signed_div_rem(a_scaled, b_fp, RANGE_CHECK_BOUND);
    return c_fp;
}

// Multiplies fixed point quantity "a_fp" by non-fixed point quantity "b"
func mul_fp_nfp{range_check_ptr}(a_fp: felt, b: felt) -> felt {
    let c_fp = a_fp * b;
    return c_fp;
}

// Divides fixed point quantity "a_fp" by non-fixed point quantity "b", with range check
func div_fp_nfp{range_check_ptr}(a_fp: felt, b: felt) -> felt {
    let (c_fp, _) = signed_div_rem(a_fp, b, RANGE_CHECK_BOUND);
    return c_fp;
}

// Finds distance between fixed point coordinate values (x0, y0) and (x, y)
func distance_two_points_fp{range_check_ptr}(x0_fp: felt, y0_fp: felt, x_fp: felt, y_fp: felt) -> felt {
    let x_diff_fp = x_fp - x0_fp;
    let y_diff_fp = y_fp - y0_fp;
    let x_diff_sq_fp = mul_fp(x_diff_fp, x_diff_fp);
    let y_diff_sq_fp = mul_fp(y_diff_fp, y_diff_fp);
    let sum_fp = x_diff_sq_fp + y_diff_sq_fp;
    let r_fp = sqrt_fp(sum_fp);
    return r_fp;
}

// Shifts fixed point theta value in radians to an equivalent value
// within range -pi <= theta_shifted <= +pi
@view
func theta_shifter_fp{range_check_ptr}(theta_fp: felt) -> (value_fp: felt) {
    alloc_locals;

    local range_check_ptr = range_check_ptr;
    local theta_abs_fp = abs_value(theta_fp);

    // if theta_abs >= pi, then shift to range theta_abs < pi
    let bool_size = is_le(PI_fp, theta_abs_fp);
    if (bool_size == 1) {
        //
        // # of cycles to shift theta by = 1 + ((theta_abs - pi)/(2*pi)) // 1
        //
        let diff_fp = theta_abs_fp - PI_fp;
        // Use signed_div_rem (instead of div_fp) next to get non-FP integer result
        let (cycles_minus_one, _) = signed_div_rem(diff_fp, TWO_PI_fp, RANGE_CHECK_BOUND);
        let cycles = 1 + cycles_minus_one;
        let shift_fp = mul_fp_nfp(TWO_PI_fp, cycles);
        tempvar theta_abs_shifted_fp = theta_abs_fp - shift_fp;
        tempvar range_check_ptr = range_check_ptr;
    } else {
        tempvar theta_abs_shifted_fp = theta_abs_fp;
        tempvar range_check_ptr = range_check_ptr;
    }

    if (theta_abs_fp == theta_fp) {
        tempvar theta_shifted_fp = theta_abs_shifted_fp;
        tempvar range_check_ptr = range_check_ptr;
    } else {
        tempvar theta_shifted_fp = -theta_abs_shifted_fp;
        tempvar range_check_ptr = range_check_ptr;
    }

    return (value_fp=theta_shifted_fp);
}

// Approximates cosine(theta) for FP theta value, using 5 terms (to 8th order) of Taylor series
@view
func cosine_8th_fp{range_check_ptr}(theta_fp: felt) -> (value_fp: felt) {
    //
    // cos(theta) ~= 1 - theta^2/2! + theta^4/4! - theta^6/6! + theta^8/8!
    //
    // This approximation works well only if theta is within range of -pi to pi
    // Shift theta as needed to equivalent angle -pi <= theta_shifted <= +pi
    let (theta_shifted_fp) = theta_shifter_fp(theta_fp);

    let theta_2_fp = mul_fp(theta_shifted_fp, theta_shifted_fp);
    let theta_4_fp = mul_fp(theta_2_fp, theta_2_fp);
    let theta_6_fp = mul_fp(theta_2_fp, theta_4_fp);
    let theta_8_fp = mul_fp(theta_2_fp, theta_6_fp);

    let theta_2_div2_fp = div_fp_nfp(theta_2_fp, 2);
    let theta_4_div24_fp = div_fp_nfp(theta_4_fp, 24);
    let theta_6_div720_fp = div_fp_nfp(theta_6_fp, 720);
    let theta_8_div40320_fp = div_fp_nfp(theta_8_fp, 40320);

    let value_fp = 1 * SCALE_FP - theta_2_div2_fp + theta_4_div24_fp - theta_6_div720_fp + theta_8_div40320_fp;

    return (value_fp=value_fp);
}
