%lang starknet

from starkware.cairo.common.cairo_builtins import HashBuiltin
from starkware.cairo.common.alloc import alloc
from starkware.cairo.common.math import abs_value, signed_div_rem, sqrt
from starkware.cairo.common.math_cmp import is_le

from contracts.constants import SCALE_FP, SCALE_FP_SQRT, RANGE_CHECK_BOUND, TWO_PI, PI


func sqrt_fp{range_check_ptr}(x: felt) -> felt {
    let x_ = sqrt(x); // notice: sqrt() now returns a single felt, not a tuple anymore (tuple is returned for cairo < 0.10)
    let y = x_ * SCALE_FP_SQRT;  // compensate for the square root
    return y;
}

func mul_fp{range_check_ptr}(a: felt, b: felt) -> felt {
    // signed_div_rem by SCALE_FP after multiplication
    tempvar product = a * b;
    let (c, _) = signed_div_rem(product, SCALE_FP, RANGE_CHECK_BOUND);
    return c;
}

func div_fp{range_check_ptr}(a: felt, b: felt) -> felt {
    // multiply by SCALE_FP before signed_div_rem
    tempvar a_scaled = a * SCALE_FP;
    let (c, _) = signed_div_rem(a_scaled, b, RANGE_CHECK_BOUND);
    return c;
}

// func mul_fp_ul{range_check_ptr}(a: felt, b_ul: felt) -> felt {
//     let c = a * b_ul;
//     return c;
// }

func div_fp_ul{range_check_ptr}(a: felt, b_ul: felt) -> felt {
    let (c, _) = signed_div_rem(a, b_ul, RANGE_CHECK_BOUND);
    return c;
}

func distance_two_points{range_check_ptr}(x0: felt, y0: felt, x: felt, y: felt) -> felt {
    let x_diff = x - x0;
    let y_diff = y - y0;
    let x_diff_sq = mul_fp(x_diff, x_diff);
    let y_diff_sq = mul_fp(y_diff, y_diff);
    let sum = x_diff_sq + y_diff_sq;
    let r = sqrt_fp(sum);
    return r;
}

@view
func cosine_8th{range_check_ptr}(theta: felt) -> felt {

    //
    // cos(theta) ~= 1 - theta^2/2! + theta^4/4! - theta^6/6! + theta^8/8!
    //

    let theta_shifted = theta_shifter(theta);

    let theta_2 = mul_fp(theta_shifted, theta_shifted);
    let theta_4 = mul_fp(theta_2, theta_2);
    let theta_6 = mul_fp(theta_2, theta_4);
    let theta_8 = mul_fp(theta_2, theta_6);

    let theta_2_div2 = div_fp_ul(theta_2, 2);
    let theta_4_div24 = div_fp_ul(theta_4, 24);
    let theta_6_div720 = div_fp_ul(theta_6, 720);
    let theta_8_div40320 = div_fp_ul(theta_8, 40320);

    let value = 1*SCALE_FP - theta_2_div2 + theta_4_div24 - theta_6_div720 + theta_8_div40320;

    return value;
}

@view
func theta_shifter{range_check_ptr}(theta: felt) -> felt {
    alloc_locals;

    //
    // shifts theta so it is in range -pi <= theta <= +pi
    //

    local range_check_ptr = range_check_ptr;
    local theta_abs = abs_value(theta);

    // if theta_abs >= pi, then shift to range theta_abs < +pi
    let bool_size = is_le(PI, theta_abs);
    if (bool_size == 1) {
        //
        // # of cycles to shift theta by = 1 + ((theta_abs - pi)/(2*pi)) // 1
        //
        let diff = theta_abs - PI;
        // should I use unsigned_div_rem here bc I know it's > 0??????????????????????????????
        let (cycles_minus_one, _) = signed_div_rem(diff, TWO_PI, RANGE_CHECK_BOUND);
        let cycles = 1 + cycles_minus_one;
        let shift = cycles * TWO_PI;
        tempvar theta_abs_shifted = theta_abs - shift;
        tempvar range_check_ptr = range_check_ptr;
    } else {
        tempvar theta_abs_shifted = theta_abs;
        tempvar range_check_ptr = range_check_ptr;
    }

    if (theta_abs == theta) {
        tempvar theta_shifted = theta_abs_shifted;
        tempvar range_check_ptr = range_check_ptr;
    } else {
        tempvar theta_shifted = -theta_abs_shifted;
        tempvar range_check_ptr = range_check_ptr;
    }

    return theta_shifted;
}