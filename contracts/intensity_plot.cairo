%lang starknet

from starkware.cairo.common.alloc import alloc
from starkware.cairo.common.cairo_builtins import HashBuiltin
from starkware.cairo.common.math_cmp import is_le, is_nn

from contracts.constants import SCALE_FP, TWO_PI_fp, v_fp, phi_1_fp, phi_2_fp, decay_exp, x_min_fp, x_max_fp, y_min_fp, y_max_fp
from contracts.math import mul_fp, div_fp, div_fp_nfp
from contracts.structs import Common_params, Indiv_params, Two_waves_params
from contracts.wave_physics import intensity_fp, wave_sum_fp, wave_function_fp

//
// Functions to fill arrays
//

// Inner loop goes through all y values, for a particular x value, 
// and calculates combined wave intensity at each position (x, y)
func intensity_fp_s_filler_inner_loop{range_check_ptr}(
    t_fp: felt, intensity_fp_s: felt*, x_fp_s: felt*, y_fp_s: felt*, num_pts_y: felt, waves: Two_waves_params
) {
    alloc_locals;

    // Return after cuing up all y values
    if (num_pts_y == 0) {
        return ();
    }

    // Recursively call to go through all y values
    intensity_fp_s_filler_inner_loop(t_fp, intensity_fp_s + 1, x_fp_s, y_fp_s + 1, num_pts_y - 1, waves);

    // After first return from recursive calls
    // calculate both wave functions at time t_fp and coordinates (x_fp_s[0], y_fp_s[0])
    let wave_fn_1_fp = wave_function_fp(t_fp, x_fp_s[0], y_fp_s[0], waves.common, waves.wave_1);
    let wave_fn_2_fp = wave_function_fp(t_fp, x_fp_s[0], y_fp_s[0], waves.common, waves.wave_2);

    // Add wave functions, find intensity, 
    // and fill member of intensity_fp_s array corresponding to (x, y)
    assert intensity_fp_s[0] = intensity_fp(wave_sum_fp(wave_fn_1_fp, wave_fn_2_fp));

    return ();
}

// Outer loop goes through x values, calling inner loop at each x value
func intensity_fp_s_filler_outer_loop{range_check_ptr}(
    t_fp: felt,
    intensity_fp_s: felt*,
    x_fp_s: felt*,
    y_fp_s: felt*,
    num_pts_x: felt,
    num_pts_y: felt,
    waves: Two_waves_params,
) {
    alloc_locals;

    // Return after going through all x values
    if (num_pts_x == 0) {
        return ();
    }

    // For current x value, call inner loop to go through y all values
    intensity_fp_s_filler_inner_loop(t_fp, intensity_fp_s, x_fp_s, y_fp_s, num_pts_y, waves);

    // Go to next x value
    intensity_fp_s_filler_outer_loop(
        t_fp, intensity_fp_s + num_pts_y, x_fp_s + 1, y_fp_s, num_pts_x - 1, num_pts_y, waves
    );

    return ();
}

// Calculate coordinate values and fill coordinate array
func coordinate_fp_s_filler{}(c_fp_s: felt*, c_min_fp: felt, delta_c_fp: felt, num_pts: felt) {

    // Return after cuing up all coordinate values
    if (num_pts == 0) {
        return ();
    }

    // Recursively call to go through all coordinate values
    coordinate_fp_s_filler(c_fp_s + 1, c_min_fp + delta_c_fp, delta_c_fp, num_pts - 1);

    // after return from 0, now c_min = c_max, num_pts = 1, begin to fill array
    assert c_fp_s[0] = c_min_fp;

    return ();
}

//
//  View function for input of num_pts, lambda, and d; then create intensity plot data
//
@view
func intensity_plot_arr{range_check_ptr}(num_pts: felt, lambda: felt, d: felt) -> (
    intensity_fp_s_len: felt, intensity_fp_s: felt*
) {
    alloc_locals;

    // Check inputs
    with_attr error_message("Check that 2 <=num_pts <= 25") {
        tempvar num_pts_minus_2 = num_pts - 2;
        assert is_nn(num_pts_minus_2) = 1;
    }
    with_attr error_message("Check that 2 <=num_pts <= 25") {
        assert is_le(num_pts, 25) = 1; 
    }
    with_attr error_message("Check that lambda >= 1") {
        tempvar lambda_minus_1 = lambda - 1;
        assert is_nn(lambda_minus_1) = 1;
    }
    with_attr error_message("Check that d >= 0") {
        assert is_nn(d) = 1;
    }

    // Scale up inputs lambda and d to be fixed point values
    local lambda_fp = lambda * SCALE_FP;
    local d_fp = d * SCALE_FP;

    // Allocate memory segments arrays
    let (x_fp_s: felt*) = alloc();
    let (y_fp_s: felt*) = alloc();
    let (intensity_fp_s: felt*) = alloc();

    // Plot size and point spacing
    let x_width_fp = x_max_fp - x_min_fp;
    let delta_x_fp = div_fp_nfp(x_width_fp, num_pts - 1);
    let y_height_fp = y_max_fp - y_min_fp;
    let delta_y_fp = div_fp_nfp(y_height_fp, num_pts - 1);

    // Fill x_fp_s and y_fp_s coordinate value arrays
    coordinate_fp_s_filler(x_fp_s, x_min_fp, delta_x_fp, num_pts);
    coordinate_fp_s_filler(y_fp_s, y_min_fp, delta_y_fp, num_pts);

    // Source positions
    let x0_1_fp = 0;
    let y0_1_fp = div_fp_nfp(d_fp, 2);
    let x0_2_fp = 0;
    let y0_2_fp = -y0_1_fp;

    //
    // Wave parameters needed for plot
    //
    // time (for now, consider only t_fp = 0)
    local t_fp = 0;
    // wave number
    local k_fp = div_fp(TWO_PI_fp, lambda_fp);
    // angular frequency
    local omega_fp = mul_fp(v_fp, k_fp);


    // Struct for waves' common physical parameters
    let common = Common_params(omega_fp=omega_fp, k_fp=k_fp, decay_exp=decay_exp);

    // Struct for each wave's individual parameters
    let wave_1 = Indiv_params(x0_fp=x0_1_fp, y0_fp=y0_1_fp, phi_fp=phi_1_fp);
    let wave_2 = Indiv_params(x0_fp=x0_2_fp, y0_fp=y0_2_fp, phi_fp=phi_2_fp);

    // Struct for two waves
    let waves = Two_waves_params(common=common, wave_1=wave_1, wave_2=wave_2);

    // Fill intensity_fp_s array
    intensity_fp_s_filler_outer_loop(t_fp, intensity_fp_s, x_fp_s, y_fp_s, num_pts, num_pts, waves);

    // Calculate intensity_fp_s_len
    let intensity_fp_s_len = num_pts * num_pts;

    return (intensity_fp_s_len=intensity_fp_s_len, intensity_fp_s=intensity_fp_s);
}
