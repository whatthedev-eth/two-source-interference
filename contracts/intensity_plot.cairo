%lang starknet

from starkware.cairo.common.alloc import alloc
from starkware.cairo.common.cairo_builtins import HashBuiltin

from contracts.constants import TWO_PI, v, phi_1, phi_2, decay_exp, x_min, x_max, y_min, y_max
from contracts.math import mul_fp, div_fp, div_fp_ul
from contracts.structs import Common_params, Indiv_params, Two_waves_params
from contracts.wave_physics import intensity, wave_sum, wave_function


//
// Functions to fill arrays
//

func intensity_s_filler_inner_loop{range_check_ptr}(t: felt, intensity_s: felt*, x_s: felt*, y_s: felt*, num_pts_y: felt, waves: Two_waves_params) { 
    alloc_locals; 

    if (num_pts_y == 0) {
         return ();
    }
    
    intensity_s_filler_inner_loop(t, intensity_s + 1, x_s, y_s + 1, num_pts_y - 1, waves);
    
    tempvar y = y_s[0];   
    tempvar x = x_s[0];

    let wave_fn_1 = wave_function(t, x, y, waves.common, waves.wave_1);
    let wave_fn_2 = wave_function(t, x, y, waves.common, waves.wave_2);

    assert intensity_s[0] = intensity(wave_sum(wave_fn_1, wave_fn_2));
    
    return();
}

func intensity_s_filler_outer_loop{range_check_ptr}(t: felt, intensity_s: felt*, x_s: felt*, y_s: felt*, num_pts_x: felt, num_pts_y: felt, waves: Two_waves_params) { 
    alloc_locals; 

    if (num_pts_x == 0) {
         return ();
    }
    
    intensity_s_filler_inner_loop(t, intensity_s, x_s, y_s, num_pts_y, waves);

    intensity_s_filler_outer_loop(t, intensity_s + num_pts_y, x_s + 1, y_s, num_pts_x - 1, num_pts_y, waves);
    
    return();
}

func coordinate_s_filler{}(c_s: felt*, c_min: felt, delta_c: felt, num_pts: felt) { 

    if (num_pts == 0) {
        // when 0 is reached, return
        return ();
    }

    // otherwise, call recursively
    coordinate_s_filler(c_s + 1, c_min + delta_c, delta_c, num_pts - 1);

    // after return from 0, now c_min = c_max, num_pts = 1, begin to fill array
    assert c_s[0] = c_min;
    
    return();
}


//
//  External function for input of f and d, and then create plot data
//
@external
func intensity_plot_arr{range_check_ptr}(num_pts: felt, f: felt, d: felt) -> (intensity_s_len: felt, intensity_s: felt*) {
    alloc_locals;

    // allocate memory segments arrays
    let (x_s: felt*) = alloc();
    let (y_s: felt*) = alloc();
    let (intensity_s: felt*) = alloc();
    
    // plot size and point spacing
    let x_width = x_max - x_min;
    let delta_x = div_fp_ul(x_width, num_pts - 1);
    let y_height = y_max - y_min;
    let delta_y = div_fp_ul(y_height, num_pts - 1);

    // fill x_s and y_s coordinate value arrays
    coordinate_s_filler(x_s, x_min, delta_x, num_pts);
    coordinate_s_filler(y_s, y_min, delta_y, num_pts);

    // source positions
    let x0_1 = 0;
    let y0_1 = d / 2;
    let x0_2 = 0;
    let y0_2 = -y0_1;

    // wave parameters needed for plot
    // time (for now, consider only t = 0)
    local t = 0;
    // angular frequency
    local omega = mul_fp(TWO_PI, f);
    // wave number
    local k = div_fp(omega, v);

    // struct for waves' common physical parameters
    let common = Common_params(omega=omega, k=k, decay_exp=decay_exp);

    // struct for each wave's individual parameters
    let wave_1 = Indiv_params(x0=x0_1, y0=y0_1, phi=phi_1);
    let wave_2 = Indiv_params(x0=x0_2, y0=y0_2, phi=phi_2);
    
    // struct for two waves
    let waves = Two_waves_params(common=common, wave_1=wave_1, wave_2=wave_2);

    // fill intensity_s array
    intensity_s_filler_outer_loop(t, intensity_s, x_s, y_s, num_pts, num_pts, waves);

    // calculate intensity_s_len
    let intensity_s_len = num_pts * num_pts;

    return(intensity_s_len=intensity_s_len, intensity_s=intensity_s);
}