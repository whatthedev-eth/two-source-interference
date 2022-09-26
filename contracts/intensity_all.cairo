%lang starknet

from starkware.cairo.common.alloc import alloc
from starkware.cairo.common.cairo_builtins import HashBuiltin
from starkware.cairo.common.math import abs_value, signed_div_rem, sqrt
from starkware.cairo.common.math_cmp import is_le


//
// Fixed point math constants
//
const SCALE_FP = 10**20;
const SCALE_FP_SQRT = 10**10;
const RANGE_CHECK_BOUND = 2 ** 120;
// these are not used yet
const PRIME = 3618502788666131213697322783095070105623107215331596699973092056135872020481;
const HALF_PRIME = 1809251394333065606848661391547535052811553607665798349986546028067936010240;

//
// Math constants
//
const TWO_PI = 6283185 * SCALE_FP / 1000000;
const PI = TWO_PI / 2;

//
// Physical parameters
//
// wave speed
const v = 3 * SCALE_FP / 10;
// phase shifts of sources
const phi_1 = 0 * SCALE_FP;
const phi_2 = 0 * SCALE_FP;
// decay exponent (power of r to show wave dissipation)
// for spherical wave: (ideal decay ->) -1 <= decay_exp <= 0 (<- no decay)
// but -1 is too strong for graphic display
const decay_exp = 0;

//
// Plot parameters
//
// number of points to plot along each axis
const num_pts = 10;
// min and max values for axes
const x_min = 0 * SCALE_FP;
const x_max = 10 * SCALE_FP;
const y_min = -5 * SCALE_FP;
const y_max = 5 * SCALE_FP;


//
//  Structs
//

// Common wave parameters
struct Common_params {
    omega: felt,
    k: felt,
    decay_exp: felt,
}

// Individual wave's parameters
struct Indiv_params {
    x0: felt,
    y0: felt,
    phi: felt,
}

// Parameters for two waves
struct Two_waves_params {
    common: Common_params,
    wave_1: Indiv_params,
    wave_2: Indiv_params,
}


//
//  Math functions
//

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
func cosine_8th{range_check_ptr}(theta: felt) -> (value:felt) {

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

    return (value);
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


//
// Wave physics functions
//

func intensity{range_check_ptr}(wave_fn: felt) -> felt {
    let intensity = mul_fp(wave_fn, wave_fn);
    return intensity;
}

func wave_sum{range_check_ptr}(wave_fn_1: felt, wave_fn_2: felt) -> felt {
    return wave_fn_1 + wave_fn_2;
}

func wave_function{range_check_ptr}(t: felt, x: felt, y: felt, common: Common_params, wave: Indiv_params) -> felt {
    // unpack common wave parameters
    let omega = common.omega;
    let k = common.k;
    let decay_exp = common.decay_exp;

    // unpack individual wave parameters
    let x0 = wave.x0;
    let y0 = wave.y0;
    let phi = wave.phi;

    // distance from source
    let r = distance_two_points(x0, y0, x, y);
    
    // wave function argument
    let k_r = mul_fp(k, r);
    let omega_t = mul_fp(omega, t);
    let theta = k_r - omega_t + phi;
                              
    // wave function undecayed
    let wave_fn = cosine_8th(theta);
    
    // decayed wave function
    //decayed_wave_fn = wave_fn * r**decay_exp;
    
    //return decayed_wave_fn;
    return wave_fn;
}


//
// Functions to fill arrays
//
func intensity_s_filler{range_check_ptr}(t: felt, intensity_s: felt*, x_s: felt*, y_s: felt*, num_pts_x: felt, num_pts_y: felt, num_pts_y_const: felt, waves: Two_waves_params) { 
    alloc_locals; 

    if (num_pts_x == 0) {
         return ();
    }
    if (num_pts_y == 0) {
         return ();
    }
    
    // otherwise, call recursively for x
    intensity_s_filler(t, intensity_s + 1, x_s + 1, y_s, num_pts_x - 1, num_pts_y, num_pts_y_const, waves);
        
    // begin at intensity_s[num_pts_x - 1], x_s[num_pts_x - 1], num_pts_x = 1 (after return from num_pts_x = 0)
    tempvar x = x_s[0];
        
    // NESTED recursive call for y, beginning with num_pts_x = 1
    intensity_s_filler(t, intensity_s + num_pts, x_s, y_s + 1, num_pts_x, num_pts_y - 1, num_pts_y_const, waves);
        
    // begin at intensity_s[num_pts_x * num_pts_y - 1], y_s[num_pts_y - 1], num_pts_y = 1
    tempvar y = y_s[0];
            
    let wave_fn_1 = wave_function(t, x, y, waves.common, waves.wave_1);
    let wave_fn_2 = wave_function(t, x, y, waves.common, waves.wave_2);

    assert intensity_s[0] = intensity(wave_sum(wave_fn_1, wave_fn_2));
    
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
func intensity_plot{range_check_ptr}(f: felt, d: felt) -> () {
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
    intensity_s_filler{}(t, intensity_s, x_s, y_s, num_pts, num_pts, num_pts, waves);

    return();
}