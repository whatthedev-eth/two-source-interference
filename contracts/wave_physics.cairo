%lang starknet

from starkware.cairo.common.cairo_builtins import HashBuiltin
from contracts.math import cosine_8th_fp, distance_two_points_fp, mul_fp
from contracts.structs import Common_params, Indiv_params, Two_waves_params

//
// Functions for wave physics
//

// Calculates wave intensity ~ (wave function)^2
func intensity_fp{range_check_ptr}(wave_fn_fp: felt) -> felt {
    let intensity_fp = mul_fp(wave_fn_fp, wave_fn_fp);
    return intensity_fp;
}

// Adds two wave functions
func wave_sum_fp{range_check_ptr}(wave_fn_1_fp: felt, wave_fn_2_fp: felt) -> felt {
    return wave_fn_1_fp + wave_fn_2_fp;
}

// Calculates wave function at time t and coordinates (x, y)
func wave_function_fp{range_check_ptr}(
    t_fp: felt, x_fp: felt, y_fp: felt, common: Common_params, wave: Indiv_params
) -> felt {
    // Unpack common wave parameters
    let omega_fp = common.omega_fp;
    let k_fp = common.k_fp;
    let decay_exp = common.decay_exp;

    // Unpack individual wave parameters
    let x0_fp = wave.x0_fp;
    let y0_fp = wave.y0_fp;
    let phi_fp = wave.phi_fp;

    // Distance from source
    let r_fp = distance_two_points_fp(x0_fp, y0_fp, x_fp, y_fp);

    //
    // wave function = amplitude*cos(k*r - omega*t + phi)
    // (amplitude = 1)
    //
    // Wave function argument
    let k_r_fp = mul_fp(k_fp, r_fp);
    let omega_t_fp = mul_fp(omega_fp, t_fp);
    let theta_fp = k_r_fp - omega_t_fp + phi_fp;

    // Wave function undecayed
    let (wave_fn_fp) = cosine_8th_fp(theta_fp);

    // TODO Decayed wave function
    // TODO Write func power_fp in math.cairo
    // r_to_power_decay_exp_fp = power_fp(r_fp, decay_exp);
    // decayed_wave_fn_fp = mul_fp(wave_fn_fp, r_to_power_decay_exp_fp);
    // return decayed_wave_fn_fp;

    return wave_fn_fp;
}
