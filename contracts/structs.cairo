%lang starknet

// Common wave parameters
struct Common_params {
    omega_fp: felt,
    k_fp: felt,
    decay_exp: felt,
}

// Individual wave's parameters
struct Indiv_params {
    x0_fp: felt,
    y0_fp: felt,
    phi_fp: felt,
}

// Parameters for two waves
struct Two_waves_params {
    common: Common_params,
    wave_1: Indiv_params,
    wave_2: Indiv_params,
}
