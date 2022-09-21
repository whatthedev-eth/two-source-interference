%lang starknet

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
