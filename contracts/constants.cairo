%lang starknet

//
// Fixed point math constants
//
const SCALE_FP = 10**20;
const SCALE_FP_SQRT = 10**10;
const RANGE_CHECK_BOUND = 2 ** 120;
// Constants for felts, not used yet
// const PRIME = 3618502788666131213697322783095070105623107215331596699973092056135872020481;
// const HALF_PRIME = 1809251394333065606848661391547535052811553607665798349986546028067936010240;

//
// Math constants
//
const TWO_PI_fp = 628318530 * SCALE_FP / 100000000;
const PI_fp = 314159265 * SCALE_FP / 100000000;

//
// Physical parameters
//
// Wave speed
const v_fp = 30 * SCALE_FP;
// Phase shifts of sources
const phi_1_fp = 0 * SCALE_FP;
const phi_2_fp = 0 * SCALE_FP;

// Decay exponent (power of r to show wave dissipation) for spherical wave: 
// (ideal decay ->) -1 <= decay_exp <= 0 (<- no decay)
// but -1 is too strong for graphic display
const decay_exp = 0;

//
// Plot parameters
//
// Min and max values for axes
const x_min_fp = 0 * SCALE_FP;
const x_max_fp = 1000 * SCALE_FP;
const y_min_fp = -500 * SCALE_FP;
const y_max_fp = 500 * SCALE_FP;