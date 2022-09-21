%lang starknet

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
const num_pts = 20;
// min and max values for axes
const x_min = 0 * SCALE_FP;
const x_max = 10 * SCALE_FP;
const y_min = -5 * SCALE_FP;
const y_max = 5 * SCALE_FP;