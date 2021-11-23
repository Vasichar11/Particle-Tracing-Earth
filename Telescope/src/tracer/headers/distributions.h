#include <random>
#include <array>
#include <algorithm>
#include <iostream>
#include "common.h"
#include "constants.h"


//"Half-normal distribution" arround mean with stdeviation
 std::array<real, Constants::aeq_dstr/2> half_norm(real mean, real std);

//"Symmetrical distribution" by mirroring half normal and shifting it.
std::array<real, Constants::aeq_dstr> symmetrical_dstr(std::array<real, Constants::aeq_dstr/2> half_dstr, real shift);


