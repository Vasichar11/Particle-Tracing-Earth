#include <random>
#include <iostream>
#include <array>
#include <algorithm>
#include "common.h"
#include "constants.h"


//"Half-normal distribution" arround halfmean with stdeviation
 std::array<real, Constants::aeq_dstr/2 + 1> half_norm(real halfmean, real std);

//"Symmetrical distribution" by mirroring half normal and shifting it.
std::array<real, Constants::aeq_dstr> symmetrical_dstr(std::array<real, Constants::aeq_dstr/2 + 1> half_dstr, real shift);


