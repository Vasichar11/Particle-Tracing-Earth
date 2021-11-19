#include <random>
#include <array>
#include <functional>
#include <iostream>
#include <algorithm>
#include <iostream>
#include "common.h"
#include "constants.h"

//"Symmetrical distribution" by mirroring half normal and shifting it.
 std::array<real, Constants::aeq_dstr> symmetrical_dstr(real halfmean, real std, real shift);

