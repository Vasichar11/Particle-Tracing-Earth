#include "is_in_packet.h"

// Check if the particle is within the wave's latitude.
// Return wave's index where there's matching latitude (WPI) and the boolean value wpi.
int is_in_packet(const real min_lat, const real max_lat, const real latitude_tmp, const int i, std::vector<real>& wave_lat)
{
    std::vector<real>::iterator ptr, requested_index;  // To find requested index
    int index = -1;  // Default index if no match is found
    real diff, min_diff = 100;  // To find minimum difference between latitudes

    if (min_lat < (latitude_tmp * Universal::R2D) && (latitude_tmp * Universal::R2D) < max_lat)
    {
        for (ptr = wave_lat.begin() + i; ptr < wave_lat.begin() + (Simulation::puls_dur + i); ptr++)
        {
            diff = std::abs(*ptr - latitude_tmp * Universal::R2D);  // Return absolute difference
            if (diff < min_diff)
            {
                requested_index = ptr;  // Current requested index
                min_diff = diff;  // Current minimum
            }
        }
        index = std::distance(wave_lat.begin(), requested_index);  // Final requested index
    }

    return index;
}