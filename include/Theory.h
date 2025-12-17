#pragma once
#include <array>
#include <cstddef>

enum class Notes { C, Cs, D, Ds, E, F, Fs, G, Gs, A, As, B };

enum class Mode { Major, Minor };

struct Note {
  Notes note;
  int octave;
};

struct Scale {
  Notes root;
  Mode mode;
};

using ScaleIntervals = std::array<int, 7>;

inline const std::array<ScaleIntervals, 2> SCALE_INTERVALS = {{
    {0, 2, 4, 5, 7, 9, 11}, // major
    {0, 2, 3, 5, 7, 8, 10}  // minor
}};

inline const auto& getScaleIntervals(Mode mode) {
    return SCALE_INTERVALS[static_cast<size_t>(mode)];
}
