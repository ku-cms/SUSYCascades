#include "CascadesTreeEncoder.hh"
uint16_t CascadesTreeEncoder::Encode(const int& A, const int& B, const int& C, const int& D, const int& E) {
  if (A < 1 || A > 8 || B < 1 || B > 3 || C < 1 || C > 3 || D < 0 || D > 3 || E < 0 || E > 3) return 0;
  else return ((A - 1) << 8) | ((B - 1) << 6) | ((C - 1) << 4) | (D << 2) | E;
}

std::tuple<int, int, int, int, int> CascadesTreeEncoder::Decode(const uint16_t& packed) {
    int E =  packed        & 0b11;
    int D = (packed >> 2)  & 0b11;
    int C = ((packed >> 4) & 0b11) + 1;
    int B = ((packed >> 6) & 0b11) + 1;
    int A = ((packed >> 8) & 0b111) + 1;
    return std::make_tuple(A, B, C, D, E);
}

