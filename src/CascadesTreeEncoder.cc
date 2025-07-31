#include "CascadesTreeEncoder.hh"
uint16_t CascadesTreeEncoder::Encode(const int& A, const int& B, const int& C, const int& D, const int& E) {
  if (A < 1 || A > 8 || B < 1 || B > 3 || C < 1 || C > 3 || D < 0 || D > 4 || E < 0 || E > 4) return 0;
  else return ((A - 1) << 10) | ((B - 1) << 8) | ((C - 1) << 6) | (D << 3) | E;
}

std::tuple<int, int, int, int, int> CascadesTreeEncoder::Decode(const uint16_t& packed) {
    int E =  packed        & 0b111; 
    int D = (packed >> 3)  & 0b111;
    int C = ((packed >> 6) & 0b11) + 1;
    int B = ((packed >> 8) & 0b11) + 1;
    int A = ((packed >> 10) & 0b111) + 1;
    return std::make_tuple(A, B, C, D, E);
}

