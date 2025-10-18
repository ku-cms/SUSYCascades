#include "CascadesTreeEncoder.hh"
uint16_t CascadesTreeEncoder::Encode(const int& A, const int& B, const int& C, const int& D, const int& E) {
  if (A < 1 || A > 8 || B < 1 || B > 3 || C < 1 || C > 3 || D < 0 || D > 5 || E < 0 || E > 5) return 0;
  else return ((A - 1) << 11) | ((B - 1) << 9) | ((C - 1) << 7) | (D << 4) | E;
}

std::tuple<int, int, int, int, int> CascadesTreeEncoder::Decode(const uint16_t& packed) {
    int E =  packed        & 0b111;       // 3 bits for E
    int D = (packed >> 3)  & 0b111;       // 3 bits for D
    int C = ((packed >> 6) & 0b11) + 1;   // 2 bits for C
    int B = ((packed >> 8) & 0b11) + 1;   // 2 bits for B
    int A = ((packed >> 10) & 0b111) + 1; // 3 bits for A
    return std::make_tuple(A, B, C, D, E);
}
