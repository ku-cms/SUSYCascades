#ifndef CascadesTreeEncoder_HH
#define CascadesTreeEncoder_HH
#include <tuple>
#include <cstdint>

class CascadesTreeEncoder {
  public:
    // Encode cascade gen prod & decay values into a single uint16_t
    // Production: 1-8, SlepSneu1st: 1-3, SlepSneu2nd: 1-3, N2_1st: 0-4, N2_2nd: 0-4
    static uint16_t Encode(const int& A, const int& B, const int& C, const int& D, const int& E);
    // Example encode: uint16_t packed = CascadesTreeEncoder::Encode(5, 2, 2, 1, 0);
    static std::tuple<int, int, int, int, int> Decode(const uint16_t& packed);
    // Example decode: auto [Production, SlepSneu1st, SlepSneu2nd, N2_1st, N2_2nd] = CascadesTreeEncoder::Decode(packed);
    //                 std::cout << "Production mode = " << Production << std::endl;
};

#endif
