#include "constants.h"

extern "C" {
    void krnl(uint32_t* input, uint8_t* output) {
        #pragma HLS INTERFACE m_axi port = input bundle = gmem0
        #pragma HLS INTERFACE m_axi port = output bundle = gmem1

        /*Insert Code Here*/
    }
}