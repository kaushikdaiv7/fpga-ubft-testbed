#include "constants.h"

/*Insert included libraries Here*/
//#include "../include_host/blake3.h"
#include <vector> 
#include <random>
#include <assert.h>
#include <string.h>
#include <stdio.h>
#include <ap_int.h>
#include <stdint.h> 

#define CHUNK_START 1 << 0
#define CHUNK_END 1 << 1
#define PARENT 1 << 2
#define ROOT 1 << 3

#define BLAKE3_OUT_LEN 32
#define BLAKE3_BLOCK_LEN 64
#define BLAKE3_CHUNK_LEN 1024


// This struct is private.
typedef struct _blake3_chunk_state {
uint32_t chaining_value[8];
uint64_t chunk_counter;
uint8_t block[BLAKE3_BLOCK_LEN];
uint8_t block_len;
uint8_t blocks_compressed;
uint32_t flags;
} _blake3_chunk_state;

// An incremental hasher that can accept any number of writes.
typedef struct blake3_hasher {
_blake3_chunk_state chunk_state;
uint32_t key_words[8];                                                // check later
uint32_t cv_stack[8 * 54]; // Space for 54 subtree chaining values:
uint8_t cv_stack_len;      // 2^54 * CHUNK_LEN = 2^64
uint32_t flags;
} blake3_hasher;

typedef struct {
    blake3_hasher updated_hasher;
    uint32_t popped_cv[8];
} hasher_pop_result;

typedef struct {
    uint8_t hash[BLAKE3_BLOCK_LEN];
} blake3_hash_output;


blake3_hasher blake3_hasher_init();

blake3_hasher blake3_hasher_update(blake3_hasher self, const void *input,
                          size_t input_len);

blake3_hash_output blake3_hasher_finalize(const blake3_hasher self);

void print_blake3_hasher(blake3_hasher hasher);

extern "C" {

    void krnl(uint8_t* input, uint8_t* output) {
        #pragma HLS INTERFACE m_axi port = input bundle = gmem0
        #pragma HLS INTERFACE m_axi port = output bundle = gmem1

        blake3_hasher hasher = blake3_hasher_init();
        
        // uint8_t data[4] = {0xAB, 0xCD, 0xEF, 0x01};
        for (size_t i = 0; i < 4; i++) {
            printf("%02x", input[i]);
        }
        hasher = blake3_hasher_update(hasher, (const uint8_t *)input, 4); // sizeof(input)
        
        // print_blake3_hasher(hasher);
        printf("\n");
        
        // hasher = blake3_hasher_update(hasher, (const uint8_t *)data1, sizeof(data1));        
        blake3_hash_output hash = blake3_hasher_finalize(hasher);

        // Print the hash in hexadecimal format

        for (size_t i = 0; i < BLAKE3_OUT_LEN; i++) {
            output[i] = hash.hash[i];
        }

        printf("\nhash value from kernel: \n");
        for (size_t i = 0; i < BLAKE3_OUT_LEN; i++) {
            printf("%02x", hash.hash[i]);
        }
        printf("\n");

        // Convert the hash to a hexadecimal string
        // char hex_string[BLAKE3_OUT_LEN * 2 + 1]; // Each byte becomes two characters, plus a null terminator
        // for (size_t i = 0; i < BLAKE3_OUT_LEN; i++) {
        //     sprintf(&hex_string[i * 2], "%02x", hash.hash[i]);
        // }

        // Compare the hexadecimal string with the target hash
        // const char* target_hash = "e90990bc73ba46c8f2d0b3392c1dc403350d0cfc020738865700b33e74654256";
        // if (strcmp(hex_string, target_hash) == 0) {
        //     printf("Working fine\n");
        // } else {
        //     printf("Ohho, error\n");
        // }

    }
}
        /*Insert Code Here*/

void print_blake3_hasher(blake3_hasher hasher) {
    printf("chunk_state.chunk_counter: %llu\n", hasher.chunk_state.chunk_counter);
    printf("chunk_state.block_len: %u\n", hasher.chunk_state.block_len);
    printf("chunk_state.blocks_compressed: %u\n", hasher.chunk_state.blocks_compressed);
    printf("chunk_state.flags: %u\n", hasher.chunk_state.flags);

    printf("chunk_state.chaining_value: ");
    for(int i = 0; i < 8; i++) {
        printf("%u ", hasher.chunk_state.chaining_value[i]);
    }
    printf("\n");

    printf("chunk_state.block: ");
    for(int i = 0; i < 64; i++) { // Assuming BLAKE3_BLOCK_LEN is 64
        printf("%02x", hasher.chunk_state.block[i]);
        if(i % 16 == 15) printf("\n"); // New line for readability
    }

    printf("key_words: ");
    for(int i = 0; i < 8; i++) {
        printf("%u ", hasher.key_words[i]);
    }
    printf("\n");

    printf("cv_stack: ");
    for(int i = 0; i < 8 * 54; i++) {
        printf("%u ", hasher.cv_stack[i]);
        if(i % 8 == 7) printf("| "); // Separator for readability
    }
    printf("\n");

    printf("cv_stack_len: %u\n", hasher.cv_stack_len);
    printf("flags: %u\n", hasher.flags);
}

static uint32_t IV[8] = {
    0x6A09E667, 0xBB67AE85, 0x3C6EF372, 0xA54FF53A,
    0x510E527F, 0x9B05688C, 0x1F83D9AB, 0x5BE0CD19,
};

static size_t MSG_PERMUTATION[16] = {2, 6,  3,  10, 7, 0,  4,  13,
                                    1, 11, 12, 5,  9, 14, 15, 8};

// used by finalize
inline static uint32_t rotate_right(uint32_t x, int n) {
    //  printf("Function rotate_right called \n");
    return (x >> n) | (x << (32 - n));
}

// used by finalize
// The mixing function, G, which mixes either a column or a diagonal.
inline static void g(uint32_t state[16], size_t a, size_t b, size_t c, size_t d,
                    uint32_t mx, uint32_t my) {
    // printf("Function g called \n");
    state[a] = state[a] + state[b] + mx;
    state[d] = rotate_right(state[d] ^ state[a], 16);
    state[c] = state[c] + state[d];
    state[b] = rotate_right(state[b] ^ state[c], 12);
    state[a] = state[a] + state[b] + my;
    state[d] = rotate_right(state[d] ^ state[a], 8);
    state[c] = state[c] + state[d];
    state[b] = rotate_right(state[b] ^ state[c], 7);
}

// used by finalize
inline static void round_function(uint32_t state[16], uint32_t m[16]) {
    // printf("Function round_function called \n");
    // Mix the columns.
    g(state, 0, 4, 8, 12, m[0], m[1]);
    g(state, 1, 5, 9, 13, m[2], m[3]);
    g(state, 2, 6, 10, 14, m[4], m[5]);
    g(state, 3, 7, 11, 15, m[6], m[7]);
    // Mix the diagonals.
    g(state, 0, 5, 10, 15, m[8], m[9]);
    g(state, 1, 6, 11, 12, m[10], m[11]);
    g(state, 2, 7, 8, 13, m[12], m[13]);
    g(state, 3, 4, 9, 14, m[14], m[15]);
}

// used by finalize
inline static void permute(uint32_t m[16]) {
    // printf("Function permute called \n");
    uint32_t permuted[16];
    for (size_t i = 0; i < 16; i++) {
        permuted[i] = m[MSG_PERMUTATION[i]];
    }
    // Replacing memcpy with a for loop
    for (size_t i = 0; i < 16; i++) {
        m[i] = permuted[i];
    }
}

// used by finalize
inline static void compress(const uint32_t chaining_value[8],
                            const uint32_t block_words[16], uint64_t counter,
                            uint32_t block_len, uint32_t flags,
                            uint32_t out[16]) {
    // printf("Function compress called \n");
    uint32_t state[16] = {
        chaining_value[0],
        chaining_value[1],
        chaining_value[2],
        chaining_value[3],
        chaining_value[4],
        chaining_value[5],
        chaining_value[6],
        chaining_value[7],
        IV[0],
        IV[1],
        IV[2],
        IV[3],
        (uint32_t)counter,
        (uint32_t)(counter >> 32),
        block_len,
        flags,
    };
    uint32_t block[16];
    // Replacing memcpy with a for loop
    for (size_t i = 0; i < 16; i++) {
        block[i] = block_words[i];
    }

    round_function(state, block); // round 1
    permute(block);
    round_function(state, block); // round 2
    permute(block);
    round_function(state, block); // round 3
    permute(block);
    round_function(state, block); // round 4
    permute(block);
    round_function(state, block); // round 5
    permute(block);
    round_function(state, block); // round 6
    permute(block);
    round_function(state, block); // round 7

    for (size_t i = 0; i < 8; i++) {
        state[i] ^= state[i + 8];
        state[i + 8] ^= chaining_value[i];
    }

    
    // Replacing memcpy with a for loop
    for (size_t i = 0; i < 16; i++) {
        out[i] = state[i];
    }

}

// used by finalize
inline static void words_from_little_endian_bytes(const uint8_t block[BLAKE3_BLOCK_LEN],
                                                size_t bytes_len,
                                                uint32_t block_words[16]) {
    assert(bytes_len % 4 == 0);

    for (size_t i = 0, j = 0; i < (bytes_len / 4); i++) {
        block_words[i] = ((uint32_t)(block[j++]));
        block_words[i] += ((uint32_t)(block[j++])) << 8;
        block_words[i] += ((uint32_t)(block[j++])) << 16;
        block_words[i] += ((uint32_t)(block[j++])) << 24;
    }
}


// Each chunk or parent node can produce either an 8-word chaining value or, by
// setting the ROOT flag, any number of final output bytes. The Output struct
// captures the state just prior to choosing between those two possibilities.
typedef struct output {
    uint32_t input_chaining_value[8];
    uint32_t block_words[16];
    uint64_t counter;
    uint32_t block_len;
    uint32_t flags;
} output;


// used by finalize
inline static void output_chaining_value(output self, uint32_t out[8]) {
    // printf("Function output_chaining_value called \n");
    uint32_t out16[16];
    compress(self.input_chaining_value, self.block_words, self.counter,
            self.block_len, self.flags, out16);
    // Replacing memcpy with a for loop
    for (size_t i = 0; i < 8; i++) {
        out[i] = out16[i];
    }
}

// used by finalize
inline static void output_root_bytes(output self, uint8_t out[], size_t out_len) {
    // printf("Function output_root_bytes called \n");
    uint64_t output_block_counter = 0;
    size_t out_index = 0;  // Use an index to replace the pointer out_u8

    while (out_len > 0) {
        uint32_t words[16];
        compress(self.input_chaining_value, self.block_words, output_block_counter,
                self.block_len, self.flags | ROOT, words);

        for (size_t word = 0; word < 16; word++) {
            for (int byte = 0; byte < 4; byte++) {
                if (out_len == 0) {
                    return;
                }
                out[out_index++] = (uint8_t)(words[word] >> (8 * byte));
                out_len--;
            }
        }
        output_block_counter++;
    }
}

// Used by init
inline static _blake3_chunk_state chunk_state_init(const uint32_t key_words[8],
                                                   uint64_t chunk_counter, uint32_t flags) {
    _blake3_chunk_state new_state;

    // Replacing memcpy with a for loop
    for (size_t i = 0; i < 8; i++) {
        new_state.chaining_value[i] = key_words[i];
    }
    new_state.chunk_counter = chunk_counter;
    // Replacing memset with a for loop
    for (size_t i = 0; i < BLAKE3_BLOCK_LEN; i++) { // Assuming BLAKE3_BLOCK_LEN is defined
        new_state.block[i] = 0;
    }
    new_state.block_len = 0;
    new_state.blocks_compressed = 0;
    new_state.flags = flags;

    return new_state;
}



inline static size_t chunk_state_len(_blake3_chunk_state self) {
    return BLAKE3_BLOCK_LEN * (size_t)self.blocks_compressed +
        (size_t)self.block_len;
}

//used by finalize
// inline static uint32_t chunk_state_start_flag(const _blake3_chunk_state *self) {
//     // printf("Function chunk_state_start_flag called \n");
//     if (self->blocks_compressed == 0) {
//         return CHUNK_START;
//     } else {
//         return 0;
//     }
// }

inline static uint32_t chunk_state_start_flag(_blake3_chunk_state self) {
    if (self.blocks_compressed == 0) {
        return CHUNK_START;
    } else {
        return 0;
    }
}


// inline static uint32_t chunk_state_start_flag(const _blake3_chunk_state self) {
//     // printf("Function chunk_state_start_flag called \n");
//     if (self.blocks_compressed == 0) {
//         return CHUNK_START;
//     } else {
//         return 0;
//     }
// }

//used by update
// inline static void chunk_state_update(_blake3_chunk_state *self,
//                                     const void *input, size_t input_len) {
//     // printf("Function chunk_state_update called \n");                                      
//     const uint8_t *input_u8 = (const uint8_t *)input;
//     while (input_len > 0) {
//         // If the block buffer is full, compress it and clear it. More input is
//         // coming, so this compression is not CHUNK_END.
//         if (self->block_len == BLAKE3_BLOCK_LEN) {
//             uint32_t block_words[16];
//             words_from_little_endian_bytes(self->block, BLAKE3_BLOCK_LEN,
//                                             block_words);
//             uint32_t out16[16];
//             compress(self->chaining_value, block_words, self->chunk_counter,
//                     BLAKE3_BLOCK_LEN, self->flags | chunk_state_start_flag(self),
//                     out16);
//             // Replacing memcpy with a for loop
//             for (size_t i = 0; i < 8; i++) {
//                 self->chaining_value[i] = out16[i];
//             }
//             self->blocks_compressed++;
//             // Replacing memset with a for loop
//             for (size_t i = 0; i < 64; i++) {
//                 self->block[i] = 0;
//             }
//             self->block_len = 0;
//         }

//         // Copy input bytes into the block buffer.
//         size_t want = BLAKE3_BLOCK_LEN - (size_t)self->block_len;
//         size_t take = want;
//         if (input_len < want) {
//             take = input_len;
//         }
//         // Replacing memcpy with a for loop
//         for (size_t i = 0; i < take; i++) {
//             self->block[(size_t)self->block_len + i] = input_u8[i];
//         }
//         self->block_len += (uint8_t)take;
//         input_u8 += take;
//         input_len -= take;
//     }
// }


inline static _blake3_chunk_state chunk_state_update(_blake3_chunk_state self,
                                                     const void *input, size_t input_len) {
    const uint8_t *input_u8 = (const uint8_t *)input;
    while (input_len > 0) {
        if (self.block_len == BLAKE3_BLOCK_LEN) {
            uint32_t block_words[16];
            words_from_little_endian_bytes(self.block, BLAKE3_BLOCK_LEN, block_words);
            uint32_t out16[16];
            compress(self.chaining_value, block_words, self.chunk_counter,
                     BLAKE3_BLOCK_LEN, self.flags | chunk_state_start_flag(self),
                     out16);
            for (size_t i = 0; i < 8; i++) {
                self.chaining_value[i] = out16[i];
            }
            self.blocks_compressed++;
            for (size_t i = 0; i < BLAKE3_BLOCK_LEN; i++) {
                self.block[i] = 0;
            }
            self.block_len = 0;
        }

        size_t want = BLAKE3_BLOCK_LEN - self.block_len;
        size_t take = (input_len < want) ? input_len : want;
        for (size_t i = 0; i < take; i++) {
            self.block[self.block_len + i] = input_u8[i];
        }
        //check here
        // printf("\nchecking block_len udpdate:" );
        // printf("Input: %zu\n", input_len);
        // printf("Want: %zu\n", want);
        // printf("Take: %zu\n", take);
        self.block_len += (uint8_t)take;
        // printf("block_len: %u\n", self.block_len);
        input_u8 += take;
        input_len -= take;
    }
    return self;
}



// inline static _blake3_chunk_state chunk_state_update(_blake3_chunk_state self,
//                                                      const void *input, size_t input_len) {
//     // printf("Function chunk_state_update called \n");
//     const uint8_t *input_u8 = (const uint8_t *)input;

//     while (input_len > 0) {
//         if (self.block_len == BLAKE3_BLOCK_LEN) {
//             uint32_t block_words[16];
//             words_from_little_endian_bytes(self.block, BLAKE3_BLOCK_LEN, block_words);
//             uint32_t out16[16];
//             compress(self.chaining_value, block_words, self.chunk_counter,
//                      BLAKE3_BLOCK_LEN, self.flags | chunk_state_start_flag(self),
//                      out16);
//             for (size_t i = 0; i < 8; i++) {
//                 self.chaining_value[i] = out16[i];
//             }
//             self.blocks_compressed++;
//             for (size_t i = 0; i < BLAKE3_BLOCK_LEN; i++) {
//                 self.block[i] = 0;
//             }
//             self.block_len = 0;
//         }

//         size_t want = BLAKE3_BLOCK_LEN - (size_t)self.block_len;
//         size_t take = (input_len < want) ? input_len : want;
//         for (size_t i = 0; i < take; i++) {
//             self.block[self.block_len + i] = input_u8[i];
//         }
//         self.block_len += (uint8_t)take;
//         input_u8 += take;
//         input_len -= take;
//     }

//     return self;
// }

// used by finalize
// inline static output chunk_state_output(const _blake3_chunk_state *self) {
//     // printf("Function chunk_state_output called \n");    
//     output ret;
//     // Replacing memcpy with a for loop
//     for (size_t i = 0; i < 8; i++) {
//         ret.input_chaining_value[i] = self->chaining_value[i];
//     }
//     words_from_little_endian_bytes(self->block, sizeof(self->block),
//                                     ret.block_words);
//     ret.counter = self->chunk_counter;
//     ret.block_len = (uint32_t)self->block_len;
//     ret.flags = self->flags | chunk_state_start_flag(self) | CHUNK_END;
//     return ret;
// }

inline static output chunk_state_output(_blake3_chunk_state self) {
    // printf("Function chunk_state_output called \n");    
    output ret;

    // Replacing memcpy with a for loop
    for (size_t i = 0; i < 8; i++) {
        ret.input_chaining_value[i] = self.chaining_value[i];
    }

    words_from_little_endian_bytes(self.block, sizeof(self.block), ret.block_words);
    ret.counter = self.chunk_counter;
    ret.block_len = (uint32_t)self.block_len;
    ret.flags = self.flags | chunk_state_start_flag(self) | CHUNK_END;

    return ret;
}




// used by finalize
inline static output parent_output(const uint32_t left_child_cv[8],
                                const uint32_t right_child_cv[8],
                                const uint32_t key_words[8],
                                uint32_t flags) {
    // printf("Function parent_output called \n");    
    output ret;
    // Replacing memcpy with a for loop
    for (size_t i = 0; i < 8; i++) {
        ret.input_chaining_value[i] = key_words[i];
    }

    // Replacing memcpy with a for loop
    for (size_t i = 0; i < 8; i++) {
        ret.block_words[i] = left_child_cv[i];
    }

    // Replacing memcpy with a for loop
    for (size_t i = 0; i < 8; i++) {
        ret.block_words[8 + i] = right_child_cv[i];
    }
    ret.counter = 0; // Always 0 for parent nodes.
    ret.block_len =
        BLAKE3_BLOCK_LEN; // Always BLAKE3_BLOCK_LEN (64) for parent nodes.
    ret.flags = PARENT | flags;
    return ret;
}

//used by update
inline static void parent_cv(const uint32_t left_child_cv[8],
                            const uint32_t right_child_cv[8],
                            const uint32_t key_words[8], uint32_t flags,
                            uint32_t out[8]) {
    // printf("Function parent_cv called \n"); 
    output o = parent_output(left_child_cv, right_child_cv, key_words, flags);
    // We only write to `out` after we've read the inputs. That makes it safe for
    // `out` to alias an input, which we do below.
    output_chaining_value(o, out);
}

// Used by init
// inline static void hasher_init_internal(blake3_hasher *self,
//                                         const uint32_t key_words[8],
//                                         uint32_t flags) {
//     // printf("Function hasher_init_internal called \n");
//     self->chunk_state = chunk_state_init(key_words, 0, flags);
    
//     // Replacing memcpy with a for loop
//     for (size_t i = 0; i < 8; i++) {
//         self->key_words[i] = key_words[i];
//     }
//     self->cv_stack_len = 0;
//     self->flags = flags;
// }

// working fine
inline static blake3_hasher hasher_init_internal(const uint32_t key_words[8], uint32_t flags) {
    // printf("Function hasher_init_internal called \n");
    blake3_hasher self;
    self.chunk_state = chunk_state_init(key_words, 0, flags);
    
    for (size_t i = 0; i < 8; i++) {
        self.key_words[i] = key_words[i];
    }
    self.cv_stack_len = 0;
    self.flags = flags;

    return self;
}


// Used by init
// Construct a new `Hasher` for the regular hash function.
// void blake3_hasher_init(blake3_hasher *self) {
//     // printf("Function blake3_hasher_init called \n");
//     hasher_init_internal(self, IV, 0);
// }

// working fine
blake3_hasher blake3_hasher_init() {
    // printf("Function blake3_hasher_init called \n");
    blake3_hasher self = hasher_init_internal(IV, 0);
    return self;
}



// inline static void hasher_push_stack(blake3_hasher *self,
//                                     const uint32_t cv[8]) {
//     // printf("Function hasher_push_stack called \n");                                    
//     // Replacing memcpy with a for loop

//     size_t start_idx = (size_t)self->cv_stack_len * 8;
//     for (size_t i = 0; i < 8; i++) {
//         self->cv_stack[start_idx + i] = cv[i];
//     }
//     self->cv_stack_len++;
// }

blake3_hasher hasher_push_stack(blake3_hasher self, const uint32_t cv[8]) {
    size_t start_idx = self.cv_stack_len * 8;
    for (size_t i = 0; i < 8; i++) {
        self.cv_stack[start_idx + i] = cv[i];
    }
    self.cv_stack_len++;
    return self;
}

// Returns a pointer to the popped CV, which is valid until the next push.
// inline static const uint32_t *hasher_pop_stack(blake3_hasher *self) {
//     // printf("Function hasher_pop_stack called \n"); 
//     self->cv_stack_len--;
//     return &self->cv_stack[(size_t)self->cv_stack_len * 8];
// }

hasher_pop_result hasher_pop_stack(blake3_hasher self) {
    self.cv_stack_len--;
    hasher_pop_result result;
    result.updated_hasher = self;

    // Use a for loop to copy the elements
    size_t start_idx = self.cv_stack_len * 8;
    for (size_t i = 0; i < 8; i++) {
        result.popped_cv[i] = self.cv_stack[start_idx + i];
    }

    return result;
}


// Ask prithviraj how to handle this?
// used by update
// Section 5.1.2 of the BLAKE3 spec explains this algorithm in more detail.
// inline static void hasher_add_chunk_cv(blake3_hasher *self, uint32_t new_cv[8],
//                                     uint64_t total_chunks) {
//     // This chunk might complete some subtrees. For each completed subtree, its
//     // left child will be the current top entry in the CV stack, and its right
//     // child will be the current value of `new_cv`. Pop each left child off the
//     // stack, merge it with `new_cv`, and overwrite `new_cv` with the result.
//     // After all these merges, push the final value of `new_cv` onto the stack.
//     // The number of completed subtrees is given by the number of trailing 0-bits
//     // in the new total number of chunks.
//     // printf("Function hasher_add_chunk_cv called \n"); 
//     while ((total_chunks & 1) == 0) {
//         parent_cv(hasher_pop_stack(self), new_cv, self->key_words, self->flags,
//                 new_cv);
//         total_chunks >>= 1;
//     }
//     hasher_push_stack(self, new_cv);
// }

blake3_hasher hasher_add_chunk_cv(blake3_hasher self, uint32_t new_cv[8], uint64_t total_chunks) {
    while ((total_chunks & 1) == 0) {
        hasher_pop_result pop_result = hasher_pop_stack(self);
        self = pop_result.updated_hasher;
        parent_cv(pop_result.popped_cv, new_cv, self.key_words, self.flags, new_cv);
        total_chunks >>= 1;
    }
    self = hasher_push_stack(self, new_cv);
    return self;
}


// Used by update
// Add input to the hash state. This can be called any number of times.
// void blake3_hasher_update(blake3_hasher *self, const void *input,
//                         size_t input_len) {
//     // printf("Function blake3_hasher_update called \n"); 
//     const uint8_t *input_u8 = (const uint8_t *)input;
//     while (input_len > 0) {
//         // If the current chunk is complete, finalize it and reset the chunk state.
//         // More input is coming, so this chunk is not ROOT.
//         if (chunk_state_len(self->chunk_state) == BLAKE3_CHUNK_LEN) {
//             output chunk_output = chunk_state_output(&self->chunk_state);
//             uint32_t chunk_cv[8];
//             output_chaining_value(chunk_output, chunk_cv);
//             uint64_t total_chunks = self->chunk_state.chunk_counter + 1;
//             hasher_add_chunk_cv(self, chunk_cv, total_chunks);
//             self->chunk_state = chunk_state_init(self->key_words, total_chunks,
//                             self->flags);
//         }

//         // Compress input bytes into the current chunk state.
//         size_t want = BLAKE3_CHUNK_LEN - chunk_state_len(self->chunk_state);
//         size_t take = want;
//         if (input_len < want) {
//             take = input_len;
//         }
//         chunk_state_update(&self->chunk_state, input_u8, take);
//         input_u8 += take;
//         input_len -= take;
//     }
// }

// Ask prithviraj how to handle this?
blake3_hasher blake3_hasher_update(blake3_hasher self, const void *input, size_t input_len) {
    const uint8_t *input_u8 = (const uint8_t *)input;
    while (input_len > 0) {
        if (chunk_state_len(self.chunk_state) == BLAKE3_CHUNK_LEN) {
            output chunk_output = chunk_state_output(self.chunk_state);
            uint32_t chunk_cv[8];
            output_chaining_value(chunk_output, chunk_cv);
            uint64_t total_chunks = self.chunk_state.chunk_counter + 1;
            self = hasher_add_chunk_cv(self, chunk_cv, total_chunks);
            self.chunk_state = chunk_state_init(self.key_words, total_chunks, self.flags);
        }

        size_t want = BLAKE3_CHUNK_LEN - chunk_state_len(self.chunk_state);
        size_t take = (input_len < want) ? input_len : want;
        self.chunk_state = chunk_state_update(self.chunk_state, input_u8, take);
        input_u8 += take;
        input_len -= take;
    }
    return self;
}

// need to work on it 
// used by finalize
// Finalize the hash and write any number of output bytes.
// void blake3_hasher_finalize(const blake3_hasher *self, void *out,
//                             size_t out_len) {
//     // Starting with the output from the current chunk, compute all the parent
//     // chaining values along the right edge of the tree, until we have the root
//     // output.
//     // printf("Function blake3_hasher_finalize called \n"); 
//     output current_output = chunk_state_output(self->chunk_state);
//     size_t parent_nodes_remaining = (size_t)self->cv_stack_len;
//     while (parent_nodes_remaining > 0) {
//         parent_nodes_remaining--;
//         uint32_t current_cv[8];
//         output_chaining_value(current_output, current_cv);
//         current_output = parent_output(&self->cv_stack[parent_nodes_remaining * 8],
//                                     current_cv, self->key_words, self->flags);
//     }
//     output_root_bytes(current_output, out, out_len);
// }


blake3_hash_output blake3_hasher_finalize(blake3_hasher self) {
    blake3_hash_output hash_output;
    output current_output = chunk_state_output(self.chunk_state);
    size_t parent_nodes_remaining = (size_t)self.cv_stack_len;
    while (parent_nodes_remaining > 0) {
        parent_nodes_remaining--;
        uint32_t current_cv[8];
        output_chaining_value(current_output, current_cv);
        current_output = parent_output(self.cv_stack + (parent_nodes_remaining * 8),
                                    current_cv, self.key_words, self.flags);
    }
    output_root_bytes(current_output, hash_output.hash, BLAKE3_BLOCK_LEN);


    return hash_output;
}

// void blake3_hasher_finalize(blake3_hasher self, uint8_t out[], size_t out_len) {
//     // printf("Function blake3_hasher_finalize called \n");
//     output current_output = chunk_state_output(self.chunk_state);
//     size_t parent_nodes_remaining = (size_t)self.cv_stack_len;
//     while (parent_nodes_remaining > 0) {
//         parent_nodes_remaining--;
//         uint32_t current_cv[8];
//         output_chaining_value(current_output, current_cv);
//         current_output = parent_output(&self.cv_stack[parent_nodes_remaining * 8],
//                                     current_cv, self.key_words, self.flags);
//     }
//     output_root_bytes(current_output, out, out_len);
// }   





        // static uint32_t IV[8] = {
        //     0x6A09E667, 0xBB67AE85, 0x3C6EF372, 0xA54FF53A,
        //     0x510E527F, 0x9B05688C, 0x1F83D9AB, 0x5BE0CD19,
        // };

        // static size_t MSG_PERMUTATION[16] = {2, 6,  3,  10, 7, 0,  4,  13,
        //                                     1, 11, 12, 5,  9, 14, 15, 8};

        // // used by finalize
        // inline static uint32_t rotate_right(uint32_t x, int n) {
        //     //  printf("Function rotate_right called \n");
        //     return (x >> n) | (x << (32 - n));
        // }

        // // used by finalize
        // // The mixing function, G, which mixes either a column or a diagonal.
        // inline static void g(uint32_t state[16], size_t a, size_t b, size_t c, size_t d,
        //                     uint32_t mx, uint32_t my) {
        //     // printf("Function g called \n");
        //     state[a] = state[a] + state[b] + mx;
        //     state[d] = rotate_right(state[d] ^ state[a], 16);
        //     state[c] = state[c] + state[d];
        //     state[b] = rotate_right(state[b] ^ state[c], 12);
        //     state[a] = state[a] + state[b] + my;
        //     state[d] = rotate_right(state[d] ^ state[a], 8);
        //     state[c] = state[c] + state[d];
        //     state[b] = rotate_right(state[b] ^ state[c], 7);
        // }

        // // used by finalize
        // inline static void round_function(uint32_t state[16], uint32_t m[16]) {
        //     // printf("Function round_function called \n");
        //     // Mix the columns.
        //     g(state, 0, 4, 8, 12, m[0], m[1]);
        //     g(state, 1, 5, 9, 13, m[2], m[3]);
        //     g(state, 2, 6, 10, 14, m[4], m[5]);
        //     g(state, 3, 7, 11, 15, m[6], m[7]);
        //     // Mix the diagonals.
        //     g(state, 0, 5, 10, 15, m[8], m[9]);
        //     g(state, 1, 6, 11, 12, m[10], m[11]);
        //     g(state, 2, 7, 8, 13, m[12], m[13]);
        //     g(state, 3, 4, 9, 14, m[14], m[15]);
        // }

        // // used by finalize
        // inline static void permute(uint32_t m[16]) {
        //     // printf("Function permute called \n");
        //     uint32_t permuted[16];
        //     for (size_t i = 0; i < 16; i++) {
        //         permuted[i] = m[MSG_PERMUTATION[i]];
        //     }
        //     // Replacing memcpy with a for loop
        //     for (size_t i = 0; i < 16; i++) {
        //         m[i] = permuted[i];
        //     }
        // }

        // // used by finalize
        // inline static void compress(const uint32_t chaining_value[8],
        //                             const uint32_t block_words[16], uint64_t counter,
        //                             uint32_t block_len, uint32_t flags,
        //                             uint32_t out[16]) {
        //     // printf("Function compress called \n");
        //     uint32_t state[16] = {
        //         chaining_value[0],
        //         chaining_value[1],
        //         chaining_value[2],
        //         chaining_value[3],
        //         chaining_value[4],
        //         chaining_value[5],
        //         chaining_value[6],
        //         chaining_value[7],
        //         IV[0],
        //         IV[1],
        //         IV[2],
        //         IV[3],
        //         (uint32_t)counter,
        //         (uint32_t)(counter >> 32),
        //         block_len,
        //         flags,
        //     };
        //     uint32_t block[16];
        //     // Replacing memcpy with a for loop
        //     for (size_t i = 0; i < 16; i++) {
        //         block[i] = block_words[i];
        //     }

        //     round_function(state, block); // round 1
        //     permute(block);
        //     round_function(state, block); // round 2
        //     permute(block);
        //     round_function(state, block); // round 3
        //     permute(block);
        //     round_function(state, block); // round 4
        //     permute(block);
        //     round_function(state, block); // round 5
        //     permute(block);
        //     round_function(state, block); // round 6
        //     permute(block);
        //     round_function(state, block); // round 7

        //     for (size_t i = 0; i < 8; i++) {
        //         state[i] ^= state[i + 8];
        //         state[i + 8] ^= chaining_value[i];
        //     }

            
        //     // Replacing memcpy with a for loop
        //     for (size_t i = 0; i < 16; i++) {
        //         out[i] = state[i];
        //     }

        // }

        
        // // inline static void words_from_little_endian_bytes(const void *bytes,
        // //                                                 size_t bytes_len,
        // //                                                 uint32_t *out) {
        // //     // printf("Function words_from_little_endian_bytes called \n");
        // //     assert(bytes_len % 4 == 0);
        // //     const uint8_t *u8_ptr = (const uint8_t *)bytes;
        // //     for (size_t i = 0; i < (bytes_len / 4); i++) {
        // //         out[i] = ((uint32_t)(*u8_ptr++));
        // //         out[i] += ((uint32_t)(*u8_ptr++)) << 8;
        // //         out[i] += ((uint32_t)(*u8_ptr++)) << 16;
        // //         out[i] += ((uint32_t)(*u8_ptr++)) << 24;
        // //     }
        // // }

        // // used by finalize
        // // changed pointer and tested - working fine
        // inline static void words_from_little_endian_bytes(uint8_t block[BLAKE3_BLOCK_LEN],
        //                                                 size_t bytes_len,
        //                                                 uint32_t block_words[16]) {
        //     assert(bytes_len % 4 == 0);

        //     for (size_t i = 0, j = 0; i < (bytes_len / 4); i++) {
        //         block_words[i] = ((uint32_t)(block[j++]));
        //         block_words[i] += ((uint32_t)(block[j++])) << 8;
        //         block_words[i] += ((uint32_t)(block[j++])) << 16;
        //         block_words[i] += ((uint32_t)(block[j++])) << 24;
        //     }
        // }


        // // Each chunk or parent node can produce either an 8-word chaining value or, by
        // // setting the ROOT flag, any number of final output bytes. The Output struct
        // // captures the state just prior to choosing between those two possibilities.
        // typedef struct output {
        //     uint32_t input_chaining_value[8];
        //     uint32_t block_words[16];
        //     uint64_t counter;
        //     uint32_t block_len;
        //     uint32_t flags;
        // } output;


        // // used by finalize
        // // inline static void output_chaining_value(const output *self, uint32_t out[8]) {
        // //     // printf("Function output_chaining_value called \n");
        // //     uint32_t out16[16];
        // //     compress(self->input_chaining_value, self->block_words, self->counter,
        // //             self->block_len, self->flags, out16);
        // //     // Replacing memcpy with a for loop
        // //     for (size_t i = 0; i < 8; i++) {
        // //         out[i] = out16[i];
        // //     }
        // // }

        // // used by finalize
        // // changed pointer and tested - working fine
        // inline static void output_chaining_value(output self, uint32_t out[8]) {
        //     // printf("Function output_chaining_value called \n");
        //     uint32_t out16[16];
        //     compress(self.input_chaining_value, self.block_words, self.counter,
        //             self.block_len, self.flags, out16);
        //     // Replacing memcpy with a for loop
        //     for (size_t i = 0; i < 8; i++) {
        //         out[i] = out16[i];
        //     }
        // }

        // // changed pointer and tested - working fine
        // inline static void output_root_bytes(output self, uint8_t out[], size_t out_len) {
        //     // printf("Function output_root_bytes called \n");
        //     uint64_t output_block_counter = 0;
        //     size_t out_index = 0;  // Use an index to replace the pointer out_u8

        //     while (out_len > 0) {
        //         uint32_t words[16];
        //         compress(self.input_chaining_value, self.block_words, output_block_counter,
        //                 self.block_len, self.flags | ROOT, words);

        //         for (size_t word = 0; word < 16; word++) {
        //             for (int byte = 0; byte < 4; byte++) {
        //                 if (out_len == 0) {
        //                     return;
        //                 }
        //                 out[out_index++] = (uint8_t)(words[word] >> (8 * byte));
        //                 out_len--;
        //             }
        //         }
        //         output_block_counter++;
        //     }
        // }

        // // Used by init
        // // changed pointer and tested - working fine
        // inline static _blake3_chunk_state chunk_state_init(const uint32_t key_words[8],
        //                                                 uint64_t chunk_counter, uint32_t flags) {
        //     _blake3_chunk_state new_state;

        //     // Replacing memcpy with a for loop
        //     for (size_t i = 0; i < 8; i++) {
        //         new_state.chaining_value[i] = key_words[i];
        //     }
        //     new_state.chunk_counter = chunk_counter;
        //     // Replacing memset with a for loop
        //     for (size_t i = 0; i < BLAKE3_BLOCK_LEN; i++) { // Assuming BLAKE3_BLOCK_LEN is defined
        //         new_state.block[i] = 0;
        //     }
        //     new_state.block_len = 0;
        //     new_state.blocks_compressed = 0;
        //     new_state.flags = flags;

        //     return new_state;
        // }

        // //used by update
        // // changed pointer and tested - working fine
        // inline static size_t chunk_state_len(_blake3_chunk_state self) {
        //     return BLAKE3_BLOCK_LEN * (size_t)self.blocks_compressed +
        //         (size_t)self.block_len;
        // }

        // //used by finalize
        // // changed pointer and tested - working fine
        // inline static uint32_t chunk_state_start_flag(_blake3_chunk_state self) {
        //     if (self.blocks_compressed == 0) {
        //         return CHUNK_START;
        //     } else {
        //         return 0;
        //     }
        // }

        // //used by update
        // // changed pointer and tested - working fine
        // inline static _blake3_chunk_state chunk_state_update(_blake3_chunk_state self,
        //                                                     const void *input, size_t input_len) {
        //     const uint8_t *input_u8 = (const uint8_t *)input;
        //     while (input_len > 0) {
        //         if (self.block_len == BLAKE3_BLOCK_LEN) {
        //             uint32_t block_words[16];
        //             words_from_little_endian_bytes(self.block, BLAKE3_BLOCK_LEN, block_words);
        //             uint32_t out16[16];
        //             compress(self.chaining_value, block_words, self.chunk_counter,
        //                     BLAKE3_BLOCK_LEN, self.flags | chunk_state_start_flag(self),
        //                     out16);
        //             for (size_t i = 0; i < 8; i++) {
        //                 self.chaining_value[i] = out16[i];
        //             }
        //             self.blocks_compressed++;
        //             for (size_t i = 0; i < BLAKE3_BLOCK_LEN; i++) {
        //                 self.block[i] = 0;
        //             }
        //             self.block_len = 0;
        //         }

        //         size_t want = BLAKE3_BLOCK_LEN - self.block_len;
        //         size_t take = (input_len < want) ? input_len : want;
        //         for (size_t i = 0; i < take; i++) {
        //             self.block[self.block_len + i] = input_u8[i];
        //         }
        //         self.block_len += (uint8_t)take;
        //         input_u8 += take;
        //         input_len -= take;
        //     }
        //     return self;
        // }

        // // used by finalize
        // inline static output chunk_state_output(_blake3_chunk_state self) {
        //     // printf("Function chunk_state_output called \n");    
        //     output ret;

        //     // Replacing memcpy with a for loop
        //     for (size_t i = 0; i < 8; i++) {
        //         ret.input_chaining_value[i] = self.chaining_value[i];
        //     }
        //     words_from_little_endian_bytes(self.block, sizeof(self.block), ret.block_words);
        //     ret.counter = self.chunk_counter;
        //     ret.block_len = (uint32_t)self.block_len;
        //     ret.flags = self.flags | chunk_state_start_flag(self) | CHUNK_END;

        //     return ret;
        // }

        // // used by finalize
        // inline static output parent_output(const uint32_t left_child_cv[8],
        //                                 const uint32_t right_child_cv[8],
        //                                 const uint32_t key_words[8],
        //                                 uint32_t flags) {
        //     // printf("Function parent_output called \n");    
        //     output ret;
        //     // Replacing memcpy with a for loop
        //     for (size_t i = 0; i < 8; i++) {
        //         ret.input_chaining_value[i] = key_words[i];
        //     }

        //     // Replacing memcpy with a for loop
        //     for (size_t i = 0; i < 8; i++) {
        //         ret.block_words[i] = left_child_cv[i];
        //     }

        //     // Replacing memcpy with a for loop
        //     for (size_t i = 0; i < 8; i++) {
        //         ret.block_words[8 + i] = right_child_cv[i];
        //     }
        //     ret.counter = 0; // Always 0 for parent nodes.
        //     ret.block_len =
        //         BLAKE3_BLOCK_LEN; // Always BLAKE3_BLOCK_LEN (64) for parent nodes.
        //     ret.flags = PARENT | flags;
        //     return ret;
        // }

        // //used by update
        // inline static void parent_cv(const uint32_t left_child_cv[8],
        //                             const uint32_t right_child_cv[8],
        //                             const uint32_t key_words[8], uint32_t flags,
        //                             uint32_t out[8]) {
        //     // printf("Function parent_cv called \n"); 
        //     output o = parent_output(left_child_cv, right_child_cv, key_words, flags);
        //     // We only write to `out` after we've read the inputs. That makes it safe for
        //     // `out` to alias an input, which we do below.
        //     output_chaining_value(o, out);
        // }

        // // Used by init
        // // inline static void hasher_init_internal(blake3_hasher *self,
        // //                                         const uint32_t key_words[8],
        // //                                         uint32_t flags) {
        // //     // printf("Function hasher_init_internal called \n");
        // //     self->chunk_state = chunk_state_init(key_words, 0, flags);
        // //     // Replacing memcpy with a for loop
        // //     for (size_t i = 0; i < 8; i++) {
        // //         self->key_words[i] = key_words[i];
        // //     }
        // //     self->cv_stack_len = 0;
        // //     self->flags = flags;
        // // }

        // // Used by init
        // // changed pointer and tested - working fine
        // inline static blake3_hasher hasher_init_internal(const uint32_t key_words[8], uint32_t flags) {
        //     // printf("Function hasher_init_internal called \n");
        //     blake3_hasher self;
        //     self.chunk_state = chunk_state_init(key_words, 0, flags);
            
        //     for (size_t i = 0; i < 8; i++) {
        //         self.key_words[i] = key_words[i];
        //     }
        //     self.cv_stack_len = 0;
        //     self.flags = flags;

        //     return self;
        // }

        

        // // Used by init
        // // Construct a new `Hasher` for the regular hash function.
        // // void blake3_hasher_init(blake3_hasher *self) {
        // //     // printf("Function blake3_hasher_init called \n");
        // //     hasher_init_internal(self, IV, 0);
        // // }


        // // Used by init
        // // changed pointer and tested - working fine
        // blake3_hasher blake3_hasher_init() {
        //     // printf("Function blake3_hasher_init called \n");
        //     blake3_hasher self = hasher_init_internal(IV, 0);
        //     return self;
        // }

        // // changed pointer and tested - working fine
        // blake3_hasher hasher_push_stack(blake3_hasher self, const uint32_t cv[8]) {
        //     size_t start_idx = self.cv_stack_len * 8;
        //     for (size_t i = 0; i < 8; i++) {
        //         self.cv_stack[start_idx + i] = cv[i];
        //     }
        //     self.cv_stack_len++;
        //     return self;
        // }

        // // Returns a pointer to the popped CV, which is valid until the next push.
        // // changed pointer and tested - working fine
        // hasher_pop_result hasher_pop_stack(blake3_hasher self) {
        //     self.cv_stack_len--;
        //     hasher_pop_result result;
        //     result.updated_hasher = self;

        //     // Use a for loop to copy the elements
        //     size_t start_idx = self.cv_stack_len * 8;
        //     for (size_t i = 0; i < 8; i++) {
        //         result.popped_cv[i] = self.cv_stack[start_idx + i];
        //     }

        //     return result;
        // }

        // // used by update
        // // Section 5.1.2 of the BLAKE3 spec explains this algorithm in more detail.
        // // changed pointer and tested - working fine
        // blake3_hasher hasher_add_chunk_cv(blake3_hasher self, uint32_t new_cv[8], uint64_t total_chunks) {

        //     // This chunk might complete some subtrees. For each completed subtree, its
        //     // left child will be the current top entry in the CV stack, and its right
        //     // child will be the current value of `new_cv`. Pop each left child off the
        //     // stack, merge it with `new_cv`, and overwrite `new_cv` with the result.
        //     // After all these merges, push the final value of `new_cv` onto the stack.
        //     // The number of completed subtrees is given by the number of trailing 0-bits
        //     // in the new total number of chunks.
        //     // printf("Function hasher_add_chunk_cv called \n");
            
        //     while ((total_chunks & 1) == 0) {
        //         hasher_pop_result pop_result = hasher_pop_stack(self);
        //         self = pop_result.updated_hasher;
        //         parent_cv(pop_result.popped_cv, new_cv, self.key_words, self.flags, new_cv);
        //         total_chunks >>= 1;
        //     }
        //     self = hasher_push_stack(self, new_cv);
        //     return self;
        // }

        // // Used by update
        // // Add input to the hash state. This can be called any number of times.
        // // changed pointer and tested - working fine
        // // need to work - make input pointer free - i dont think input should be made pointer free coz we get input as a pointer
        // blake3_hasher blake3_hasher_update(blake3_hasher self, const void *input, size_t input_len) {
        //     const uint8_t *input_u8 = (const uint8_t *)input;
        //     while (input_len > 0) {
        //         if (chunk_state_len(self.chunk_state) == BLAKE3_CHUNK_LEN) {
        //             output chunk_output = chunk_state_output(self.chunk_state);
        //             uint32_t chunk_cv[8];
        //             output_chaining_value(chunk_output, chunk_cv);
        //             uint64_t total_chunks = self.chunk_state.chunk_counter + 1;
        //             self = hasher_add_chunk_cv(self, chunk_cv, total_chunks);
        //             self.chunk_state = chunk_state_init(self.key_words, total_chunks, self.flags);
        //         }

        //         size_t want = BLAKE3_CHUNK_LEN - chunk_state_len(self.chunk_state);
        //         size_t take = (input_len < want) ? input_len : want;
        //         self.chunk_state = chunk_state_update(self.chunk_state, input_u8, take);
        //         input_u8 += take;
        //         input_len -= take;
        //     }
        //     return self;
        // }

        // // used by finalize
        // // Finalize the hash and write any number of output bytes.
        // // changed pointer and tested - working fine
        // blake3_hash_output blake3_hasher_finalize(blake3_hasher self) {
        //     blake3_hash_output hash_output;
        //     output current_output = chunk_state_output(self.chunk_state);
        //     size_t parent_nodes_remaining = (size_t)self.cv_stack_len;
        //     while (parent_nodes_remaining > 0) {
        //         parent_nodes_remaining--;
        //         uint32_t current_cv[8];
        //         output_chaining_value(current_output, current_cv);
        //         current_output = parent_output(self.cv_stack + (parent_nodes_remaining * 8),
        //                                     current_cv, self.key_words, self.flags);
        //     }
        //     output_root_bytes(current_output, hash_output.hash, BLAKE3_BLOCK_LEN);


        //     return hash_output;
        // }
