#ifndef _BLAKE3_REFERENCE_IMPL_H
#define _BLAKE3_REFERENCE_IMPL_H

#include <stddef.h>
#include <stdint.h>

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

// void blake3_hasher_init(blake3_hasher *self);

blake3_hasher blake3_hasher_init();

// void blake3_hasher_update(blake3_hasher *self, const void *input,
//                           size_t input_len);

blake3_hasher blake3_hasher_update(blake3_hasher self, const void *input,
                          size_t input_len);

// void blake3_hasher_finalize(const blake3_hasher *self, void *out,
//                             size_t out_len);

blake3_hash_output blake3_hasher_finalize(const blake3_hasher self);

void print_blake3_hasher(blake3_hasher hasher);

#endif // _BLAKE3_REFERENCE_IMPL_H
