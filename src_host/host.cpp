#include "host.h"
#include "constants.h"
#include <vector> 
#include <random>
#include <assert.h>
#include <string.h>
#include <iostream>
#include <iomanip>

// blake3 library constants, structs and function declarations

#define CHUNK_START 1 << 0
#define CHUNK_END 1 << 1
#define PARENT 1 << 2
#define ROOT 1 << 3

#define BLAKE3_OUT_LEN 32
#define BLAKE3_BLOCK_LEN 64
#define BLAKE3_CHUNK_LEN 1024

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

int main(int argc, char** argv) {

    if (argc != 2) {
        std::cout << "Usage: " << argv[0] << " <XCLBIN File>" << std::endl;
        return EXIT_FAILURE;
    }
    
    clock_t htod, dtoh, comp; 

    /*====================================================CL===============================================================*/

    std::string binaryFile = argv[1];
    cl_int err;
    cl::Context context;
    cl::Kernel krnl1, krnl2;
    cl::CommandQueue q;
    
    auto devices = get_xil_devices();
    auto fileBuf = read_binary_file(binaryFile);
    cl::Program::Binaries bins{{fileBuf.data(), fileBuf.size()}};
    bool valid_device = false;
    for (unsigned int i = 0; i < devices.size(); i++) {
        auto device = devices[i];
        OCL_CHECK(err, context = cl::Context(device, nullptr, nullptr, nullptr, &err));
        OCL_CHECK(err, q = cl::CommandQueue(context, device, 0, &err));
        std::cout << "Trying to program device[" << i << "]: " << device.getInfo<CL_DEVICE_NAME>() << std::endl;
        cl::Program program(context, {device}, bins, nullptr, &err);
        if (err != CL_SUCCESS) {
            std::cout << "Failed to program device[" << i << "] with xclbin file!\n";
        } else {
            std::cout << "Device[" << i << "]: program successful!\n";
            std::cout << "Setting CU(s) up..." << std::endl; 
            OCL_CHECK(err, krnl1 = cl::Kernel(program, "krnl", &err));
            valid_device = true;
            break; // we break because we found a valid device
        }
    }
    if (!valid_device) {
        std::cout << "Failed to program any device found, exit!\n";
        exit(EXIT_FAILURE);
    }

    /*====================================================INIT INPUT/OUTPUT VECTORS===============================================================*/
    //ubft uses blake3 for datatype - uint8_t
    // This implementation is tested on vector of uint8_t of size 4
    std::vector<uint8_t, aligned_allocator<uint8_t> > input(4);
    std::vector<uint8_t, aligned_allocator<uint8_t> > hash_hw(64);
    uint8_t *hash_sw = (uint8_t*) malloc(sizeof(uint8_t) * 1);

    /*====================================================SW VERIFICATION===============================================================*/

    /* Add software Blake*/

    // Initialize the Blake3 hasher
    blake3_hasher hasher = blake3_hasher_init();

    // The data to hash
    const uint8_t data[4] = {0xAB, 0xCD, 0xEF, 0x01};
    input = {0xAB, 0xCD, 0xEF, 0x01};
    printf("\n");

    // add input data to hasher
    hasher = blake3_hasher_update(hasher, (const uint8_t *)data, sizeof(data));
    blake3_hash_output hash = blake3_hasher_finalize(hasher);

    printf("hash value from host: \n");
    // Print the hash in hexadecimal format
    for (size_t i = 0; i < BLAKE3_OUT_LEN; i++) {
        printf("%02x", hash.hash[i]);
    }
    printf("\n");

    /*====================================================Setting up kernel I/O===============================================================*/

    /* INPUT BUFFERS */
    OCL_CHECK(err, cl::Buffer buffer_input(context, CL_MEM_USE_HOST_PTR | CL_MEM_READ_ONLY, sizeof(uint8_t) * 4, input.data(), &err)); 

    /* OUTPUT BUFFERS */
    OCL_CHECK(err, cl::Buffer buffer_output(context, CL_MEM_USE_HOST_PTR | CL_MEM_WRITE_ONLY, sizeof(uint8_t) * 64, hash_hw.data(), &err)); 

    /* SETTING INPUT PARAMETERS */
    OCL_CHECK(err, err = krnl1.setArg(0, buffer_input));
    OCL_CHECK(err, err = krnl1.setArg(1, buffer_output));


    /*====================================================KERNEL===============================================================*/
    /* HOST -> DEVICE DATA TRANSFER*/
    std::cout << "HOST -> DEVICE" << std::endl; 
    htod = clock(); 
    OCL_CHECK(err, err = q.enqueueMigrateMemObjects({buffer_input}, 0 /* 0 means from host*/));
    q.finish();
    htod = clock() - htod; 
    
    /*STARTING KERNEL(S)*/
    std::cout << "STARTING KERNEL(S)" << std::endl; 
    comp = clock(); 
	OCL_CHECK(err, err = q.enqueueTask(krnl1));
    q.finish(); 
    comp = clock() - comp;
    std::cout << "KERNEL(S) FINISHED" << std::endl; 

    /*DEVICE -> HOST DATA TRANSFER*/
    std::cout << "HOST <- DEVICE" << std::endl; 
    dtoh = clock();
    OCL_CHECK(err, err = q.enqueueMigrateMemObjects({buffer_output}, CL_MIGRATE_MEM_OBJECT_HOST));
    q.finish();
    dtoh = clock() - dtoh;

    /*====================================================VERIFICATION & TIMING===============================================================*/

    // printf("Host -> Device : %lf ms\n", 1000.0 * htod/CLOCKS_PER_SEC);
    // printf("Device -> Host : %lf ms\n", 1000.0 * dtoh/CLOCKS_PER_SEC);
    // printf("Computation : %lf ms\n",  1000.0 * comp/CLOCKS_PER_SEC);
    
    bool match = true;

    printf("hash value from kernel: \n");
    for (size_t i = 0; i < BLAKE3_OUT_LEN; i++) {
        printf("%02x", hash_hw[i]);
    }
    printf("\n");

    for (size_t i = 0; i < BLAKE3_OUT_LEN; i++) {
        if (hash_hw[i] != hash.hash[i]){
            match = false;
            break;
        }
    }

    if (match){
        printf("PASSED - Hash value from kernel and host match");
    } else {
        printf("FAILED - Hash value from kernel and host do not match");
    }
    printf("\n");


    free(hash_sw);
    // std::cout << "TEST " << (match ? "PASSED" : "FAILED") <<  std::endl;
    // return (match ? EXIT_SUCCESS : EXIT_FAILURE);
    
}

// blake3 function implementations

static uint32_t IV[8] = {
    0x6A09E667, 0xBB67AE85, 0x3C6EF372, 0xA54FF53A,
    0x510E527F, 0x9B05688C, 0x1F83D9AB, 0x5BE0CD19,
};

static size_t MSG_PERMUTATION[16] = {2, 6,  3,  10, 7, 0,  4,  13,
                                    1, 11, 12, 5,  9, 14, 15, 8};

inline static uint32_t rotate_right(uint32_t x, int n) {
    return (x >> n) | (x << (32 - n));
}

// debugger function to print state of blake3 hasher object
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
    for(int i = 0; i < 64; i++) { 
        printf("%02x", hasher.chunk_state.block[i]);
        if(i % 16 == 15) printf("\n"); 
    }

    printf("key_words: ");
    for(int i = 0; i < 8; i++) {
        printf("%u ", hasher.key_words[i]);
    }
    printf("\n");

    printf("cv_stack: ");
    for(int i = 0; i < 8 * 54; i++) {
        printf("%u ", hasher.cv_stack[i]);
        if(i % 8 == 7) printf("| ");
    }
    printf("\n");

    printf("cv_stack_len: %u\n", hasher.cv_stack_len);
    printf("flags: %u\n", hasher.flags);
}

// The mixing function, G, which mixes either a column or a diagonal.
inline static void g(uint32_t state[16], size_t a, size_t b, size_t c, size_t d,
                    uint32_t mx, uint32_t my) {
    state[a] = state[a] + state[b] + mx;
    state[d] = rotate_right(state[d] ^ state[a], 16);
    state[c] = state[c] + state[d];
    state[b] = rotate_right(state[b] ^ state[c], 12);
    state[a] = state[a] + state[b] + my;
    state[d] = rotate_right(state[d] ^ state[a], 8);
    state[c] = state[c] + state[d];
    state[b] = rotate_right(state[b] ^ state[c], 7);
}

inline static void round_function(uint32_t state[16], uint32_t m[16]) {
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

inline static void permute(uint32_t m[16]) {
    uint32_t permuted[16];
    for (size_t i = 0; i < 16; i++) {
        permuted[i] = m[MSG_PERMUTATION[i]];
    }
    for (size_t i = 0; i < 16; i++) {
        m[i] = permuted[i];
    }
}

inline static void compress(const uint32_t chaining_value[8],
                            const uint32_t block_words[16], uint64_t counter,
                            uint32_t block_len, uint32_t flags,
                            uint32_t out[16]) {
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

    
    for (size_t i = 0; i < 16; i++) {
        out[i] = state[i];
    }

}

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


inline static void output_chaining_value(output self, uint32_t out[8]) {
    uint32_t out16[16];
    compress(self.input_chaining_value, self.block_words, self.counter,
            self.block_len, self.flags, out16);
    for (size_t i = 0; i < 8; i++) {
        out[i] = out16[i];
    }
}

inline static void output_root_bytes(output self, uint8_t out[], size_t out_len) {
    uint64_t output_block_counter = 0;

    size_t out_index = 0; 
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

inline static _blake3_chunk_state chunk_state_init(const uint32_t key_words[8],
                                                   uint64_t chunk_counter, uint32_t flags) {
    _blake3_chunk_state new_state;

    for (size_t i = 0; i < 8; i++) {
        new_state.chaining_value[i] = key_words[i];
    }
    new_state.chunk_counter = chunk_counter;

    for (size_t i = 0; i < BLAKE3_BLOCK_LEN; i++) {
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

inline static uint32_t chunk_state_start_flag(_blake3_chunk_state self) {
    if (self.blocks_compressed == 0) {
        return CHUNK_START;
    } else {
        return 0;
    }
}

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
        self.block_len += (uint8_t)take;
        input_u8 += take;
        input_len -= take;
    }
    return self;
}

inline static output chunk_state_output(_blake3_chunk_state self) { 
    output ret;
    for (size_t i = 0; i < 8; i++) {
        ret.input_chaining_value[i] = self.chaining_value[i];
    }
    words_from_little_endian_bytes(self.block, sizeof(self.block), ret.block_words);
    ret.counter = self.chunk_counter;
    ret.block_len = (uint32_t)self.block_len;
    ret.flags = self.flags | chunk_state_start_flag(self) | CHUNK_END;
    return ret;
}

inline static output parent_output(const uint32_t left_child_cv[8],
                                const uint32_t right_child_cv[8],
                                const uint32_t key_words[8],
                                uint32_t flags) {    
    output ret;
    for (size_t i = 0; i < 8; i++) {
        ret.input_chaining_value[i] = key_words[i];
    }
    for (size_t i = 0; i < 8; i++) {
        ret.block_words[i] = left_child_cv[i];
    }
    for (size_t i = 0; i < 8; i++) {
        ret.block_words[8 + i] = right_child_cv[i];
    }
    ret.counter = 0; // Always 0 for parent nodes.
    ret.block_len =
        BLAKE3_BLOCK_LEN; // Always BLAKE3_BLOCK_LEN (64) for parent nodes.
    ret.flags = PARENT | flags;
    return ret;
}

inline static void parent_cv(const uint32_t left_child_cv[8],
                            const uint32_t right_child_cv[8],
                            const uint32_t key_words[8], uint32_t flags,
                            uint32_t out[8]) {
    output o = parent_output(left_child_cv, right_child_cv, key_words, flags);
    // We only write to `out` after we've read the inputs. That makes it safe for
    // `out` to alias an input, which we do below.
    output_chaining_value(o, out);
}

inline static blake3_hasher hasher_init_internal(const uint32_t key_words[8], uint32_t flags) {
    blake3_hasher self;
    self.chunk_state = chunk_state_init(key_words, 0, flags);
    
    for (size_t i = 0; i < 8; i++) {
        self.key_words[i] = key_words[i];
    }
    self.cv_stack_len = 0;
    self.flags = flags;
    return self;
}

blake3_hasher blake3_hasher_init() {
    blake3_hasher self = hasher_init_internal(IV, 0);
    return self;
}

blake3_hasher hasher_push_stack(blake3_hasher self, const uint32_t cv[8]) {
    size_t start_idx = self.cv_stack_len * 8;
    for (size_t i = 0; i < 8; i++) {
        self.cv_stack[start_idx + i] = cv[i];
    }
    self.cv_stack_len++;
    return self;
}

hasher_pop_result hasher_pop_stack(blake3_hasher self) {
    self.cv_stack_len--;
    hasher_pop_result result;
    result.updated_hasher = self;
    size_t start_idx = self.cv_stack_len * 8;
    for (size_t i = 0; i < 8; i++) {
        result.popped_cv[i] = self.cv_stack[start_idx + i];
    }
    return result;
}

// Section 5.1.2 of the BLAKE3 spec explains this algorithm in more detail.
blake3_hasher hasher_add_chunk_cv(blake3_hasher self, uint32_t new_cv[8], uint64_t total_chunks) {
    // This chunk might complete some subtrees. For each completed subtree, its
    // left child will be the current top entry in the CV stack, and its right
    // child will be the current value of `new_cv`. Pop each left child off the
    // stack, merge it with `new_cv`, and overwrite `new_cv` with the result.
    // After all these merges, push the final value of `new_cv` onto the stack.
    // The number of completed subtrees is given by the number of trailing 0-bits
    // in the new total number of chunks.
    while ((total_chunks & 1) == 0) {
        hasher_pop_result pop_result = hasher_pop_stack(self);
        self = pop_result.updated_hasher;
        parent_cv(pop_result.popped_cv, new_cv, self.key_words, self.flags, new_cv);
        total_chunks >>= 1;
    }
    self = hasher_push_stack(self, new_cv);
    return self;
}

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