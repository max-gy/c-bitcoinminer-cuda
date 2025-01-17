#include <stdint.h>

int buildCUDAJobs(ssize_t hash);
void run_hash (uint32_t * input, int bitlength, uint32_t * outputlocation);
void __hashblock(uint32_t nonce, char* version, char* prevhash, char* merkle_root, char* time, char* nbits);