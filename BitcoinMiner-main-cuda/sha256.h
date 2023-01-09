typedef unsigned char BYTE;
int hash(uint32_t *input, int bitlength, uint32_t *outputlocation);
int hash_CUDA(uint32_t input, unsigned long bitlength, uint32_t *outputlocation);