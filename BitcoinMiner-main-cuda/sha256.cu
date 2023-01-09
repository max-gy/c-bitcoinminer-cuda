#include <iostream>
#include <stdio.h>

#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <dirent.h>
#include <ctype.h>


#include "/usr/local/cuda-11.4/targets/aarch64-linux/include/cuda_runtime.h"
#include "/usr/local/cuda-11.4/targets/aarch64-linux/include/device_launch_parameters.h"



uint32_t rotateInt(uint32_t inputWord, int numberOfBitsToRotate) 
{
    int bitWidth = sizeof(inputWord) * 8;
    // Rotating 32 bits on a 32-bit integer is the same as rotating 0 bits;
    //   33 bits -> 1 bit; etc.
    numberOfBitsToRotate = numberOfBitsToRotate % bitWidth;

    uint32_t tempWord = inputWord;

    // Rotate input to the right
    inputWord = inputWord >> numberOfBitsToRotate;

    // Build mask for carried over bits
    tempWord = tempWord << (bitWidth - numberOfBitsToRotate);

    return inputWord | tempWord;
}

int Ch(int x, int y, int z)
{
    return ((x & y) ^ (~x & z));
}

int Maj(int x, int y, int z)
{
    return ((x & y) ^ (x & z) ^ (y & z));
}

int Sig0f(int x)
{
    return(rotateInt(x, 2) ^ rotateInt(x, 13) ^ rotateInt(x, 22));
}

int Sig1f(int x)
{
    return(rotateInt(x, 6) ^ rotateInt(x, 11) ^ rotateInt(x,25));
}

uint32_t sig0(uint32_t x)
{
    return(rotateInt(x, 7) ^ rotateInt(x, 18) ^ (x >> 3));
}

uint32_t sig1(uint32_t x)
{
    return(rotateInt(x, 17) ^ rotateInt(x, 19) ^ (x >> 10));
}

void hash(uint32_t *input, int bitlength, uint32_t *outputlocation)
{
    uint32_t H_0[8] = { 0x6a09e667, 0xbb67ae85, 0x3c6ef372, 0xa54ff53a, 0x510e527f, 0x9b05688c, 0x1f83d9ab, 0x5be0cd19 };

    uint32_t K[64] = {0x428a2f98, 0x71374491, 0xb5c0fbcf, 0xe9b5dba5,
    0x3956c25b, 0x59f111f1, 0x923f82a4, 0xab1c5ed5,
    0xd807aa98, 0x12835b01, 0x243185be, 0x550c7dc3,
    0x72be5d74, 0x80deb1fe, 0x9bdc06a7, 0xc19bf174,
    0xe49b69c1, 0xefbe4786, 0x0fc19dc6, 0x240ca1cc,
    0x2de92c6f, 0x4a7484aa, 0x5cb0a9dc, 0x76f988da,
    0x983e5152, 0xa831c66d, 0xb00327c8, 0xbf597fc7,
    0xc6e00bf3, 0xd5a79147, 0x06ca6351, 0x14292967,
    0x27b70a85, 0x2e1b2138, 0x4d2c6dfc, 0x53380d13,
    0x650a7354, 0x766a0abb, 0x81c2c92e, 0x92722c85,
    0xa2bfe8a1, 0xa81a664b, 0xc24b8b70, 0xc76c51a3,
    0xd192e819, 0xd6990624, 0xf40e3585, 0x106aa070,
    0x19a4c116, 0x1e376c08, 0x2748774c, 0x34b0bcb5,
    0x391c0cb3, 0x4ed8aa4a, 0x5b9cca4f, 0x682e6ff3,
    0x748f82ee, 0x78a5636f, 0x84c87814, 0x8cc70208,
    0x90befffa, 0xa4506ceb, 0xbef9a3f7, 0xc67178f2};

    int wordlength = bitlength / 32 + 1;
    int k = (512 * 512 - bitlength - 1) % 512;
    uint32_t message[10000] = {};

    for(int i = 0; i < wordlength; i++)
        message[i] = input[i];

    if(bitlength % 32 != 0)
        message[bitlength / 32] = message[bitlength / 32] | (1 << (32 - 1 - bitlength % 32));
    else
        message[bitlength / 32] = 1 << 31;
    
    uint32_t rounds;

    // Assuming our data isn't bigger than 2^32 bits long... which it won't be for a block hash.
    if(wordlength % 16 == 0 || wordlength % 16 == 15)
    {
        message[wordlength + 15 + 16 - wordlength % 16] = bitlength;
        rounds = wordlength / 16 + 2;
    }
    else
    {
        message[wordlength + 15 - wordlength % 16] = bitlength;
        rounds = wordlength / 16 + 1;
    }
        
    uint32_t M[32][16];

    for(int i = 0; i < 16; i++)
        for(int j = 0; j <= rounds; j++)
            M[j][i] = message[i + j * 16];
    
    uint32_t H[32][8];

    for(int i = 0; i < 8; i++)
        H[0][i] = H_0[i];

    // Here our hash function rounds actually start.
    for(int i = 1; i <= rounds; i++)
    {
        uint32_t a = H[i-1][0];
        uint32_t b = H[i-1][1];
        uint32_t c = H[i-1][2];
        uint32_t d = H[i-1][3];
        uint32_t e = H[i-1][4];
        uint32_t f = H[i-1][5];
        uint32_t g = H[i-1][6];
        uint32_t h = H[i-1][7];

        uint32_t W[64];

        for(int j = 0; j < 64; j++)
        {
            uint32_t ch = Ch(e, f, g);
            uint32_t maj = Maj(a, b, c);
            uint32_t Sig0 = Sig0f(a);
            uint32_t Sig1 = Sig1f(e);

            if(j < 16)
                W[j] = M[i-1][j];
            else
                W[j] = sig1(W[j-2]) + W[j-7] + sig0(W[j-15]) + W[j-16];
            
            uint32_t T1 = h + Sig1 + ch + K[j] + W[j];
            uint32_t T2 = Sig0 + maj;
            h = g;
            g = f;
            f = e;
            e = d + T1;
            d = c;
            c = b;
            b = a;
            a = T1 + T2;
        }

        H[i][0] = a + H[i-1][0];
        H[i][1] = b + H[i-1][1];
        H[i][2] = c + H[i-1][2];
        H[i][3] = d + H[i-1][3];
        H[i][4] = e + H[i-1][4];
        H[i][5] = f + H[i-1][5];
        H[i][6] = g + H[i-1][6];
        H[i][7] = h + H[i-1][7];
    }

    for(int i = 0; i < 8; i++)
        outputlocation[i] = H[rounds][i];
}

// CUDA PART

#define checkCudaErrors(x) \
{ \
    cudaGetLastError(); \
    x; \
    cudaError_t err = cudaGetLastError(); \
    if (err != cudaSuccess) \
        printf("GPU: cudaError %d (%s)\n", err, cudaGetErrorString(err)); \
}

#define ROTLEFT(a,b) (((a) << (b)) | ((a) >> (32-(b))))
#define ROTRIGHT(a,b) (((a) >> (b)) | ((a) << (32-(b))))

#define CH(x,y,z) (((x) & (y)) ^ (~(x) & (z)))
#define MAJ(x,y,z) (((x) & (y)) ^ ((x) & (z)) ^ ((y) & (z)))
#define EP0(x) (ROTRIGHT(x,2) ^ ROTRIGHT(x,13) ^ ROTRIGHT(x,22))
#define EP1(x) (ROTRIGHT(x,6) ^ ROTRIGHT(x,11) ^ ROTRIGHT(x,25))
#define SIG0(x) (ROTRIGHT(x,7) ^ ROTRIGHT(x,18) ^ ((x) >> 3))
#define SIG1(x) (ROTRIGHT(x,17) ^ ROTRIGHT(x,19) ^ ((x) >> 10))


typedef unsigned char BYTE;             // 8-bit byte
typedef uint32_t  WORD;             // 32-bit word, change to "long" for 16-bit machines

typedef struct JOB {
	BYTE * data;
	unsigned long long size;
	BYTE digest[64];
	char fname[128];
}JOB;


typedef struct {
	BYTE data[64];
	WORD datalen;
	unsigned long long bitlen;
	WORD state[8];
} SHA256_CTX;

__constant__ WORD dev_k[64];

static const WORD host_k[64] = {
	0x428a2f98,0x71374491,0xb5c0fbcf,0xe9b5dba5,0x3956c25b,0x59f111f1,0x923f82a4,0xab1c5ed5,
	0xd807aa98,0x12835b01,0x243185be,0x550c7dc3,0x72be5d74,0x80deb1fe,0x9bdc06a7,0xc19bf174,
	0xe49b69c1,0xefbe4786,0x0fc19dc6,0x240ca1cc,0x2de92c6f,0x4a7484aa,0x5cb0a9dc,0x76f988da,
	0x983e5152,0xa831c66d,0xb00327c8,0xbf597fc7,0xc6e00bf3,0xd5a79147,0x06ca6351,0x14292967,
	0x27b70a85,0x2e1b2138,0x4d2c6dfc,0x53380d13,0x650a7354,0x766a0abb,0x81c2c92e,0x92722c85,
	0xa2bfe8a1,0xa81a664b,0xc24b8b70,0xc76c51a3,0xd192e819,0xd6990624,0xf40e3585,0x106aa070,
	0x19a4c116,0x1e376c08,0x2748774c,0x34b0bcb5,0x391c0cb3,0x4ed8aa4a,0x5b9cca4f,0x682e6ff3,
	0x748f82ee,0x78a5636f,0x84c87814,0x8cc70208,0x90befffa,0xa4506ceb,0xbef9a3f7,0xc67178f2
};

/*********************** FUNCTION DECLARATIONS **********************/
char * print_sha(BYTE * buff);
__device__ void sha256_init(SHA256_CTX *ctx);
__device__ void sha256_update(SHA256_CTX *ctx, const BYTE data[], size_t len);
__device__ void sha256_final(SHA256_CTX *ctx, BYTE hash[]);


char * hash_to_string(BYTE * buff) {
	char * string = (char *)malloc(70);
	int k, i;
	for (i = 0, k = 0; i < 32; i++, k+= 2)
	{
		sprintf(string + k, "%.2x", buff[i]);
		//printf("%02x", buff[i]);
	}
	string[64] = 0;
	return string;
}

void print_job(JOB * j){
	printf("%s  %s\n", hash_to_string(j->digest), j->fname);
}

void print_jobs(JOB ** jobs, int n) {
	for (int i = 0; i < n; i++)
	{
        print_job(jobs[i]);
		// printf("@ %p JOB[%i] \n", jobs[i], i);
		// printf("\t @ 0x%p data = %x \n", jobs[i]->data, (jobs[i]->data == 0)? 0 : jobs[i]->data[0]);
		// printf("\t @ 0x%p size = %llu \n", &(jobs[i]->size), jobs[i]->size);
		// printf("\t @ 0x%p fname = %s \n", &(jobs[i]->fname), jobs[i]->fname);
		// printf("\t @ 0x%p digest = %s \n------\n", jobs[i]->digest, hash_to_string(jobs[i]->digest));
	}
}

__device__ void mycpy12(uint32_t *d, const uint32_t *s) {
#pragma unroll 3
    for (int k=0; k < 3; k++) d[k] = s[k];
}

__device__ void mycpy16(uint32_t *d, const uint32_t *s) {
#pragma unroll 4
    for (int k=0; k < 4; k++) d[k] = s[k];
}

__device__ void mycpy32(uint32_t *d, const uint32_t *s) {
#pragma unroll 8
    for (int k=0; k < 8; k++) d[k] = s[k];
}

__device__ void mycpy44(uint32_t *d, const uint32_t *s) {
#pragma unroll 11
    for (int k=0; k < 11; k++) d[k] = s[k];
}

__device__ void mycpy48(uint32_t *d, const uint32_t *s) {
#pragma unroll 12
    for (int k=0; k < 12; k++) d[k] = s[k];
}

__device__ void mycpy64(uint32_t *d, const uint32_t *s) {
#pragma unroll 16
    for (int k=0; k < 16; k++) d[k] = s[k];
}

__device__ void sha256_transform(SHA256_CTX *ctx, const BYTE data[])
{
	WORD a, b, c, d, e, f, g, h, i, j, t1, t2, m[64];
    WORD S[8];

    //mycpy32(S, ctx->state);

    #pragma unroll 16
	for (i = 0, j = 0; i < 16; ++i, j += 4)
		m[i] = (data[j] << 24) | (data[j + 1] << 16) | (data[j + 2] << 8) | (data[j + 3]);

    #pragma unroll 64
	for (; i < 64; ++i)
		m[i] = SIG1(m[i - 2]) + m[i - 7] + SIG0(m[i - 15]) + m[i - 16];

	a = ctx->state[0];
	b = ctx->state[1];
	c = ctx->state[2];
	d = ctx->state[3];
	e = ctx->state[4];
	f = ctx->state[5];
	g = ctx->state[6];
	h = ctx->state[7];

    #pragma unroll 64
	for (i = 0; i < 64; ++i) {
		t1 = h + EP1(e) + CH(e, f, g) + dev_k[i] + m[i];
		t2 = EP0(a) + MAJ(a, b, c);
		h = g;
		g = f;
		f = e;
		e = d + t1;
		d = c;
		c = b;
		b = a;
		a = t1 + t2;
	}

	ctx->state[0] += a;
	ctx->state[1] += b;
	ctx->state[2] += c;
	ctx->state[3] += d;
	ctx->state[4] += e;
	ctx->state[5] += f;
	ctx->state[6] += g;
	ctx->state[7] += h;
}

__device__ void sha256_init(SHA256_CTX *ctx)
{
	ctx->datalen = 0;
	ctx->bitlen = 0;
	ctx->state[0] = 0x6a09e667;
	ctx->state[1] = 0xbb67ae85;
	ctx->state[2] = 0x3c6ef372;
	ctx->state[3] = 0xa54ff53a;
	ctx->state[4] = 0x510e527f;
	ctx->state[5] = 0x9b05688c;
	ctx->state[6] = 0x1f83d9ab;
	ctx->state[7] = 0x5be0cd19;
}

__device__ void sha256_update(SHA256_CTX *ctx, const BYTE data[], size_t len)
{
	WORD i;

	// for each byte in message
	for (i = 0; i < len; ++i) {
		// ctx->data == message 512 bit chunk
		ctx->data[ctx->datalen] = data[i];
		ctx->datalen++;
		if (ctx->datalen == 64) {
			sha256_transform(ctx, ctx->data);
			ctx->bitlen += 512;
			ctx->datalen = 0;
		}
	}
}

__device__ void sha256_final(SHA256_CTX *ctx, BYTE hash[])
{
	WORD i;

	i = ctx->datalen;

	// Pad whatever data is left in the buffer.
	if (ctx->datalen < 56) {
		ctx->data[i++] = 0x80;
		while (i < 56)
			ctx->data[i++] = 0x00;
	}
	else {
		ctx->data[i++] = 0x80;
		while (i < 64)
			ctx->data[i++] = 0x00;
		sha256_transform(ctx, ctx->data);
		memset(ctx->data, 0, 56);
	}

	// Append to the padding the total message's length in bits and transform.
	ctx->bitlen += ctx->datalen * 8;
	ctx->data[63] = ctx->bitlen;
	ctx->data[62] = ctx->bitlen >> 8;
	ctx->data[61] = ctx->bitlen >> 16;
	ctx->data[60] = ctx->bitlen >> 24;
	ctx->data[59] = ctx->bitlen >> 32;
	ctx->data[58] = ctx->bitlen >> 40;
	ctx->data[57] = ctx->bitlen >> 48;
	ctx->data[56] = ctx->bitlen >> 56;
	sha256_transform(ctx, ctx->data);

	// Since this implementation uses little endian byte ordering and SHA uses big endian,
	// reverse all the bytes when copying the final state to the output hash.
	for (i = 0; i < 4; ++i) {
		hash[i] = (ctx->state[0] >> (24 - i * 8)) & 0x000000ff;
		hash[i + 4] = (ctx->state[1] >> (24 - i * 8)) & 0x000000ff;
		hash[i + 8] = (ctx->state[2] >> (24 - i * 8)) & 0x000000ff;
		hash[i + 12] = (ctx->state[3] >> (24 - i * 8)) & 0x000000ff;
		hash[i + 16] = (ctx->state[4] >> (24 - i * 8)) & 0x000000ff;
		hash[i + 20] = (ctx->state[5] >> (24 - i * 8)) & 0x000000ff;
		hash[i + 24] = (ctx->state[6] >> (24 - i * 8)) & 0x000000ff;
		hash[i + 28] = (ctx->state[7] >> (24 - i * 8)) & 0x000000ff;
	}

}

__global__ void sha256_cuda(JOB ** jobs, int n) {
	int i = blockIdx.x * blockDim.x + threadIdx.x;
	// perform sha256 calculation here
	if (i < n){
		SHA256_CTX ctx;
		sha256_init(&ctx);
		sha256_update(&ctx, jobs[i]->data, jobs[i]->size);
		sha256_final(&ctx, jobs[i]->digest);
	}
}



void runJobs(JOB ** jobs, int n){
	int blockSize = 4;
	int numBlocks = (n + blockSize - 1) / blockSize;
	sha256_cuda <<< numBlocks, blockSize >>> (jobs, n);
}


JOB * JOB_init(BYTE * data, long size, char * fname) {
	JOB * j;
	checkCudaErrors(cudaMallocManaged(&j, sizeof(JOB)));	//j = (JOB *)malloc(sizeof(JOB));
	checkCudaErrors(cudaMallocManaged(&(j->data), size));
	j->data = data;
	j->size = size;
	for (int i = 0; i < 64; i++)
	{
		j->digest[i] = 0xff;
	}
	std::cout << "JOB_init inside" << std::endl;
	strcpy(j->fname, fname);
	return j;
}

unsigned char *uint32_to_char_array(uint32_t n)
{
    unsigned char *a;

   // a = wrap_calloc(4, sizeof(unsigned char));

    a[0] = (n >> 24) & 0xff;  /* high-order (leftmost) byte: bits 24-31 */
    a[1] = (n >> 16) & 0xff;  /* next byte, counting from left: bits 16-23 */
    a[2] = (n >>  8) & 0xff;  /* next byte, bits 8-15 */
    a[3] = n         & 0xff;  /* low-order byte: bits 0-7 */

    return a;
}

void hash_CUDA (uint32_t input, unsigned long bitlength, uint32_t * outputlocation) {
    //    checkCudaErrors(cudaMallocManaged(&jobs, n * sizeof(JOB *)));

		// iterate over file list - non optional arguments
			BYTE * inputBytes[4];
			char * outputlocationBytes[4];

			//uint32_to_char_array(input);
			uint8_t byteval[4];
			
			for(int i = 0; i < 4; i++) byteval[i] = input >> (i*8);
			for(int i = 0; i < 4; i++) byteval[i] = input >> (i*8);


			std::cout << "JOB_init 0" << input << std::endl;

			//for (int i=0; i<4 ;++i)
			//	inputBytes[i] = &((BYTE*)*input)[3-i];

			JOB * jobs[4];

			std::cout << "JOB_init 1" << byteval[0] << std::endl;

			jobs[0] = JOB_init(&byteval[0], bitlength, outputlocationBytes[0]);

			std::cout << "JOB_init 2" << std::endl;

			jobs[1] = JOB_init(&byteval[1], bitlength, outputlocationBytes[1]);
			std::cout << "JOB_init 3" << std::endl;

			jobs[2] = JOB_init(&byteval[2], bitlength, outputlocationBytes[2]);
			std::cout << "JOB_init 4" << std::endl;

			jobs[3] = JOB_init(&byteval[3], bitlength, outputlocationBytes[3]);

			std::cout << "JOB_init after" << std::endl;

			
			outputlocation = ((uint32_t*) &outputlocationBytes);

			std::cout << "Before runJobs " << std::endl;

//        }

       // pre_sha256();
        runJobs(jobs, 1);
  //  }
}
   	

