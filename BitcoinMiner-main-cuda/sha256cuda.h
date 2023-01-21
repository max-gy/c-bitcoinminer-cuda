#ifndef SHA256CUDA_H
#define SHA256CUDA_H


/****************************** MACROS ******************************/
/*
#define SHA256_BLOCK_SIZE 32            // SHA256 outputs a 32 byte digest

#define ROTLEFT(a,b) (((a) << (b)) | ((a) >> (32-(b))))
#define ROTRIGHT(a,b) (((a) >> (b)) | ((a) << (32-(b))))

#define CH(x,y,z) (((x) & (y)) ^ (~(x) & (z)))
#define MAJ(x,y,z) (((x) & (y)) ^ ((x) & (z)) ^ ((y) & (z)))
#define EP0(x) (ROTRIGHT(x,2) ^ ROTRIGHT(x,13) ^ ROTRIGHT(x,22))
#define EP1(x) (ROTRIGHT(x,6) ^ ROTRIGHT(x,11) ^ ROTRIGHT(x,25))
#define SIG0(x) (ROTRIGHT(x,7) ^ ROTRIGHT(x,18) ^ ((x) >> 3))
#define SIG1(x) (ROTRIGHT(x,17) ^ ROTRIGHT(x,19) ^ ((x) >> 10))

*/
#define checkCudaErrors(x) \
{ \
    cudaGetLastError(); \
    x; \
    cudaError_t err = cudaGetLastError(); \
    if (err != cudaSuccess) \
        printf("GPU: cudaError %d (%s)\n", err, cudaGetErrorString(err)); \
}
/**************************** DATA TYPES ****************************/
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


__constant__ uint32_t dev_H[32][8];
__constant__ uint32_t dev_M[32][16];


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

void print_job(JOB  j){
	//printf("%s  %s\n", hash_to_string(j->digest), j->fname);
    printf("%s  \n", (j.data));

}

void print_jobs(JOB * jobs, int n) {
	printf("\n");
	for (int i = 0; i < n; i++)
	{
        //print_job(jobs[i]);
		// printf("@ %p JOB[%i] \n", jobs[i], i);
		 //printf("\t @ 0x%p data = %x \n", jobs[i]->data, (jobs[i]->data == 0)? 0 : jobs[i]->data[0]);
		// printf("\t @ 0x%p size = %llu \n", &(jobs[i]->size), jobs[i]->size);
		// printf("\t @ 0x%p fname = %s \n", &(jobs[i]->fname), jobs[i]->fname);
		// printf("\t @ 0x%p digest = %s \n------\n", jobs[i]->digest, hash_to_string(jobs[i]->digest));
	}
}


__device__ uint32_t rotateInt(uint32_t inputWord, int numberOfBitsToRotate) 
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

__device__ int Ch(int x, int y, int z)
{
    return ((x & y) ^ (~x & z));
}

__device__  int Maj(int x, int y, int z)
{
    return ((x & y) ^ (x & z) ^ (y & z));
}

__device__  int Sig0f(int x)
{
    return(rotateInt(x, 2) ^ rotateInt(x, 13) ^ rotateInt(x, 22));
}

__device__  int Sig1f(int x)
{
    return(rotateInt(x, 6) ^ rotateInt(x, 11) ^ rotateInt(x,25));
}

__device__ uint32_t sig0(uint32_t x)
{
    return(rotateInt(x, 7) ^ rotateInt(x, 18) ^ (x >> 3));
}

__device__ uint32_t sig1(uint32_t x)
{
    return(rotateInt(x, 17) ^ rotateInt(x, 19) ^ (x >> 10));
}

/*
void cudahash(uint32_t H[][8], uint32_t M[][16], uint32_t rounds, uint32_t *outputlocation)
{
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
*/

#endif   // SHA256_H
