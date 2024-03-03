// cd /home/hork/cuda-workspace/CudaSHA256/Debug/files
// time ~/Dropbox/FIIT/APS/Projekt/CpuSHA256/a.out -f ../file-list
// time ../CudaSHA256 -f ../file-list

#include <iostream>

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
//#include <cuda_runtime.h>
#include "sha256cuda.h"
#include "util.h"
#include "miner.h"
#include <dirent.h>
#include <ctype.h>
#include <chrono>

#include "/usr/local/cuda-11/targets/aarch64-linux/include/cuda_runtime.h"
#include "/usr/local/cuda-11/targets/aarch64-linux/include/device_launch_parameters.h"


__global__ void sha256_cuda_try(char* a, char* b) {
	int i = blockIdx.x * blockDim.x + threadIdx.x;
	// perform sha256 calculation here
	b[i] = a[i];

	return;
}


void pre_sha256() {
	// compy symbols
//	checkCudaErrors(cudaMemcpyToSymbol(dev_k, host_k, sizeof(host_k), 0, cudaMemcpyHostToDevice));
}

__global__ void sha256_cuda_hash(uint32_t* result_H) {

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

	int index = blockIdx.x * blockDim.x + threadIdx.x;
	// perform sha256 calculation here

	
	/*
	result_H[i*8+0] = i*8+0;
	result_H[i*8+1] = i*8+1;
	result_H[i*8+2] = i*8+2;
	result_H[i*8+3] = i*8+3;
	result_H[i*8+4] = i*8+4;
	result_H[i*8+5] = i*8+5;
	result_H[i*8+6] = i*8+6;
	result_H[i*8+7] = i*8+7;
*/

	int i = index+1;
	if (i>0) { 
		
		uint32_t a = dev_H[i-1][0];
		uint32_t b = dev_H[i-1][1];
		uint32_t c = dev_H[i-1][2];
		uint32_t d = dev_H[i-1][3];
		uint32_t e = dev_H[i-1][4];
		uint32_t f = dev_H[i-1][5];
		uint32_t g = dev_H[i-1][6];
		uint32_t h = dev_H[i-1][7];

		
		uint32_t W[64];

		for(int j = 0; j < 64; j++)
		{
			uint32_t ch = Ch(e, f, g);
			uint32_t maj = Maj(a, b, c);
			uint32_t Sig0 = Sig0f(a);
			uint32_t Sig1 = Sig1f(e);

			if(j < 16)
				W[j] = dev_M[i-1][j];
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

		result_H[i*8+0] = a + dev_H[i-1][0];
		result_H[i*8+1] = b + dev_H[i-1][1];
		result_H[i*8+2] = c + dev_H[i-1][2];
		result_H[i*8+3] = d + dev_H[i-1][3];
		result_H[i*8+4] = e + dev_H[i-1][4];
		result_H[i*8+5] = f + dev_H[i-1][5];
		result_H[i*8+6] = g + dev_H[i-1][6];
		result_H[i*8+7] = h + dev_H[i-1][7];


	}

	//outputlocation[i] = outputlocation[i]*outputlocation[i];

	return;
}

__global__ void sha256_cuda_hash_full(uint64_t * nonce_found) {

	//int index = blockIdx.x * blockDim.x + threadIdx.x;


    //calculates unique thread ID in the block
    int t= (blockDim.x*blockDim.y)*threadIdx.z
		   +(threadIdx.y*blockDim.x)+(threadIdx.x);
	//calculates unique block ID in the grid
    int b= (gridDim.x*gridDim.y)*blockIdx.z
		   +(blockIdx.y*gridDim.x)+(blockIdx.x);
	//block size (this is redundant though) 
	int T= blockDim.x*blockDim.y*blockDim.z;
	//grid size (this is redundant though)
	int B= gridDim.x*gridDim.y*gridDim.z;
	

	int index= b * T + t;

	uint32_t H_0[8] = { 0x6a09e667, 0xbb67ae85, 0x3c6ef372, 0xa54ff53a, 0x510e527f, 0x9b05688c, 0x1f83d9ab, 0x5be0cd19 };

	//result_H[index] = index;


	uint32_t blockheader[20];
	for (int bi=0;bi<19;bi++) blockheader[bi] = dev_blockheader[bi];
	blockheader[19] = dev_nonce[0] + index;
    //*(blockheader + 19) = dev_nonce[0] + index;

    //print_bytes((unsigned char*)blockheader, 80);



    for(int i = 0; i < 20; i++)
        blockheader[i] = __Reverse32(blockheader[i]);

	
    uint32_t hash0[8];

	uint32_t outputlocation[8];
	uint32_t H[32][8];

	for (int hash_round = 1;hash_round<3;hash_round++) 
	{
		
		
		int bitlength = 640;
		uint32_t* input = blockheader;
		if (hash_round == 2) {
			bitlength = 256;
			input = outputlocation;
		}	

		

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
				
				
				uint32_t T1 = h + Sig1 + ch + dev_K[j] + W[j];
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

		for(int oi = 0; oi < 8; oi++)
			outputlocation[oi] = H[rounds][oi];

		
	}

	nonce_found[0] = 0;
	bool solved = false;
	for(int i = 0; i < 8; i++)
	{
		if(outputlocation[7-i] < dev_difficulty[i])
		{
			solved = true;
			nonce_found[0] = index;
		}
		else if(outputlocation[7-i] > dev_difficulty[i])
			break;
		// And if they're equal, we keep going!
	}
	nonce_found[0] = index;
//	result_H[index] = solved;

	return;

 
}


void __hashblock(uint32_t _nonce, char* version, char* prevhash, 
	char* merkle_root, char* time, char* nbits)
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
		
	

	uint32_t nonce[1];
	nonce[0] = _nonce;

	uint32_t difficulty[8];
    uint32_t bits[1];
    hexstr_to_intarray(nbits, bits);
    bits_to_difficulty(*bits, difficulty);

	checkCudaErrors(cudaMemcpyToSymbol(dev_K, K, sizeof(K), 0, cudaMemcpyHostToDevice));
	checkCudaErrors(cudaMemcpyToSymbol(dev_difficulty, difficulty, sizeof(difficulty), 0, cudaMemcpyHostToDevice));

	int solved = 0;
	while (!solved) {

		uint32_t blockheader[20];

		hexstr_to_intarray(version, blockheader);
		hexstr_to_intarray(prevhash, blockheader + 1);
		hexstr_to_intarray(merkle_root, blockheader + 9);
		hexstr_to_intarray(time, blockheader + 17);
		hexstr_to_intarray(nbits, blockheader + 18);
		*(blockheader + 19) = nonce[0];

		const uint64_t SIZE_BUFFER = 1024*1024;
		uint64_t * cuda_result_nonce_1 = 0;
		uint64_t * cuda_result_nonce_2 = 0;

		uint64_t result_nonce[1];


		checkCudaErrors(cudaMemcpyToSymbol(dev_nonce, nonce, sizeof(nonce), 0, cudaMemcpyHostToDevice));
		checkCudaErrors(cudaMemcpyToSymbol(dev_blockheader, blockheader, sizeof(blockheader), 0, cudaMemcpyHostToDevice));

		std::chrono::steady_clock::time_point start = std::chrono::steady_clock::now();


		cudaMalloc(&cuda_result_nonce_1, sizeof(result_nonce));

		sha256_cuda_hash_full <<< 1024, 1024 >>> (cuda_result_nonce_1);

		//cudaMalloc(&cuda_result_nonce_2, sizeof(result_nonce));

		//sha256_cuda_hash_full <<< (64,64,64), 1024 >>> (cuda_result_nonce_2);

		
		cudaMemcpy(result_nonce,cuda_result_nonce_1,sizeof(result_nonce),cudaMemcpyDeviceToHost); 

		nonce[0]+=result_nonce[0];

		std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();
		uint64_t duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
		uint64_t hashrate = result_nonce[0] / ((uint64_t)duration);
		std::cout << "Currently mining at " << hashrate << "000 hashes / second" << std::endl;
			
		
		std::cout << "finished cuda \n";
		std::cout << result_nonce[0] << " ";// << nonce[0] << "\n";
		

	}



}

void run_hash (uint32_t * input, int bitlength, uint32_t * outputlocation) {

	uint32_t H_0[8] = { 0x6a09e667, 0xbb67ae85, 0x3c6ef372, 0xa54ff53a, 0x510e527f, 0x9b05688c, 0x1f83d9ab, 0x5be0cd19 };

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
    

	
	//cudahash(H, M, rounds, outputlocation);

	uint32_t* cuda_result = 0;

	uint32_t* cuda_message= 0;

	//checkCudaErrors(cudaMemcpyToSymbol(dev_message, message, sizeof(message), 0, cudaMemcpyHostToDevice));
	checkCudaErrors(cudaMemcpyToSymbol(dev_H, H, sizeof(H), 0, cudaMemcpyHostToDevice));
	checkCudaErrors(cudaMemcpyToSymbol(dev_M, M, sizeof(M), 0, cudaMemcpyHostToDevice));

	cudaMallocManaged(&cuda_result, 32*8*sizeof(uint32_t));

	cudaMemcpy(cuda_message,message,sizeof(M),cudaMemcpyHostToDevice); 

	sha256_cuda_hash <<< 1, rounds >>> (cuda_result);

	uint32_t __H[32*8];

	cudaMemcpy(__H,cuda_result,sizeof(__H),cudaMemcpyDeviceToHost); 


	for(int i = 0; i < 8; i++)
		__H[i] = H_0[i];
		

	for(int i = 0; i < 8; i++)
		outputlocation[i] = __H[rounds*8+i];

    //std::cout << "\n Run hash input\n";
	//for (int y=0;y<32*8;y++) 
        	//std::cout << __H[y] << " ";
    
}


void runJobs_workingedition(JOB ** jobs, int n){

	char a[4] = "sav";
	char b[sizeof(a)/sizeof(char)];
	char *ca = 0;
	char *cb = 0;

	

	cudaMalloc(&ca, sizeof(a));
	cudaMalloc(&cb, sizeof(b));
	cudaMemcpy(ca, a, sizeof(a), cudaMemcpyHostToDevice);

	
	sha256_cuda_try <<< 1, sizeof(a)/sizeof(char) >>> (ca, cb);

	cudaMemcpy(b, cb, sizeof(b), cudaMemcpyDeviceToHost);

	std::cout << "hhelo";
	std::cout << b;
	std::cout << 10;
}

