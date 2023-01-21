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

#include <dirent.h>
#include <ctype.h>

#include "/usr/local/cuda-11.4/targets/aarch64-linux/include/cuda_runtime.h"
#include "/usr/local/cuda-11.4/targets/aarch64-linux/include/device_launch_parameters.h"

char * trim(char *str){
    size_t len = 0;
    char *frontp = str;
    char *endp = NULL;

    if( str == NULL ) { return NULL; }
    if( str[0] == '\0' ) { return str; }

    len = strlen(str);
    endp = str + len;

    /* Move the front and back pointers to address the first non-whitespace
     * characters from each end.
     */
    while( isspace((unsigned char) *frontp) ) { ++frontp; }
    if( endp != frontp )
    {
        while( isspace((unsigned char) *(--endp)) && endp != frontp ) {}
    }

    if( str + len - 1 != endp )
            *(endp + 1) = '\0';
    else if( frontp != str &&  endp == frontp )
            *str = '\0';

    /* Shift the string so that it starts at str so that if it's dynamically
     * allocated, we can still free it on the returned pointer.  Note the reuse
     * of endp to mean the front of the string buffer now.
     */
    endp = str;
    if( frontp != str )
    {
            while( *frontp ) { *endp++ = *frontp++; }
            *endp = '\0';
    }


    return str;
}

__global__ void sha256_cuda(JOB * jobs, JOB * results) {
	int i = blockIdx.x * blockDim.x + threadIdx.x;
	// perform sha256 calculation here

	//results[i].data = jobs[i].data; 
	if (0 < 1){
		SHA256_CTX ctx;
	}

	 //reinterpret_cast<unsigned char*>('2');
	return;
}

__global__ void sha256_cuda_try(char* a, char* b) {
	int i = blockIdx.x * blockDim.x + threadIdx.x;
	// perform sha256 calculation here
	b[i] = a[i];

	return;
}

__global__ void sha256_cuda_hash_640(uint32_t *input, uint32_t *outputlocation) {
	int i = blockIdx.x * blockDim.x + threadIdx.x;
	// perform sha256 calculation here
	cudahash(&input[i], 640, &outputlocation[i]);

	//outputlocation[i] = outputlocation[i]*outputlocation[i];

	return;
}

__global__ void sha256_cuda_hash_256(uint32_t *input, uint32_t *outputlocation) {
	int i = blockIdx.x * blockDim.x + threadIdx.x;
	// perform sha256 calculation here
	cudahash(&input[i], 256, &outputlocation[i]);

	//outputlocation[i] = outputlocation[i]*outputlocation[i];

	return;
}



void pre_sha256() {
	// compy symbols
//	checkCudaErrors(cudaMemcpyToSymbol(dev_k, host_k, sizeof(host_k), 0, cudaMemcpyHostToDevice));
}


void runJobs(JOB ** jobs, int n){
	int blockSize = 4;
	int numBlocks = (n + blockSize - 1) / blockSize;

	JOB* cuda_jobs = 0;
	JOB* cuda_jobs_r = 0;

	BYTE * buff;

	
	JOB js[1];
	char a[9] = "hellohey";
	unsigned char* d = reinterpret_cast<unsigned char*>(a);
	buff = reinterpret_cast<unsigned char*>(a);
	
	std::cout << "Init Job";
	//checkCudaErrors(cudaMallocManaged(&js, sizeof(JOB)));	//j = (JOB *)malloc(sizeof(JOB));
	//checkCudaErrors(cudaMallocManaged(&(js.data), sizeof(a)));
	js[0].data = buff;
	js[0].size = sizeof(a);
	for (int i = 0; i < 64; i++)
	{
		js[0].digest[i] = 0xff;
	}
	strcpy(js[0].fname, "1");

	JOB jsr[sizeof(js)/sizeof(JOB)];

	checkCudaErrors(cudaMalloc(&cuda_jobs, sizeof(js)));
	checkCudaErrors(cudaMalloc(&cuda_jobs_r, sizeof(jsr)));

	checkCudaErrors(cudaMemcpy(cuda_jobs, js, sizeof(js), cudaMemcpyHostToDevice));
	
	sha256_cuda <<< 1, sizeof(js)/sizeof(JOB) >>> (cuda_jobs, cuda_jobs_r);

	checkCudaErrors(cudaMemcpy(jsr, cuda_jobs_r, sizeof(jsr), cudaMemcpyDeviceToHost));
	
	std::cout << jsr[0].digest;
	//print_jobs(js, 1);

}

void run_hash (uint32_t * input, int bitlength, uint32_t * outputlocation) {


	uint32_t *inputCuda = 0;
	uint32_t *outputlocationCuda = 0;

	cudaMalloc(&inputCuda, sizeof(uint32_t)*20);
	cudaMalloc(&outputlocationCuda, sizeof(uint32_t)*8);

	cudaMemcpy(inputCuda, input, sizeof(uint32_t)*20, cudaMemcpyHostToDevice);

	if (bitlength == 640) sha256_cuda_hash_640 <<< 1, sizeof(uint32_t)*8 >>> (inputCuda, outputlocationCuda);

	if (bitlength == 256) sha256_cuda_hash_256 <<< 1, sizeof(outputlocation) >>> (inputCuda, outputlocationCuda);

	cudaMemcpy(outputlocation, outputlocationCuda, sizeof(outputlocation), cudaMemcpyDeviceToHost);


	//std::cout << &input;
	//std::cout << &outputlocation;
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



JOB * JOB_init(BYTE * data, long size, char * fname) {
	JOB * j;
	std::cout << "Init Job";
	checkCudaErrors(cudaMallocManaged(&j, sizeof(JOB)));	//j = (JOB *)malloc(sizeof(JOB));
	checkCudaErrors(cudaMallocManaged(&(j->data), size));
	j->data = data;
	j->size = size;
	for (int i = 0; i < 64; i++)
	{
		j->digest[i] = 0xff;
	}
	strcpy(j->fname, fname);
	return j;
}

int buildCUDAJobs(ssize_t hash) {

	std::cout << "JOB???";

	BYTE * buff;
	JOB ** jobs;

	//unsigned char myString [] = "This is my string";
	//BYTE * byte = &myString[0];
	checkCudaErrors(cudaMallocManaged(&jobs, 1 * sizeof(JOB *)));

	char a[9] = "hellohey";
	unsigned char* d = reinterpret_cast<unsigned char*>(a);
	buff = reinterpret_cast<unsigned char*>(a);
	
	jobs[0] = JOB_init(buff, sizeof(a), "1");
	
	//print_jobs(jobs, 1);


	pre_sha256();
	
	std::cout << "run JOB???";

	runJobs(jobs, 1);

	std::cout << "after run JOB???";


	cudaDeviceSynchronize();
	cudaDeviceReset();
	return 10;
}