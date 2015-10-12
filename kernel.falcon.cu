#include "utils.h"
#include <thrust/host_vector.h>
#include <thrust/device_vector.h>
#include <iostream>
#include <map>
#include <cuda_runtime.h>
#include "timer.h"
using namespace std;
using namespace thrust;
//#include "cuda_runtime.h"
#include "device_launch_parameters.h"

#include <stdio.h>

const int sharedMemSizeStatic = 2000; // make sure to not go over the Shared memory limit of 48 KBytes
__shared__ int localSpectrum[sharedMemSizeStatic];
cudaError_t addWithCuda(int *c, const int *a, const int *b, size_t size);


// Cuda kernels and Device functions
__global__ void addKernel(int *c, const int *a, const int *b)
{
    int i = threadIdx.x;
    c[i] = a[i] + b[i];
}


__device__ bool compareKmer(int* currentSpectrum, int currentSpectrumKmerAddress, int* comparedSpectrum, int comparedSpectrumKmerAddress, int fset)
{
	for (int i = 0; i < fset; i++)
	{
		if (currentSpectrum[currentSpectrumKmerAddress + i] !=  comparedSpectrum[comparedSpectrumKmerAddress + i])
			return false;
	}

	return true;
}




__device__ int compareSpectrumPair(int* currentSpectrum, int currentSpectrumSize, int currentSpectrumStart,
									int* comparedSpectrum, int comparedSpectrumSize, int comparedSpectrumStart,
									int* histogramData, int globalCurrentSpectrumIdx, int globalComparedSpectrumIdx,
									int numberOfBins, int fsetSize)
{
	//Serially go over the individual bins for fset comparisons
	//Improve the double for-loop later for faster execution
	int global = blockDim.x * blockIdx.x + threadIdx.x;
	//printf("Thread number: %i\t", global);
	int runningSumCurrent = 0;
	int runningSumCompared = 0;
	int score = 0;
	for (int i=0; i < (numberOfBins ); i++)
	{
		/*if current_spectrum_histogram[b] != 0 and compared_spectrum_histogram[b] != 0:

                    for ii in range(running_sum_i,running_sum_i+ current_spectrum_histogram[b]):
                        for jj in range(running_sum_j,running_sum_j+ compared_spectrum_histogram[b]):
                            # print(kmer_data[key_list[i]][ii],kmer_data[key_list[j]][jj])
                            if kmer_data[key_list[i]][ii] == kmer_data[key_list[j]][jj]:
                                current_comaparison_score += 1
                                # print("Kmer Match")*/

		
		int CurrentSpectrumBinVal = histogramData[globalCurrentSpectrumIdx*numberOfBins + i];
		int ComparedSpectrumBinVal = histogramData[globalComparedSpectrumIdx*numberOfBins + i];

		if (CurrentSpectrumBinVal != 0 && ComparedSpectrumBinVal != 0 ) //ensure both bins have non-zero number of kmers
		{
			//printf("%d \t %d \t %d \n",threadIdx.x, CurrentSpectrumBinVal, ComparedSpectrumBinVal);
			for (int ii = runningSumCurrent; ii< (runningSumCurrent + CurrentSpectrumBinVal); ii++)
			{
				for (int jj = runningSumCompared; jj < (runningSumCompared + ComparedSpectrumBinVal); jj++)
				{
					int currentSpectrumKmerAddress = currentSpectrumStart + ii;
					int comparedSpectrumKmerAddress = comparedSpectrumStart + jj;
					bool result = compareKmer(currentSpectrum,currentSpectrumKmerAddress, comparedSpectrum,comparedSpectrumKmerAddress,fsetSize); 
					if (result == true)
						score+=1;
				}
			}
		
		
		}

		runningSumCurrent += CurrentSpectrumBinVal;
		runningSumCompared += ComparedSpectrumBinVal;



	}
	
	
	return score;
}

__global__ void compareSpectra(int* inputData, int* histogramData, int* spectraSizeArray, int* spectraScanArray,
							   int numberOfSpectra, int totalElements, int TotalHistogramElements, int numberOfBins,
								const int fsetSize, int loopIteration,	
								int* d_test_array)
{
	//const int sharedMemSizeStatic = 1000;
	//Initialize localSpectrum to all zeros???
	
	

	int globalIdx = blockDim.x * blockIdx.x + threadIdx.x;
	int localIdx = threadIdx.x;
	int globalCurrentSpectrumIdx = blockDim.x * blockIdx.x + loopIteration;
	int currentSpectrumSize = spectraSizeArray[globalCurrentSpectrumIdx];
	int globalComparedSpectrumIdx =  blockDim.x * blockIdx.x + localIdx; //strided array access for corresponding threads after the current thread
	
	// Slow method of creating local copy of current spectrum for comparison
	if (localIdx == loopIteration)
	{
		
		//if (globalIdx == 2)
			//printf("Come here\n");


		for (int i = 0; i < sharedMemSizeStatic ; i++)
		{
			if (i < currentSpectrumSize && globalIdx == 0){
				localSpectrum[i] = inputData[i];  //Simple lookup in global data array
				//d_test_array[i] = inputData[i];
			}
			else if (i < currentSpectrumSize && globalIdx != 0){
				localSpectrum[i] = 	inputData[spectraScanArray[globalIdx] + i]; // need to know where to start reading the current spectrum
				//d_test_array[i] = inputData[spectraScanArray[globalIdx] + i];
			}

			else {
				localSpectrum[i] = 0; // Just set the extra elements in the shared array to zero for n
				//d_test_array[i] = 0;
			}
		}
	}

	/*if (localIdx == loopIteration)
	{
		for (int i = 0; i < sharedMemSizeStatic ; i++)
		{
			d_test_array[i + blockIdx.x * sharedMemSizeStatic] = localSpectrum[i];
		}
	}*/
	__syncthreads();

	// Compare current spectrum with the other (blockDim.x - i) spectra to the right. Each thread responsible for one comparison
	if (localIdx > loopIteration)
	{
		 
		int re = compareSpectrumPair(localSpectrum,currentSpectrumSize,0,inputData,spectraSizeArray[globalIdx],spectraScanArray[globalIdx],histogramData,
			globalCurrentSpectrumIdx, globalComparedSpectrumIdx ,numberOfBins,fsetSize);
		//printf("\n%d \t %d\n",globalIdx, re);

	}
	
	__syncthreads();

}




// Helper function for using CUDA to add vectors in parallel.
cudaError_t addWithCuda(int *c, const int *a, const int *b, size_t size)
{
    int *dev_a = 0;
    int *dev_b = 0;
    int *dev_c = 0;
    cudaError_t cudaStatus;

    // Choose which GPU to run on, change this on a multi-GPU system.
    cudaStatus = cudaSetDevice(0);
    if (cudaStatus != cudaSuccess) {
        fprintf(stderr, "cudaSetDevice failed!  Do you have a CUDA-capable GPU installed?");
        goto Error;
    }

    // Allocate GPU buffers for three vectors (two input, one output)    .
    cudaStatus = cudaMalloc((void**)&dev_c, size * sizeof(int));
    if (cudaStatus != cudaSuccess) {
        fprintf(stderr, "cudaMalloc failed!");
        goto Error;
    }

    cudaStatus = cudaMalloc((void**)&dev_a, size * sizeof(int));
    if (cudaStatus != cudaSuccess) {
        fprintf(stderr, "cudaMalloc failed!");
        goto Error;
    }

    cudaStatus = cudaMalloc((void**)&dev_b, size * sizeof(int));
    if (cudaStatus != cudaSuccess) {
        fprintf(stderr, "cudaMalloc failed!");
        goto Error;
    }

    // Copy input vectors from host memory to GPU buffers.
    cudaStatus = cudaMemcpy(dev_a, a, size * sizeof(int), cudaMemcpyHostToDevice);
    if (cudaStatus != cudaSuccess) {
        fprintf(stderr, "cudaMemcpy failed!");
        goto Error;
    }

    cudaStatus = cudaMemcpy(dev_b, b, size * sizeof(int), cudaMemcpyHostToDevice);
    if (cudaStatus != cudaSuccess) {
        fprintf(stderr, "cudaMemcpy failed!");
        goto Error;
    }

    // Launch a kernel on the GPU with one thread for each element.
    addKernel<<<1, size>>>(dev_c, dev_a, dev_b);

    // cudaDeviceSynchronize waits for the kernel to finish, and returns
    // any errors encountered during the launch.
    cudaStatus = cudaDeviceSynchronize();
    if (cudaStatus != cudaSuccess) {
        fprintf(stderr, "cudaDeviceSynchronize returned error code %d after launching addKernel!\n", cudaStatus);
        goto Error;
    }

    // Copy output vector from GPU buffer to host memory.
    cudaStatus = cudaMemcpy(c, dev_c, size * sizeof(int), cudaMemcpyDeviceToHost);
    if (cudaStatus != cudaSuccess) {
        fprintf(stderr, "cudaMemcpy failed!");
        goto Error;
    }

Error:
    cudaFree(dev_c);
    cudaFree(dev_a);
    cudaFree(dev_b);
    
    return cudaStatus;
}


//Helper function to print a device array in the console
void printDeviceArray(string msg, int* input, int size)
{
	int* aux = new int[size];
	checkCudaErrors(cudaMemcpy(aux, input, size*sizeof(int), cudaMemcpyDeviceToHost));

	cout<< msg << endl;
	for (int i = 0; i < size; i++)
		cout<< aux[i] << "\t" ;
	cout<<endl;
}

//Helper function to print a host array in the console
void printHostArray(string msg, int* input, int size)
{
	cout<< msg << endl;
	for (int i = 0; i < size; i++)
		cout<< input[i] << "\t" ;
	cout<<endl;
}

static inline void loadbar(unsigned int x, unsigned int n, unsigned int w = 50)
{
    
	if (x == (n-1))
	{
		float ratio = 1;
		int   c      =  ratio * w;
 
		cout << setw(3) << (int)(ratio*100) << "% [";
		for (int x=0; x<c; x++) cout << "=";
		for (int x=c; x<w; x++) cout << " ";
		cout << "]\r" << flush;

		return;
	}
	
	else
	{
		//if ( (x != n) && (x % (n/100+1) != 0) )
			//return;
 
		float ratio  =  x/(float)n;
	
		int   c      =  ratio * w;
 
		cout << setw(3) << (int)(ratio*100) << "% [";
		for (int x=0; x<c; x++) cout << "=";
		for (int x=c; x<w; x++) cout << " ";
		cout << "]\r" << flush;

		return;
	}
}

// Main host function to compute fset scores
void fset(int* inputData, int* inputHistogram, int* spectraSizes, int totalHistogramElements, int totalElements, int totalSpectra, int fsetSize, 
			int blocks, int threads, int* outputScores)
{
	
		// Caluclate partial sum of spectra Sizes
	int* spectraScanArray = new int[totalSpectra];
	int temp = 0;
	for (int i = 0; i < totalSpectra; i++)
	{
		spectraScanArray[i] = temp;
		temp+= spectraSizes[i];
	}


	int *d_inputData;
	int *d_inputHistogram;
	int *d_spectraSizeArray;
	int *d_spectraScanArray;
	int *d_numberOfBins;
	checkCudaErrors(cudaMalloc((void**)&d_inputData, totalElements*sizeof(int)));
	checkCudaErrors(cudaMemcpy(d_inputData, inputData, totalElements*sizeof(int), cudaMemcpyHostToDevice));
	checkCudaErrors(cudaMalloc((void**)&d_inputHistogram, totalHistogramElements*sizeof(int)));
	checkCudaErrors(cudaMemcpy(d_inputHistogram, inputHistogram, totalHistogramElements*sizeof(int), cudaMemcpyHostToDevice));
	checkCudaErrors(cudaMalloc((void**)&d_spectraSizeArray, totalSpectra*sizeof(int)));
	checkCudaErrors(cudaMemcpy(d_spectraSizeArray, spectraSizes, totalSpectra*sizeof(int), cudaMemcpyHostToDevice));
	checkCudaErrors(cudaMalloc((void**)&d_spectraScanArray, totalSpectra*sizeof(int)));
	checkCudaErrors(cudaMemcpy(d_spectraScanArray, spectraScanArray, totalSpectra*sizeof(int), cudaMemcpyHostToDevice));
	checkCudaErrors(cudaMalloc((void**)&d_numberOfBins, 1*sizeof(int)));
	checkCudaErrors(cudaMemset(d_numberOfBins, totalHistogramElements/totalSpectra, 1*sizeof(int)));

	
	//printHostArray("Size array:",spectraSizes,totalSpectra);
	//printHostArray("Scanned Size array:",spectraScanArray,totalSpectra);



	// Call to main kernel
	
	cout << "Number of blocks = "<<blocks <<"\t" << "Number of threads per block = " << threads <<endl;
	GpuTimer timer;
	timer.Start();
	for (int i = 0; i < threads; i++)
	{
		/*__global__ void compareSpectra(int* inputData, int* histogramData, int* spectraSizeArray, int* spectraScanArray,
							   int numberOfSpectra, int totalElements, int TotalHistogramElements,
								const int fset_size, int loopIteration,	
								int* d_test_array)
		*/

		//cout<< "Iteration number: " << i <<endl;
		int *d_test_array;
		checkCudaErrors(cudaMalloc((void**)&d_test_array, 4*sizeof(int)));
		compareSpectra<<<blocks, threads>>>(d_inputData,d_inputHistogram,d_spectraSizeArray, d_spectraScanArray, totalSpectra,totalElements,totalHistogramElements,totalHistogramElements/totalSpectra,fsetSize,i,d_test_array);
		cudaDeviceSynchronize(); 
		cudaError_t  code = cudaGetLastError();
		if (code != cudaSuccess) 
			printf ("Cuda error -- %s\n", cudaGetErrorString(code)); 
		
		//loadbar(i,threads); //Takes extra time to print, so maybe not include it after all


	}

	timer.Stop();
	printf("\nGPU code ran in: %f seconds\n", timer.Elapsed()/1000);

	return;

}

