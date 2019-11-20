#include <cuda_runtime.h>
#include <stdio.h>
#include <ctype.h>
#include <cstdlib>
#include <fstream>
#include <vector>
#include <iostream>
#include <string>
#include <sstream>
#include <chrono>
#include <time.h>

#include "Image.h"
__device__ int medHelper(int* arr, int winSize, int height, int col, int row) {
	int g = 0;
	int rad = winSize / 2;
	int* mini = (int*)malloc(sizeof(int) * winSize * winSize);
	//printf("col: %d srow: %d h-rad: %d rad: %d \n", col, sRow, (height - rad), rad);
	if (col >= rad && col < (height - rad) && row >= rad && row < (height - rad)) {
		for (int k = col - rad; k <= col + rad; k++) {
			for (int j = row - rad; j <= row + rad; j++) {
				if (((j * height) + k) > (512 * 512)) {
					printf("error here \n");
					return;
				}
				///printf("col: %d srow: %d k: %d j: %d pos: %d \n", col, sRow, k, j, (j * height) + k);
				//printf("at %d arr: %d g: %d \n", (j * height) + k, arr[(j * height) + k], g);
				mini[g] = arr[(j * height) + k];
				g++;
			}
		}
		int key;
		int key2;
		for (int i = 1; i < (winSize * winSize); i++) {
			key = mini[i];
			key2 = i - 1;
			while (key2 >= 0 && mini[key2] > key) {
				mini[key2 + 1] = mini[key2];
				key2 = key2 - 1;
			}
			mini[key2 + 1] = key;
		}
		int median = mini[(winSize * winSize) / 2];//4
		//free(mini);
		return median;
	}
	return 0;
}
__global__ void medianFilter(int* arr, int* copy, int winSize, int height, int width, int config) {
	//same thing but 1d array
	int row = 0;
	int col = 0;
	int rad = winSize / 2;
		//rad: can be 1, 3, 5, or 7
		if (blockIdx.x == 0 && threadIdx.x == 0) {
			for (int k = 0; k < rad; k++) {
				for (int i = 0; i < height; i++) {
					//copy border to handle edge case
					copy[i + (k * height)] = arr[i + (k * height)];
					copy[k + (i * height)] = arr[k + (i * height)];
					copy[(height - 1 - k) + (i * height)] = arr[(height - 1 - k) + (i * height)];
					copy[i + ((height - 1 - k) * height)] = arr[i + ((height - 1 - k) * height)];
				}
			}
		}
		
		int start;
		if (config == 1) {
			start = ((blockIdx.x + 1) * 256) + (threadIdx.x * 4096) - 256;
		}
		else if (config == 2) {
			start = (blockIdx.x * 4096) + ((threadIdx.x + 1) * 256) - 256;
		}
		else if (config == 3) {
			start = (blockIdx.x * 2048) + ((threadIdx.x + 1) * 256) - 256; //8 threads, 128 blocks
		}
		else if (config == 4) {
			start = ((blockIdx.x + 1) * 256) + (threadIdx.x * 2048) - 256; //128 threads, 8 blocks
		}
		int end = start + 256;
		int sRow = start / 512;
		int sCol = start % 512;
		int sc = sCol + 256;
		//printf("start: %d rad: %d sCol: %d sRow: %d \n", start, rad, sCol, sRow);
		int* mini = (int*)malloc(sizeof(int) * winSize * winSize);
		//int mini[225];
		for (col = sCol; col < sc; col++) {
			int g = 0;
			//printf("col: %d srow: %d h-rad: %d rad: %d \n", col, sRow, (height - rad), rad);
			if (col >= rad && col < (height - rad) && sRow >= rad && sRow < (height - rad)) {
				for (int k = col - rad; k <= col + rad; k++) {
					for (int j = sRow - rad; j <= sRow + rad; j++) {
						if (((j * height) + k) > (512 * 512)) {
							printf("error here \n");
							return;
						}
						///printf("col: %d srow: %d k: %d j: %d pos: %d \n", col, sRow, k, j, (j * height) + k);
						//printf("at %d arr: %d g: %d \n", (j * height) + k, arr[(j * height) + k], g);
						mini[g] = arr[(j * height) + k];
						g++;
					}
				}
				int key;
				int key2;
				for (int i = 1; i < (winSize * winSize); i++) {
					key = mini[i];
					key2 = i - 1;
					while (key2 >= 0 && mini[key2] > key) {
						mini[key2 + 1] = mini[key2];
						key2 = key2 - 1;
					}
					mini[key2 + 1] = key;
				}
				int median = mini[(winSize * winSize) / 2];//4
				//printf("c: %d r: %d med: %d \n", col, sRow, median);
				
				if (((sRow * height) + col) >= (512 * 512)) {
					printf("prob \n");
					return;
				}
				//int median = medHelper(arr, winSize, height, col, row);
				copy[(sRow * height) + col] = median;
			}
		}
		free(mini);	
}
using namespace std;
int readImage(char fname[], Image& image);
int readImageHeader(char fname[], int& N, int& M, int& Q, bool& type);
int writeImage(char fname[], Image& image);
int main(int argc, char** argv) {
	if (argc != 4) {
		cout << "error: needs to be argc 3, format should be <int filter> <input file> <output file> \n";
		cout << "argc: " << argc << "\n";
		return 0;
	}
	int filterSize = atoi(argv[1]);
	if (filterSize != 3 && filterSize != 7 && filterSize != 11 && filterSize != 15) {
		cout << "error: filter size must be 3, 7, 11 or 15 \n";
		return 0;
	}
	string inputFile = argv[2];
	string outputFile = argv[3];
	int col = 0;
	int row = 0;
	int M, N, Q; // rows, cols, grayscale
	bool type;
	// read image header
	readImageHeader(argv[2], N, M, Q, type);
	cout << "header: N " << N << " M " << M << " Q " << Q << "\n";
	Image image(N, M, Q);
	readImage(argv[2], image);
	int configVer = 0;
	//make 1d array with cudamalloc
	//convert 2d array to 1d array
	int* bigOne = (int*) malloc(sizeof(int) * N * M);
	//size_t pitch;
	//try converting to 1d array first
	int z = 0;
	int* gold1 = (int*)malloc(sizeof(int) * N * M);
	int* gold2 = (int*)malloc(sizeof(int) * N * M);
	int* neutral = (int*)malloc(sizeof(int) * N * M);
	for (int i = 0; i < N; i++) {
		for (int j = 0; j < M; j++) {
			*(bigOne+z) = image.getPixelVal(j, i);
			*(gold1 + z) = image.getPixelVal(j, i);
			*(neutral + z) = image.getPixelVal(j, i);
			z++;
		}
	}
	//test
	int* bigDev;
	cudaMalloc(&bigDev, (sizeof(int) * M * N));
	cudaMemcpy(bigDev, bigOne, (sizeof(int) * M * N), cudaMemcpyHostToDevice);
	int* bigDev2;
	//cudaError_t cudStat;
	cudaMalloc(&bigDev2, sizeof(int) * M * N);
	//my gpu has maxwell architecture, meaning that its compute cabability is 5.2 and onwards, can run at most 1024 threads per block and approx 16 active blocks

	int rad = filterSize / 2;
	int* cp = (int*)malloc(sizeof(int) * N * M);
	int threadsPerBlock = 1024;
	//int blocksPerGrid = (512 * 512 + threadsPerBlock - 1) / threadsPerBlock;
	int* mini = (int*)malloc(sizeof(int) * filterSize * filterSize);
	//dim3 threadBlocks = (64, 64);
	//dim3 blockTot = (16, 16); //each block has 256 to deal with
	dim3 tb1 = (64, 64);
	dim3 bt1 = (16, 16);
	dim3 tb2 = (16, 16);
	dim3 bt2 = (64, 64);
	dim3 tb3 = (8, 8);
	dim3 bt3 = (128, 128);
	dim3 tb4 = (128, 128);
	dim3 bt4 = (8, 8);
	clock_t s1 = clock();
	//pre
	/*for (int k = 0; k < rad; k++) {
		for (int i = 0; i < M; i++) {
			//copy border to handle edge case
			//copy[i + (k * height)] = arr[i + (k * height)];
			cudaMemcpy((bigDev2+ i + (k * N)), (bigDev + i + (k * N)), sizeof(int), cudaMemcpyDeviceToDevice);
			//copy[k + (i * height)] = arr[k + (i * height)];
			cudaMemcpy((bigDev2 + k + (i * N)), (bigDev + k + (i * N)), sizeof(int), cudaMemcpyDeviceToDevice);
			//copy[(height - 1 - k) + (i * height)] = arr[(height - 1 - k) + (i * height)];
			cudaMemcpy((bigDev2 + (N - 1 - k) + (i * N)), (bigDev + (N - 1 - k) + (i * N)), sizeof(int), cudaMemcpyDeviceToDevice);
			//copy[i + ((height - 1 - k) * height)] = arr[i + ((height - 1 - k) * height)];
			cudaMemcpy((bigDev2 + i + (N - 1 - k) * N), (bigDev + i + (N - 1 - k) * N), sizeof(int), cudaMemcpyDeviceToDevice);
		}
	}*/
	//cudaError_t err2 = cudaDeviceGetLimit();
	//cudaError_t err2 = cudaDeviceSetLimit(cudaLimitMallocHeapSize, 264217728);
	//cout << "err2: " << cudaGetErrorString(err2) << "\n";
	cudaError_t err;
	medianFilter << <bt1, tb1 >> > (bigDev, bigDev2, filterSize, M, N, 1);
	err = cudaDeviceSynchronize();
	cout << "err " << cudaGetErrorString(err) << "\n";
	cudaMemcpy(bigOne, bigDev2, sizeof(int) * N * M, cudaMemcpyDeviceToHost);
	clock_t e1 = clock();
	clock_t s2 = clock();
	medianFilter << <bt2, tb2 >> > (bigDev, bigDev2, filterSize, M, N, 2);
	err = cudaDeviceSynchronize();
	cout << "err " << cudaGetErrorString(err) << "\n";
	cudaMemcpy(bigOne, bigDev2, sizeof(int) * N * M, cudaMemcpyDeviceToHost);
	clock_t e2 = clock();
	clock_t s3 = clock();
	medianFilter << <bt3, tb3 >> > (bigDev, bigDev2, filterSize, M, N, 3);
	err = cudaDeviceSynchronize();
	cout << "err " << cudaGetErrorString(err) << "\n";
	cudaMemcpy(bigOne, bigDev2, sizeof(int) * N * M, cudaMemcpyDeviceToHost);
	clock_t e3 = clock();
	clock_t s4 = clock();
	medianFilter << <bt4, tb4 >> > (bigDev, bigDev2, filterSize, M, N, 4);
	err = cudaDeviceSynchronize();
	cout << "err " << cudaGetErrorString(err) << "\n";
	cudaMemcpy(bigOne, bigDev2, sizeof(int) * N * M, cudaMemcpyDeviceToHost);
	clock_t e4 = clock();

	for (int k = 0; k < rad; k++) {
		for (int i = 0; i < M; i++) {
			//copy border
			gold2[i + (k * M)] = gold1[i + (k * M)];
			gold2[k + (i * M)] = gold1[k + (i * M)];
			gold2[(M - 1 - k) + (i * M)] = gold1[(N - 1 - k) + (i * M)];
			gold2[i + ((M - 1 - k) * M)] = gold1[i + ((N - 1 - k) * M)];
		}
	}
	for (row = rad; row < M - rad; row++) {
		for (col = rad; col < N - rad; col++) {
			int g = 0;
			for (int k = row - rad; k <= row + rad; k++) {
				for (int j = col - rad; j <= col + rad; j++) {
					mini[g] = gold1[k + (j * N)];
					g++;
				}
			}
			//sort mini via insertion sort
			int key;
			int key2;
			for (int i = 1; i < (filterSize * filterSize); i++) {
				key = mini[i];
				key2 = i - 1;
				while (key2 >= 0 && mini[key2] > key) {
					mini[key2 + 1] = mini[key2];
					key2 = key2 - 1;
				}
				mini[key2 + 1] = key;
			}
			//the resulting array length of the window will always be odd so can just take size/2 to be the median after sort
			int median = mini[(filterSize * filterSize) / 2];//4
			gold2[row + (col * M)] = median;
		}
	}
	free(mini);
	Image cp2(N, M, Q);
	for (int i = 0; i < M; i++) {
		for (int k = 0; k < N; k++) {
			cp2.setPixelVal(i, k, bigOne[i + (N * k)]);
		}
	}
	writeImage(argv[3], cp2);
	cout << "filter " << filterSize << " N " << N << " M " << M << " rad " << rad << "\n";
	float correctC = 0;
	float totC = N * M;
	for (int i = 0; i < M; i++) {
		for (int k = 0; k < N ; k++) {
			//cout << "at " << (i + (k * N)) << " b " << bigOne[i + (k * N)] << " g " << gold2[i + (k * N)] << " n " << neutral[i + (k * N)] << "\n";
			if (bigOne[i + (k * N)] == gold2[i + (k * N)]) {
				correctC++;
			}
			else {
				//cout << "at " << (i + (k * N)) << " b " << bigOne[i + (k * N)] << " g " << gold2[i + (k * N)] << " n " << neutral[i + (k * N)] << "\n";
			}
		}
	}
	float percent = correctC / totC;
	cout << "percent correct (version 4): " << percent * 100 << "\n";
	double time1 = (double)(e1 - s1);
	time1 = time1 / CLOCKS_PER_SEC;
	double time2 = (double)(e2 - s2);
	time2 = time2 / CLOCKS_PER_SEC;
	double time3 = (double)(e3 - s3);
	time3 = time3 / CLOCKS_PER_SEC;
	double time4 = (double)(e4 - s4);
	time4 = time4 / CLOCKS_PER_SEC;
	cout << "note: each thread processes medians for all 4 versions \n";
	cout << "time for kernel to run with 64 blocks and 16 threads per block: " << time1 << " seconds \n";
	cout << "time for kernel to run with 16 blocks and 64 threads per block: " << time2 << " seconds \n";
	cout << "time for kernel to run with 8 blocks and 128 threads per block: " << time3 << " seconds \n";
	cout << "time for kernel to run with 128 blocks and 8 threads per block: " << time4 << " seconds \n";
	return 0;
}

int readImage(char fname[], Image& image)
{
	int i, j;
	int N, M, Q;
	unsigned char* charImage;
	char header[100], * ptr;
	ifstream ifp;

	ifp.open(fname, ios::in | ios::binary);

	if (!ifp)
	{
		cout << "Can't read image: " << fname << endl;
		exit(1);
	}

	// read header
	ifp.getline(header, 100, '\n');
	if ((header[0] != 80) || (header[1] != 53))
	{
		cout << "Image " << fname << " is not PGM" << endl;
		exit(1);
	}

	ifp.getline(header, 100, '\n');
	while (header[0] == '#')
		ifp.getline(header, 100, '\n');

	M = strtol(header, &ptr, 0);
	N = atoi(ptr);

	ifp.getline(header, 100, '\n');
	Q = strtol(header, &ptr, 0);

	charImage = (unsigned char*) new unsigned char[M * N];

	ifp.read(reinterpret_cast<char*>(charImage), (M * N) * sizeof(unsigned char));

	if (ifp.fail())
	{
		cout << "Image " << fname << " has wrong size" << endl;
		exit(1);
	}

	ifp.close();

	int val;
	for (i = 0; i < N; i++)
		for (j = 0; j < M; j++)
		{
			val = (int)charImage[i * M + j];
			image.setPixelVal(i, j, val);
		}

	delete[] charImage;
	return (1);
}
int readImageHeader(char fname[], int& N, int& M, int& Q, bool& type)
{
	int i, j;
	unsigned char* charImage;
	char header[100], * ptr;
	ifstream ifp;

	ifp.open(fname, ios::in | ios::binary);

	if (!ifp)
	{
		cout << "Can't read image: " << fname << endl;
		exit(1);
	}

	// read header

	type = false; // PGM

	ifp.getline(header, 100, '\n');
	if ((header[0] == 80) && (header[1] == 53))
	{
		type = false;
	}
	else if ((header[0] == 80) && (header[1] == 54))
	{
		type = true;
	}
	else
	{
		cout << "Image " << fname << " is not PGM or PPM" << endl;
		exit(1);
	}

	ifp.getline(header, 100, '\n');
	while (header[0] == '#')
		ifp.getline(header, 100, '\n');

	M = strtol(header, &ptr, 0);
	N = atoi(ptr);

	ifp.getline(header, 100, '\n');

	Q = strtol(header, &ptr, 0);

	ifp.close();

	return(1);
}
int writeImage(char fname[], Image& image)
{
	int i, j;
	int N, M, Q;
	unsigned char* charImage;
	ofstream ofp;

	image.getImageInfo(N, M, Q);

	charImage = (unsigned char*) new unsigned char[M * N];
	int val;
	for (i = 0; i < N; i++)
	{
		for (j = 0; j < M; j++)
		{
			val = image.getPixelVal(i, j);
			charImage[i * M + j] = (unsigned char)val;
		}
	}

	ofp.open(fname, ios::out | ios::binary);
	if (!ofp)
	{
		cout << "Can't open file: " << fname << endl;
		exit(1);
	}

	ofp << "P5" << endl;
	ofp << M << " " << N << endl;
	ofp << Q << endl;
	ofp.write(reinterpret_cast<char*>(charImage), (M * N) * sizeof(unsigned char));
	if (ofp.fail())
	{
		cout << "Can't write image " << fname << endl;
		exit(0);
	}
	ofp.close();
	delete[] charImage;
	return(1);
}