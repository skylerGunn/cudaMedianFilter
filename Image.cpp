#include "Image.h"
#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <cstdlib>
using namespace std;
/**
 * Creates an Image 0x0
 */
//image class used to help with reading in the image and writing pgm file, taken from https://kamrablog.wordpress.com/2013/01/02/read-process-and-save-a-grayscale-pgm-image/
Image::Image() {
	N = 0;
	M = 0;
	Q = 0;
	pixelVal = 0;
}

/**
 * Creates an Image of numRows x numCols and creates the arrays for it
 *
 *
 * \param numRows number of rows
 * \param numCols number of cols
 * \param grayLevels grey level.
 * \author Florence Tupin
 * \date 10 oct. 2012
 * \sa Image()
 */
Image::Image(int numRows, int numCols, int grayLevels)
{
	N = numRows;
	M = numCols;
	Q = grayLevels;
	pixelVal = new int* [N];
	for (int i = 0; i < N; i++)
	{
		pixelVal[i] = new int[M];
		for (int j = 0; j < M; j++)
			pixelVal[i][j] = 0;
	}
}

/**
 * destroy image
 */
Image::~Image()
{
	N = 0;
	M = 0;
	Q = 0;
	for (int i = 0; i < N; i++)
		delete pixelVal[N];
	delete pixelVal;
	pixelVal = 0;
}

/**
 * copies oldImage into new Image object
 *
 * @param oldImage object Image to copie.
 *
 */
Image::Image(const Image& oldImage)
{
	N = oldImage.N;
	M = oldImage.M;
	Q = oldImage.Q;

	pixelVal = new int* [N];
	for (int i = 0; i < N; i++)
	{
		pixelVal[i] = new int[M];
		for (int j = 0; j < M; j++)
			pixelVal[i][j] = oldImage.pixelVal[i][j];
	}
}

/**
 * copies oldImage into whatever you = it to
 *
 * @param oldImage object Image
 */
void Image::operator=(const Image& oldImage)
{
	N = oldImage.N;
	M = oldImage.M;
	Q = oldImage.Q;

	pixelVal = new int* [N];
	for (int i = 0; i < N; i++)
	{
		pixelVal[i] = new int[M];
		for (int j = 0; j < M; j++)
			pixelVal[i][j] = oldImage.pixelVal[i][j];
	}
}

/**
 * sets the number of rows, columns and graylevels
 *
 */
void Image::setImageInfo(int numRows, int numCols, int maxVal)
{
	N = numRows;
	M = numCols;
	Q = maxVal;
}

/**
 * returns the number of rows, columns and gray levels
 */
void Image::getImageInfo(int& numRows, int& numCols, int& maxVal)
{
	numRows = N;
	numCols = M;
	maxVal = Q;
}

/**
 * returns the gray value of a specific pixel
 */
int Image::getPixelVal(int row, int col)
{
	return pixelVal[row][col];
}

/**
 * sets the gray value of a specific pixel
 */
void Image::setPixelVal(int row, int col, int value)
{
	pixelVal[row][col] = value;
}

/**
 * checks to see if a pixel is within the image, returns true or false
 */
bool Image::inBounds(int row, int col)
{
	if (row >= N || row < 0 || col >= M || col < 0)
		return false;
	//else
	return true;
}

/**
 * negates image
 */
void Image::negateImage(Image& oldImage)
{
	int rows, cols, gray;
	rows = N;
	cols = M;
	gray = Q;

	Image tempImage(N, M, Q);

	for (int i = 0; i < rows; i++)
	{
		for (int j = 0; j < cols; j++)
			tempImage.pixelVal[i][j] = -(pixelVal[i][j]) + 255;
	}

	oldImage = tempImage;
}

/**
 * based on users input and rotates it around the center of the image
 */
void Image::rotateImage(int theta, Image& oldImage)
{
	int r0, c0;
	int r1, c1;
	int rows, cols;
	rows = oldImage.N;
	cols = oldImage.M;
	Image tempImage(rows, cols, oldImage.Q);

	float rads = (theta * 3.14159265) / 180.0;

	r0 = rows / 2;
	c0 = cols / 2;

	for (int r = 0; r < rows; r++)
	{
		for (int c = 0; c < cols; c++)
		{
			r1 = (int)(r0 + ((r - r0) * cos(rads)) - ((c - c0) * sin(rads)));
			c1 = (int)(c0 + ((r - r0) * sin(rads)) + ((c - c0) * cos(rads)));

			if (inBounds(r1, c1))
			{
				tempImage.pixelVal[r1][c1] = oldImage.pixelVal[r][c];
			}
		}
	}

	for (int i = 0; i < rows; i++)
	{
		for (int j = 0; j < cols; j++)
		{
			if (tempImage.pixelVal[i][j] == 0)
				tempImage.pixelVal[i][j] = tempImage.pixelVal[i][j + 1];
		}
	}
	oldImage = tempImage;
}
/*int Image::writeImage(char fname[], Image& image)
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
int Image::readImageHeader(char fname[], int& N, int& M, int& Q, bool& type)
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
int Image::readImage(char fname[], Image& image)
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
}*/