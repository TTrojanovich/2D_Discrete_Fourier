#define DISCRETE_FOURIER_1
#ifdef DISCRETE_FOURIER_1
#define _USE_MATH_DEFINES

#include <iostream>
#include <fstream>
#include <cmath>
#include <algorithm>
#include "header.h"

using namespace std;

const int N = 512; // size of a square image, power of 2


void DFT1D(const int N, const complex* in, complex* out)
{
	for (int k = 0; k < N; ++k)
	{
		out[k].real = 0; out[k].imag = 0;
		for (int n = 0; n < N; ++n)
		{
			out[k].real += in[n].real * (cos(2 * M_PI * n * k / (float)N)) - in[n].imag * (-sin(2 * M_PI * n * k / (float)N));
			out[k].imag += in[n].real * (-sin(2 * M_PI * n * k / (float)N)) + in[n].imag * (cos(2 * M_PI * n * k / (float)N));
		}
	}
}


void iDFT1D(const int N, const complex* in, complex* out)
{
	for (int k = 0; k < N; ++k)
	{
		out[k].real = 0; out[k].imag = 0;
		for (int n = 0; n < N; ++n)
		{
			out[k].real += in[n].real * (cos(2 * M_PI * n * k / (float)N)) - in[n].imag * (sin(2 * M_PI * n * k / (float)N));
			out[k].imag += in[n].real * (sin(2 * M_PI * n * k / (float)N)) + in[n].imag * (cos(2 * M_PI * n * k / (float)N));
		}
		out[k].real /= (float)N;
		out[k].imag /= (float)N;
	}
}


template <typename T>
void DFT2D(const int N, const complex* in, complex* out, T fourier_type)
{
	// process the rows
	for (int i = 0; i < N; i++)
	{
		fourier_type(N, in + i * (long long)N, out + i * (long long)N);
	}

	// process the columns
	complex *ca = new complex[N], *cb = new complex[N];
	for (int i = 0; i < N; i++)
	{
		for (int j = 0; j < N; j++)
		{
			// extract column
			ca[j].real = out[j * N + i].real;
			ca[j].imag = out[j * N + i].imag;
		}
		
		// perform 1D DFT on this column 
		fourier_type(N, ca, cb);

		for (int j = 0; j < N; j++)
		{
			// store result back in the array 
			out[j * N + i].real = cb[j].real;
			out[j * N + i].imag = cb[j].imag;
		}
	}

	delete[] ca; delete[] cb;
}


int main()
{
	Image imageData = readPPM("lena.ppm");

	cout << imageData.w << " " << imageData.h << endl;
	cout << "PROCEEDING..." << endl;

	complex* in_r = new complex[N * (long long)N];
	complex* in_g = new complex[N * (long long)N];
	complex* in_b = new complex[N * (long long)N];

	complex* out_r = new complex[N * (long long)N];
	complex* out_g = new complex[N * (long long)N];
	complex* out_b = new complex[N * (long long)N];

	for (int i = 0; i < N; ++i)
		for (int j = 0; j < N; ++j)
		{
			in_r[i * N + j] = complex(imageData.pixels[i * N + j].r, 0);
			in_g[i * N + j] = complex(imageData.pixels[i * N + j].g, 0);
			in_b[i * N + j] = complex(imageData.pixels[i * N + j].b, 0);
		}

	cout << "EXECUTION: Forward DFT" << endl;
	DFT2D(N, in_r, out_r, DFT1D);
	DFT2D(N, in_g, out_g, DFT1D);
	DFT2D(N, in_b, out_b, DFT1D);
	cout << "EXECUTION: Inverse DFT" << endl;
	DFT2D(N, out_r, in_r, iDFT1D);
	DFT2D(N, out_g, in_g, iDFT1D);
	DFT2D(N, out_b, in_b, iDFT1D);

	cout << "EXECUTION: Creation image for DFT " << endl;
	savePPM("fourier_optical_out.ppm", N, out_r, out_g, out_b, false);

	cout << "EXECUTION: Create image for inverse DFT" << endl;
	savePPM("fourier_out.ppm", N, in_r, in_g, in_b, true);

	delete[] in_r; delete[] in_g; delete[] in_b;
	delete[] out_r; delete[] out_g; delete[] out_b;
}


#endif
