
#include <cstdlib>
#include <cstdio>
#include <cstring>
#include <fstream>
#include <cassert>
#include <exception>
#include <algorithm>

using namespace std;


float clamp(const float& inf, const float& sup, const float& value)
{
	return max(inf, min(sup, value));
}


class complex
{
public:
	float real;
	float imag;

	complex() : real(0), imag(0) {}

	complex(float _real, float _imag)
	{
		real = _real;
		imag = _imag;
	}

	float modul() const
	{
		return sqrtf(real * real + imag * imag);
	}
};


class Image
{
public:

	struct Rgb
	{
		float r, g, b;

		Rgb() : r(0), g(0), b(0) {}
		Rgb(float rr) : r(rr), g(rr), b(rr) {}
		Rgb(float rr, float gg, float bb) : r(rr), g(gg), b(bb) {}
	};

	unsigned int w, h;
	Rgb* pixels;

	Image() : w(0), h(0), pixels(nullptr) {}
	Image(const unsigned int& _w, const unsigned int& _h) : w(_w), h(_h), pixels(nullptr)
	{
		pixels = new Rgb[w * h];
		for (int i = 0; i < w * h; ++i) 
			pixels[i] = 0;
	}

	~Image() 
	{
		if (pixels != nullptr) delete[] pixels;
	}
};


Image readPPM(const char* filename)
{
	ifstream ifs(filename, ios::binary);
	
	Image img;
	try 
	{
		if (ifs.fail()) 
			throw("Can't open input file"); 
		
		string header;
		int w, h, b;
		ifs >> header;

		if (strcmp(header.c_str(), "P6") != 0) 
			throw("Can't read input file");
		
		ifs >> w >> h >> b;
		img.w = w; img.h = h;
		img.pixels = new Image::Rgb[w * h];
		ifs.ignore(256, '\n');
		unsigned char pix[3];
		
		for (int i = 0; i < w * h; ++i)
		{
			ifs.read(reinterpret_cast<char*>(pix), 3);
			img.pixels[i].r = pix[0] / 255.f;
			img.pixels[i].g = pix[1] / 255.f;
			img.pixels[i].b = pix[2] / 255.f;
		}
		ifs.close();
	}

	catch (const char* err)
	{
		fprintf(stderr, "%s\n", err);
		ifs.close();
	}

	return img;
}


void savePPM(string filename, int N, const complex* arr_r, const complex* arr_g, const complex* arr_b, const bool flag)
{
	ofstream ofs(filename, ios::out | ios::binary);
	ofs << "P6\n" << N << " " << N << "\n255\n";

	float max1 = -INFINITY;
	float max2 = -INFINITY;
	float max3 = -INFINITY;
	for (int i = 0; i < N * N; ++i)
	{
		if (arr_r[i].modul() > max1) max1 = arr_r[i].modul();
		if (arr_g[i].modul() > max2) max2 = arr_g[i].modul();
		if (arr_b[i].modul() > max3) max3 = arr_b[i].modul();
	}
	float scale1 = 255 / log(1 + max1);
	float scale2 = 255 / log(1 + max2);
	float scale3 = 255 / log(1 + max3);
	
	if (!flag)
	{
		for (int i = N / 2; i < N; ++i)
		{
			for (int j = N / 2; j < N; ++j)
			{
				unsigned char r = (unsigned char)(scale1 * log(1 + arr_r[i * N + j].modul()));
				unsigned char g = (unsigned char)(scale2 * log(1 + arr_g[i * N + j].modul()));
				unsigned char b = (unsigned char)(scale3 * log(1 + arr_b[i * N + j].modul()));
				ofs << r << g << b;
			}

			for (int j = 0; j < N / 2; ++j)
			{
				unsigned char r = (unsigned char)(scale1 * log(1 + arr_r[i * N + j].modul()));
				unsigned char g = (unsigned char)(scale2 * log(1 + arr_g[i * N + j].modul()));
				unsigned char b = (unsigned char)(scale3 * log(1 + arr_b[i * N + j].modul()));
				ofs << r << g << b;
			}
		}

		for (int i = 0; i < N / 2; ++i)
		{
			for (int j = N / 2; j < N; ++j)
			{
				unsigned char r = (unsigned char)(scale1 * log(1 + arr_r[i * N + j].modul()));
				unsigned char g = (unsigned char)(scale2 * log(1 + arr_g[i * N + j].modul()));
				unsigned char b = (unsigned char)(scale3 * log(1 + arr_b[i * N + j].modul()));
				ofs << r << g << b;
			}

			for (int j = 0; j < N / 2; ++j)
			{
				unsigned char r = (unsigned char)(scale1 * log(1 + arr_r[i * N + j].modul()));
				unsigned char g = (unsigned char)(scale2 * log(1 + arr_g[i * N + j].modul()));
				unsigned char b = (unsigned char)(scale3 * log(1 + arr_b[i * N + j].modul()));
				ofs << r << g << b;
			}
		}
	}

	else
	{
		for (int i = 0; i < N * N; ++i)
		{
			unsigned char r = (unsigned char)(255 * clamp(0, 1, arr_r[i].modul()));
			unsigned char g = (unsigned char)(255 * clamp(0, 1, arr_g[i].modul()));
			unsigned char b = (unsigned char)(255 * clamp(0, 1, arr_b[i].modul()));
			ofs << r << g << b;
		}
	}

	ofs.close();
}