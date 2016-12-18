// TODO: add pragma once to all headers
#include <stddef.h>

#define _FFT_USE_C99_COMPLEX
//#define __BitReversedCounter_noopt
#include <fft/abandoned/ctfft2/ctfft2.h>
#include <iostream>
#include <iomanip>
#include <limits>
#include <bitset>

#include <fft/experimental/ctfft/ctfft.h>

#ifdef _FFT_USE_C99_COMPLEX
std::ostream& operator <<(std::ostream& stream, _Complex float c)
{
	stream << crealf(c);
	if (cimagf(c) != 0)
		stream << " + " << cimagf(c) << 'j';
	return stream;
}
std::ostream& operator <<(std::ostream& stream, _Complex double c)
{
	stream << creal(c);
	if (cimag(c) != 0)
		stream << " + " << cimag(c) << 'j';
	return stream;
}
std::ostream& operator <<(std::ostream& stream, _Complex long double c)
{
	stream << creall(c);
	if (cimagl(c) != 0)
		stream << " + " << cimagl(c) << 'j';
	return stream;
}
#endif // _FFT_USE_C99_COMPLEX

#define fft_size (1*(32*1024*1024))
#define fft_in_type double
#define fft_int_part_type double
#define fft_int_type std::complex<fft_int_part_type>
#define fft_out_part_type double
#define fft_out_type std::complex<fft_out_part_type>
typedef fft::experimental::ctfft2::CTFFT2<fft_size,fft_int_type > ctfft2;

ctfft2* ctfft;
fft_in_type in_data[] = { 1, 2, 3, 2, 1 };
fft_in_type* in;
fft_out_type* out;

struct __dummy_ctor_class
{
	__dummy_ctor_class()
	{
		in = new fft_in_type[fft_size];
		out = new fft_out_type[fft_size];
		std::copy(&in_data[0], &in_data[sizeof(in_data) * sizeof(fft_in_type)], in);
	}
} __dummy_ctor_obj;

#include <chrono>

#define __optimize_barrier() __asm__ __volatile__ ("":::"memory")

fft::experimental::CTFFT<unsigned int> c;

void __fft2(const std::complex<double>* in, size_t n, size_t s, std::complex<double>* out)
{
	if (n == 1)
		*out = *in;
	else
	{
		size_t n2 = n / 2;
		__fft2(&in[0], n2, s * 2, &out[0]);
		__fft2(&in[s], n2, s * 2, &out[n2]);
		for (size_t k = 0; k < n2; k++)
		{
			std::complex<double> e = out[k];
			double d = k * 2 * M_PI / n;
			std::complex<double> o = out[k + n2] * std::conj(std::complex<double>(cos(d), sin(d)));
			out[k] = e + o;
			out[k + n2] = e - o;
		}
	}
}

void _dft(const std::complex<double>* in, size_t m, std::complex<double>* out)
{
	for (size_t k = 0; k < m; k++)
	{
		out[k] = 0;
		for (size_t n = 0; n < m; n++)
		{
			double d = (k * n) * (2 * M_PI) / m;
			out[k] += in[n] * conj(std::complex<double>(cos(d), sin(d)));
		}
	}
}

void __fft(const std::complex<double>* in, size_t n, size_t s, std::complex<double>* out, const size_t* factList)
{
	if (n == 1)
		*out = *in;
	else
	{
		size_t n2 = n / *factList;
		for (size_t k = 0; k < *factList; k++)
			__fft(&in[s * k], n2, s * *factList, &out[n2 * k], &factList[1]);
		for (size_t i = 1; i < *factList; i++)
			for (size_t k = 0; k < n2; k++)
			{
				double d = (i * k) * (2 * M_PI) / n;
				out[i * n2 + k] *= std::conj(std::complex<double>(cos(d), sin(d)));
			}
		std::complex<double>* temp = new std::complex<double>[*factList * 2];
		for (size_t k = 0; k < n2; k++)
		{
			for (size_t i = 0; i < *factList; i++)
				temp[i] = out[k + i * n2];

			_dft(temp, *factList, &temp[*factList]);

			for (size_t i = 0; i < *factList; i++)
				out[k + i * n2] = temp[*factList + i];
		}
		delete[] temp;
	}
}
void _fft(const std::complex<double>* in, size_t n, std::complex<double>* out)
{
	typedef fft::internal::PrimeFactorize<size_t> fact;
	__fft(in, n, 1, out, &fact::expand(fact::factorize(n))[0]);
}

int main(int argc, char* argv[])
{
	std::cout << std::fixed << std::setprecision(std::numeric_limits<fft_out_part_type>::digits10 + 2);

	std::complex<double>* in = new std::complex<double>[16];
	std::copy(&in_data[0], &in_data[5], in);
	std::complex<double>* out0 = new std::complex<double>[16];
	_fft(in, 15, out0);
	std::complex<double>* out1 = new std::complex<double>[16];
	_dft(in, 15, out1);
	for (int i = 0; i < 15; i++)
		if (std::abs(out0[i] - out1[i]) > 0.000000001)
			std::cout << i << ": " << out0[i] << '|' << out1[i] << std::endl;

	/*std::chrono::time_point<std::chrono::steady_clock> start, end;
	__optimize_barrier();
	ctfft = new ctfft2();
	__optimize_barrier();
	for (long double i = (long double)0x1FFFFFFF; i-- > 0;)
        __optimize_barrier();
	__optimize_barrier();
	start = std::chrono::steady_clock::now();
	__optimize_barrier();
	ctfft->fft(&in[0], &out[0]);
	__optimize_barrier();
	end = std::chrono::steady_clock::now();
	__optimize_barrier();
	delete ctfft;
	std::cout << std::chrono::duration_cast<std::chrono::microseconds>(end - start).count() / 1000000.0 << std::endl;*/

	/*for (size_t i = fft_size; i--;)
		out[i] = std::conj(out[i]);
	ctfft->fft(&out[0], &out[0]);
	for (size_t i = fft_size; i--;)
		out[i] = std::conj(out[i]) / (fft_out_part_type)fft_size;

	for (int i = 0; i < fft_size; i++)
		std::cout << out[i] << std::endl;*/
	return 0;
}
