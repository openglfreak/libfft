#define _FFT_USE_C99_COMPLEX
//#define __BitReversedCounter_noopt
#include "../../include/fft/ctfft2/ctfft2.h"
#include <iostream>
#include <iomanip>
#include <limits>
#include <bitset>

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

int main(int argc, char* argv[])
{
	std::cout << std::fixed << std::setprecision(std::numeric_limits<fft_out_part_type>::digits10 + 2);

	std::chrono::time_point<std::chrono::steady_clock> start, end;
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
	std::cout << std::chrono::duration_cast<std::chrono::microseconds>(end - start).count() / 1000000.0 << std::endl;

	/*for (size_t i = fft_size; i--;)
		out[i] = std::conj(out[i]);
	ctfft->fft(&out[0], &out[0]);
	for (size_t i = fft_size; i--;)
		out[i] = std::conj(out[i]) / (fft_out_part_type)fft_size;

	for (int i = 0; i < fft_size; i++)
		std::cout << out[i] << std::endl;*/
	return 0;
}
