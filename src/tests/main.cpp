#define _FFT_USE_C99_COMPLEX
#define __BitReversedCounter_noopt
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

#define fft_size 32//(32*1024*1024)
#define fft_in_type long double
#define fft_int_part_type long double
#define fft_int_type std::complex<fft_int_part_type>
#define fft_out_part_type long double
#define fft_out_type std::complex<fft_out_part_type>
typedef fft::experimental::ctfft2::CTFFT2<fft_size,fft_int_type > ctfft2;

ctfft2 ctfft;
fft_in_type in[fft_size] = { 1, 2, 3, 2, 1 };
fft_out_type out[fft_size];

int main(int argc, char* argv[])
{
	std::cout << std::fixed << std::setprecision(std::numeric_limits<fft_out_part_type>::digits10 + 2);

	ctfft.fft(&in[0], &out[0]);

	for (size_t i = fft_size; i--;)
		out[i] = std::conj(out[i]);
	ctfft.fft(&out[0], &out[0]);
	for (size_t i = fft_size; i--;)
		out[i] = std::conj(out[i]) / (fft_out_part_type)fft_size;

	for (int i = 0; i < fft_size; i++)
		std::cout << out[i] << std::endl;
	return 0;
}
