#define _FFT_USE_C99_COMPLEX
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

/*#define fft_size 32//(32*1024*1024)
typedef fft::experimental::ctfft2::CTFFT2<fft_size,std::complex<double> > ctfft2;

ctfft2 ctfft;
double in[fft_size] = { 1, 2, 3, 2, 1 };
std::complex<double> out[fft_size];*/

int main(int argc, char* argv[])
{
	//std::cout << std::fixed << std::setprecision(std::numeric_limits<long double>::digits10 + 2);
	#if 0
	ctfft.fft(&in[0], &out[0]);
	for (int i = 0; i < fft_size; i++)
		//std::cout << (int)(*(unsigned char*)&ctfft.buffer[i]) << std::endl;
		std::cout << ctfft.buffer[i] << std::endl;
	#else
	//fft::internal::_UnitCircle<_Complex long double,32*1024*1024,0,2048*128,-1>* uc;
	fft::internal::UnitCircle<_Complex long double,32*1024*1024,0,2048*128>* uc;
	for (size_t i = 0x01; i--;)
		uc = new fft::internal::UnitCircle<_Complex long double,32*1024*1024,0,2048*128>();
	//new (&uc) fft::internal::_UnitCircle<_Complex long double,2048,0,128,-1>();
	#endif
	return 0;
}
