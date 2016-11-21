#ifndef __FFT_INTERNAL_PI_H__
#define __FFT_INTERNAL_PI_H__

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif // M_PI

#ifndef M_PIl
#define M_PIl 3.141592653589793238462643383279502884L
#endif // M_PIl

namespace fft
{
	namespace internal
	{
		struct Pi
		{
			template<typename T>
			static T get()
			{
				return T(3.1415926535897932384626433832795028841971693993751058209749445923078164062862089986280348253421170679L);
			}
		};
	}
}

#endif // __FFT_INTERNAL_PI_H__
