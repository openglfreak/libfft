#pragma once
#ifndef __FFT_INTERNAL_POW_H__
#define __FFT_INTERNAL_POW_H__

namespace fft
{
	namespace internal
	{
		template<long long B, unsigned long long E>
		struct Pow
		{
			enum { value = Pow<B,E-1>::value * B };
		};
		template<long long B>
		struct Pow<B,1>
		{
			enum { value = B };
		};
		template<long long B>
		struct Pow<B,0>
		{
			enum { value = 1 };
		};
	}
}

#endif // __FFT_INTERNAL_POW_H__
