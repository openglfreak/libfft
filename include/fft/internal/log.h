#ifndef __FFT_INTERNAL_LOG2_H__
#define __FFT_INTERNAL_LOG2_H__

namespace fft
{
	namespace internal
	{
		template<long long B, unsigned long long N>
		struct Log
		{
			enum { value = Log<B,N/B>::value + 1 };
		};
		template<long long B>
		struct Log<B,1>
		{
			enum { value = 0 };
		};
		template<long long B>
		struct Log<B,0> {};
	}
}

#endif // __FFT_INTERNAL_LOG2_H__
