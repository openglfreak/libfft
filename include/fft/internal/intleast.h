#ifndef __FFT_INTERNAL_INTLEAST_H__
#define __FFT_INTERNAL_INTLEAST_H__

#include <stdint.h>

namespace fft
{
	namespace internal
	{
		template<int>
		struct UIntLeast {};
		template<>
		struct UIntLeast<1>
		{
			typedef uint_least8_t type;
		};
		template<>
		struct UIntLeast<2>
		{
			typedef uint_least16_t type;
		};
		template<>
		struct UIntLeast<4>
		{
			typedef uint_least32_t type;
		};
		template<>
		struct UIntLeast<3> : public UIntLeast<4> {};
		template<>
		struct UIntLeast<8>
		{
			typedef uint_least64_t type;
		};
		template<>
		struct UIntLeast<5> : public UIntLeast<8> {};
		template<>
		struct UIntLeast<6> : public UIntLeast<8> {};
		template<>
		struct UIntLeast<7> : public UIntLeast<8> {};

		template<int>
		struct IntLeast {};
		template<>
		struct IntLeast<1>
		{
			typedef int_least8_t type;
		};
		template<>
		struct IntLeast<2>
		{
			typedef int_least16_t type;
		};
		template<>
		struct IntLeast<4>
		{
			typedef int_least32_t type;
		};
		template<>
		struct IntLeast<3> : public IntLeast<4> {};
		template<>
		struct IntLeast<8>
		{
			typedef int_least64_t type;
		};
		template<>
		struct IntLeast<5> : public IntLeast<8> {};
		template<>
		struct IntLeast<6> : public IntLeast<8> {};
		template<>
		struct IntLeast<7> : public IntLeast<8> {};
	}
}

#endif // __FFT_INTERNAL_INTLEAST_H__
