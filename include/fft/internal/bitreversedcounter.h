#ifndef __FFT_INTERNAL_BITREVERSEDCOUNTER_H__
#define __FFT_INTERNAL_BITREVERSEDCOUNTER_H__

#include <intrin.h>
#include "../internal/intleast.h"

namespace fft
{
	namespace internal
	{
		template<unsigned char _Bits>
		struct BitReversedCounter
		{
			enum { bits = _Bits, bytes = (_Bits+7)/8 };
			typedef typename UIntLeast<bytes>::type type;
			type value;

			BitReversedCounter() : value(0) {}
			BitReversedCounter(type v) : value(v & _make_mask()) {}

			type get()
			{
				return value;
			}
			operator type const&()
			{
				return value;
			}
			void set(type v)
			{
				value = v & _make_mask();
			}
			BitReversedCounter& operator=(type v)
			{
				value = v & _make_mask();
				return *this;
			}

			BitReversedCounter& operator++()
			{
				_inc();
				return *this;
			}
			// so I don't accidentally use this
			/*BitReversedCounter operator++(int)
			{
				BitReversedCounter temp = *this;
				this->operator++();
				return temp;
			}*/

			BitReversedCounter& operator--()
			{
				_dec();
				return *this;
			}
			// ^
			/*BitReversedCounter operator--(int)
			{
				BitReversedCounter temp = *this;
				this->operator--();
				return temp;
			}*/
		private:
			inline type _make_mask()
			#ifdef __GNUC__
			__attribute__ ((always_inline))
			#endif
			{
				// mask to remove upper bits
				return ((type)1 << _Bits) - 1;
			}
			#if defined x__GNUC__
			void _inc()
			{
				if (value == _make_mask())
					value = 0;
				else
				{
					// this is black magic, don't even try to understand it!
					// (it basically finds the most significant zero, shifts 3 to that position, and adds that to the original value)
					// example:
					//   10110000b
					//    ^ most significant zero bit (bit #6)
					//
					//   3<<6
					// = 11000000b
					//
					//   11000000b
					// + 10110000b
					// = 01110000b
					//
					// 01110000b is the new value
					//
					// IT WORKS, DON'T TOUCH IT!!!
					value += (type)3 << ((sizeof(unsigned long long) * CHAR_BIT - 1) - __builtin_clzll((unsigned long long)~value & _make_mask()));
					value &= _make_mask();
				}
			}
			void _dec()
			{
				if (value == 0)
					value = _make_mask();
				else
				{
					value -= (type)3 << ((sizeof(unsigned long long) * CHAR_BIT - 1) - __builtin_clzll(value));
					value &= _make_mask();
				}
			}
			#elif defined(x_MSC_VER)
			#pragma intrinsic(_BitScanReverse)
			#ifdef _WIN64
			#pragma intrinsic(_BitScanReverse64)
			#define __bsr64 _BitScanReverse64
			#else
			void __bsr64(unsigned long* o, __int64 i)
			{
				if (i > 0xFFFFFFFF)
				{
					_BitScanReverse(o, (unsigned long)(i >> 32));
					*o += 32;
				}
				else
					_BitScanReverse(o, (unsigned long)i);
			}
			#endif

			void _inc()
			{
				if (value == _make_mask())
					value = 0;
				else
				{
					unsigned long i;
					if (_Bits > 32)
						__bsr64(&i, (__int64)~value & _make_mask());
					else
						_BitScanReverse(&i, (unsigned long)~value & _make_mask());
					value += (type)3 << i;
					value &= _make_mask();
				}
			}
			void _dec()
			{
				if (value == 0)
					value = _make_mask();
				else
				{
					unsigned long i;
					if (_Bits > 32)
						__bsr64(&i, value);
					else
						_BitScanReverse(&i, value);
					value += (type)3 << i;
					value &= _make_mask();
				}
			}
			#else
			void _inc()
			{
				type mask = ((type)1 << (_Bits - 1));
				for (; mask && (value & mask); mask >>= 1) value ^= mask;
				value ^= mask;
			}
			void _dec()
			{
				type mask = ((type)1 << (_Bits - 1));
				for (; mask && !(value & mask); mask >>= 1) value ^= mask;
				value ^= mask;
			}
			#endif
		};
	}
}

#endif // __FFT_INTERNAL_BITREVERSEDCOUNTER_H__
