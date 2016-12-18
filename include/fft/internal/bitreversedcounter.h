#pragma once
#ifndef __FFT_INTERNAL_BITREVERSEDCOUNTER_H__
#define __FFT_INTERNAL_BITREVERSEDCOUNTER_H__

#include <stdint.h>
#ifdef _MSC_VER
#include <intrin.h>
#endif
#include <limits.h>

#include <fft/internal/intleast.h>

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
			#ifdef __GNUC__
				__attribute__((always_inline))
			#elif defined(_MSC_VER)
				__forceinline
			#else
				inline
			#endif
			inline type _make_mask()
			{
				// mask to remove upper bits
				return ((type)1 << _Bits) - 1;
			}

			#ifdef __GNUC__
				static unsigned char _bsr64(int64_t i) __attribute__ ((always_inline))
				{
					return 63 - (sizeof(unsigned int) == 8 ? __builtin_clz(i) : (sizeof(unsigned long) == 8 ? __builtin_clzl(i) : __builtin_clzll(i)));
				}

				static unsigned char _bsr32(int32_t i) __attribute__ ((always_inline))
				{
					return 31 - (sizeof(unsigned int) == 4 ? __builtin_clz(i) : (sizeof(unsigned long) == 4 ? __builtin_clzl(i) : __builtin_clzll(i)));
				}
			#elif defined(_MSC_VER)
				#pragma intrinsic(_BitScanReverse)

				#ifdef _WIN64
					#pragma intrinsic(_BitScanReverse64)

					static __forceinline unsigned char _bsr64(int64_t i)
					{
						unsigned long ret;
						_BitScanReverse64(&ret, (__int64)i);
						return ret;
					}
				#else
					static __forceinline unsigned char _bsr64(int64_t i)
					{
						unsigned long ret;
						if (i > 0xFFFFFFFF)
						{
							_BitScanReverse(&ret, (unsigned long)(i >> 32));
							ret += 32;
						}
						else
							_BitScanReverse(&ret, (unsigned long)i);
						return ret;
					}
				#endif

				static __forceinline unsigned char _bsr32(int32_t i)
				{
					unsigned long ret;
					_BitScanReverse(&ret, (unsigned long)i);
					return ret;
				}
			#else
				#define __BitReversedCounter_noopt
			#endif

			#ifdef __BitReversedCounter_noopt
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
			#else
			void _inc()
			{
				if (value == _make_mask())
					value = 0;
				else
				{
					unsigned char i;
					if (_Bits > 32)
						i = _bsr64((uint64_t)~value & _make_mask());
					else
						i = _bsr32((uint32_t)~value & _make_mask());
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
						i = _bsr64(value);
					else
						i = _bsr32(value);
					value += (type)3 << i;
					value &= _make_mask();
				}
			}
			#endif
		};
	}
}

#endif // __FFT_INTERNAL_BITREVERSEDCOUNTER_H__
