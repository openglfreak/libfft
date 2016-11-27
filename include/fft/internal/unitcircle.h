#ifndef __FFT_INTERNAL_UNITCIRCLE_H__
#define __FFT_INTERNAL_UNITCIRCLE_H__

#include <complex>
#include <algorithm>
#include <math.h>
#include "../internal/pi.h"
#include "../internal/num.h"

namespace fft
{
	namespace internal
	{
		template<typename T, size_t _PartLength, size_t _Start, size_t _Length, int>
		class _UnitCircle
		{
		private:
			enum {
				copy_quarter = (_Length % 4 == 0) && (_Length / 4 < _PartLength),
				copy_half = (_Length % 2 == 0) && (_Length / 2 < _PartLength),
				copy_full = _Length < _PartLength
			};
		public:
			typedef T type;
			enum { part_length = _PartLength, start = _Start, length = _Length };
			T values[_PartLength];

			_UnitCircle()
			{
				size_t i;
				for (i = copy_quarter ? _Length / 4 : (copy_half ? _Length / 2 : _PartLength); i--;)
					values[i] = Num<T>::polar((i + _Start) * ((2 * Pi::get<typename Num<T>::part_type>()) / _Length));
				if (copy_quarter)
				{
					// If I'd just copy and multiply by i, real values would become zero, so I backwards-copy and flip the sign of the real part (f(x) = -conj(f(x-PI/2)))
					T* p = &values[_Length / 4];
					*p = Num<T>::make(0, 1);
					for (i = _Length / 4 + 1; i < (copy_half ? _Length / 2 : _PartLength); i++)
					{
						values[i] = -Num<T>::conj(*(--p));
						if (p == &values[0])
							p = &values[_Length / 4];
					}
				}
				if (copy_half)
					for (i = _Length / 2; i < (copy_full ? _Length : _PartLength); i++)
						values[i] = -values[i - _Length / 2];
				if (copy_full)
				{
					// could do it exponentially, but nah
					for (i = _Length; i < _PartLength - (_PartLength % _Length); i += _Length)
						std::copy(&values[0], &values[_Length], &values[i]);
					if (_PartLength % _Length != 0)
						std::copy(&values[0], &values[_PartLength % _Length], &values[_PartLength - (_PartLength % _Length)]);
				}
			}
		};
		template<typename T, size_t _PartLength, size_t _Start, size_t _Length>
		class _UnitCircle<T,_PartLength,_Start,_Length,0>
		{
		public:
			typedef T type;
			enum { part_length = _PartLength, start = _Start, length = _Length };
			T values[_PartLength];

			_UnitCircle()
			{
				values[0] = 1;
				signed char v = 1;
				for (size_t i = _Length / 2 - _Start % (_Length / 2); i < _PartLength; i += _Length / 2)
					values[i] = v = -v;
			}
		};

		template<typename T, size_t _PartLength, size_t _Start = 0, size_t _Length = _PartLength>
		class UnitCircle : public _UnitCircle<T,_PartLength,_Start,_Length,Num<T>::number_type> {};
	}
}

#endif // __FFT_INTERNAL_UNITCIRCLE_H__
