#pragma once
#ifndef __FFT_INTERNAL_UNITCIRCLE_H__
#define __FFT_INTERNAL_UNITCIRCLE_H__

#include <stddef.h>

#include <complex>
#include <algorithm>
#include <iterator>
#include <math.h>

#include <fft/internal/pi.h>
#include <fft/internal/num.h>

namespace fft
{
	namespace internal
	{
		// THIS CLASS IS SHIT!!! (at least the implementation is fine, though not really commented)
		template<typename T, size_t _PartLength, size_t _Start, size_t _Length, int>
		class _UnitCircle
		{
		public:
			typedef T type;
			enum { part_length = _PartLength, start = _Start, length = _Length };
			T values[_PartLength];

			_UnitCircle()
			{
				_calc<T*>(&values[0], std::random_access_iterator_tag(), 0);
			}

			template<typename _OutIter>
			static void calc(_OutIter out)
			{
				_calc<_OutIter>(out, typename std::iterator_traits<_OutIter>::iterator_category(), (T)(typename std::iterator_traits<_OutIter>::value_type)T());
			}
		private:
			enum {
				copy_quarter = (_Length % 4 == 0) && (_Length / 4 < _PartLength),
				copy_half = (_Length % 2 == 0) && (_Length / 2 < _PartLength),
				copy_full = _Length < _PartLength
			};

			template<typename _OutIter>
			static void _calc(_OutIter out, std::random_access_iterator_tag, T)
			{
				size_t i;
				for (i = copy_quarter ? _Length / 4 : (copy_half ? _Length / 2 : _PartLength); i--;)
					out[i] = Num<T>::polar((i + _Start) * ((2 * Pi::get<typename Num<T>::part_type>()) / _Length));
				if (copy_quarter)
				{
					// If I'd just copy and multiply by i, real values would become zero, so I backwards-copy and flip the sign of the real part (f(x) = -conj(f(x-PI/2)))
					T* p = &out[_Length / 4];
					*p = Num<T>::make(0, 1);
					for (i = _Length / 4 + 1; i < (copy_half ? _Length / 2 : _PartLength); i++)
					{
						out[i] = -Num<T>::conj(*(--p));
						if (p == &out[0])
							p = &out[_Length / 4];
					}
				}
				if (copy_half)
					for (i = _Length / 2; i < (copy_full ? _Length : _PartLength); i++)
						out[i] = -out[i - _Length / 2];
				if (copy_full)
				{
					// could do it exponentially, but nah
					for (i = _Length; i < _PartLength - (_PartLength % _Length); i += _Length)
						std::copy(&out[0], &out[_Length], &out[i]);
					if (_PartLength % _Length != 0)
						std::copy(&out[0], &out[_PartLength % _Length], &out[_PartLength - (_PartLength % _Length)]);
				}
			}

			// TODO: implement _calc for non-random-access iterators (note: probably not gonna happen...)
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
				signed char v = 1 - ((_Start / (_Length / 2)) % 2) * 2;
				for (size_t i = _Start % (_Length / 2); i < _PartLength; i += _Length / 2, v = -v)
					values[i] = Num<T>::make(v);
			}

			template<typename _OutIter>
			static void calc(_OutIter out)
			{
				_calc<_OutIter>(out, typename std::iterator_traits<_OutIter>::iterator_category(), (T)(typename std::iterator_traits<_OutIter>::value_type)T());
			}
		private:
			template<typename _OutIter>
			static void _calc(_OutIter out, std::random_access_iterator_tag, T)
			{
				std::fill(out, &out[_PartLength], Num<T>::make(0));
				signed char v = 1 - ((_Start / (_Length / 2)) % 2) * 2;
				for (size_t i = _Start % (_Length / 2); i < _PartLength; i += _Length / 2, v = -v)
					out[i] = Num<T>::make(v);
			}
		};

		template<typename T, size_t _PartLength, size_t _Start, size_t _Length = _PartLength>
		class UnitCircle : public _UnitCircle<T,_PartLength,_Start,_Length,Num<T>::number_type> {};
	}
}

#endif // __FFT_INTERNAL_UNITCIRCLE_H__
