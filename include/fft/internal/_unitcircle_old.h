#ifndef __FFT_INTERNAL_UNITCIRCLE_H__
#define __FFT_INTERNAL_UNITCIRCLE_H__

#include <complex>
#include <algorithm>
#include <math.h>
#ifdef _FFT_USE_C99_COMPLEX
#include <complex.h>
#endif // _FFT_USE_C99_COMPLEX
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
				copy_full = _Length < _PartLength && false // because copy_full is not implemented
			};
		public:
			T values[_PartLength];

			_UnitCircle()
			{
				size_t i;
				for (i = copy_quarter ? _Length / 4 : (copy_half ? _Length / 2 : _PartLength); i--;)
				{
					typename Num<T>::part_type x = (i + _Start) * ((2 * fft::internal::Pi::get<typename Num<T>::part_type>()) / _Length);
					values[i] = Num<T>::make(Num<typename Num<T>::part_type>::cos(x), Num<typename Num<T>::part_type>::sin(x));
				}
				if (copy_quarter)
					for (i = _Length / 4; i < (copy_half ? _Length / 2 : _PartLength); i++)
					{
						T& v = values[i - _Length / 4];
						values[i] = Num<T>::make(-Num<T>::imag(v), Num<T>::real(v));
					}
				if (copy_half)
					for (i = _Length / 2; i < (copy_full ? _Length : _PartLength); i++)
						values[i] = -values[i - _Length / 2];
				if (copy_full)
				{
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
			T values[_PartLength];

			_UnitCircle()
			{
				// values[x] = (T)cosf((x + _Start) * ((2 * M_PI) / _Length)) = [1,0,-1]
				// values[i*_Length/2-_Start] = (T)cosf(i * M_PI)
				values[0] = 1;
				signed char v = 1;
				for (size_t i = _Length / 2 - _Start % (_Length / 2); i < _PartLength; i += _Length / 2)
					values[i] = v = -v;
			}
		};
		template<typename T, size_t _PartLength, size_t _Start, size_t _Length>
		class _UnitCircle<T,_PartLength,_Start,_Length,1>
		{
		public:
			T values[_PartLength];

			_UnitCircle()
			{
				for (size_t i = _PartLength; i--;)
					values[i] = (T)cosl((i + _Start) * ((2 * M_PIl) / _Length));
			}
		};
		template<size_t _PartLength, size_t _Start, size_t _Length>
		class _UnitCircle<float,_PartLength,_Start,_Length,1>
		{
		public:
			float values[_PartLength];

			_UnitCircle()
			{
				for (size_t i = _PartLength; i--;)
					values[i] = cosf((i + _Start) * ((2 * M_PI) / _Length));
			}
		};
		template<size_t _PartLength, size_t _Start, size_t _Length>
		class _UnitCircle<double,_PartLength,_Start,_Length,1>
		{
		public:
			double values[_PartLength];

			_UnitCircle()
			{
				for (size_t i = _PartLength; i--;)
					values[i] = cos((i + _Start) * ((2 * M_PI) / _Length));
			}
		};
		template<size_t _PartLength, size_t _Start, size_t _Length>
		class _UnitCircle<long double,_PartLength,_Start,_Length,1>
		{
		public:
			long double values[_PartLength];

			_UnitCircle()
			{
				for (size_t i = _PartLength; i--;)
					values[i] = cosl((i + _Start) * ((2 * M_PIl) / _Length));
			}
		};
		#ifdef _FFT_USE_C99_COMPLEX
		template<typename T, size_t _PartLength, size_t _Start, size_t _Length>
		class _UnitCircle<T,_PartLength,_Start,_Length,2>
		{
		public:
			T values[_PartLength];

			_UnitCircle()
			{
				for (size_t i = _PartLength; i--;)
				{
					long double x = (i + _Start) * ((2 * M_PIl) / _Length);
					values[i] = (T)(cosl(x) + sinl(x) * I);
				}
			}
		};
		template<size_t _PartLength, size_t _Start, size_t _Length>
		class _UnitCircle<_Complex float,_PartLength,_Start,_Length,2>
		{
		public:
			_Complex float values[_PartLength];

			_UnitCircle()
			{
				for (size_t i = _PartLength; i--;)
				{
					float x = (i + _Start) * ((2 * M_PI) / _Length);
					values[i] = cosf(x) + sinf(x) * I;
				}
			}
		};
		template<size_t _PartLength, size_t _Start, size_t _Length>
		class _UnitCircle<_Complex double,_PartLength,_Start,_Length,2>
		{
		public:
			_Complex double values[_PartLength];

			_UnitCircle()
			{
				for (size_t i = _PartLength; i--;)
				{
					double x = (i + _Start) * ((2 * M_PI) / _Length);
					values[i] = cos(x) + sin(x) * I;
				}
			}
		};
		template<size_t _PartLength, size_t _Start, size_t _Length>
		class _UnitCircle<_Complex long double,_PartLength,_Start,_Length,2>
		{
		public:
			_Complex long double values[_PartLength];

			_UnitCircle()
			{
				for (size_t i = _PartLength; i--;)
				{
					long double x = (i + _Start) * ((2 * M_PIl) / _Length);
					values[i] = cosl(x) + sinl(x) * I;
				}
			}
		};
		#endif // _FFT_USE_C99_COMPLEX
		template<typename T, size_t _PartLength, size_t _Start, size_t _Length>
		class _UnitCircle<T,_PartLength,_Start,_Length,3>
		{
		public:
			T values[_PartLength];

			_UnitCircle()
			{
				for (size_t i = _PartLength; i--;)
				{
					long double x = (i + _Start) * ((2 * M_PIl) / _Length);
					values[i] = T(cosl(x), sinl(x));
				}
			}
		};
		template<size_t _PartLength, size_t _Start, size_t _Length>
		class _UnitCircle<std::complex<float>,_PartLength,_Start,_Length,3>
		{
		public:
			std::complex<float> values[_PartLength];

			_UnitCircle()
			{
				for (size_t i = _PartLength; i--;)
				{
					float x = (i + _Start) * ((2 * M_PI) / _Length);
					values[i] = std::complex<float>(cosf(x), sinf(x));
				}
			}
		};
		template<size_t _PartLength, size_t _Start, size_t _Length>
		class _UnitCircle<std::complex<double>,_PartLength,_Start,_Length,3>
		{
		public:
			std::complex<double> values[_PartLength];

			_UnitCircle()
			{
				for (size_t i = _PartLength; i--;)
				{
					double x = (i + _Start) * ((2 * M_PI) / _Length);
					values[i] = std::complex<double>(cos(x), sin(x));
				}
			}
		};
		template<size_t _PartLength, size_t _Start, size_t _Length>
		class _UnitCircle<std::complex<long double>,_PartLength,_Start,_Length,3>
		{
		public:
			std::complex<long double> values[_PartLength];

			_UnitCircle()
			{
				for (size_t i = _PartLength; i--;)
				{
					long double x = (i + _Start) * ((2 * M_PIl) / _Length);
					values[i] = std::complex<long double>(cosl(x), sinl(x));
				}
			}
		};

		template<typename T, size_t _PartLength, size_t _Start = 0, size_t _Length = _PartLength>
		class UnitCircle : public _UnitCircle<T,_PartLength,_Start,_Length,Num<T>::number_type> {};
	}
}

#endif // __FFT_INTERNAL_UNITCIRCLE_H__
