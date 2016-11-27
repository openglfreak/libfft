#ifndef __FFT_INTERNAL_NUM_H__
#define __FFT_INTERNAL_NUM_H__

#include <complex>
#ifdef _FFT_USE_C99_COMPLEX
#include <complex.h>
#endif // _FFT_USE_C99_COMPLEX
#include "../internal/pi.h"

namespace fft
{
	namespace internal
	{
		template<typename T>
		struct Num
		{
			typedef T type;
			typedef T part_type;
			enum {
				holds_integer = 1,
				holds_real = 0,
				holds_complex = 0
			};
			enum {
				is_integer = 1,
				is_real = 0,
				is_complex = 0
			};
			enum {
				number_type = 0
			};
			typedef void holds_integer_t;
			typedef void is_integer_t;

			static T real(T v)
			{
				return v;
			}
			static T imag(T v)
			{
				return 0;
			}

			static T conj(T v)
			{
				return v;
			}
			template<typename _P>
			static T _sin(_P const& p)
			{
				return (T)sinf((float)Num<_P>::real(p));
			}
			template<typename _P>
			static T _cos(_P const& p)
			{
				return (T)cosf((float)Num<_P>::real(p));
			}

			template<typename _R>
			static T make(_R const& v)
			{
				return (T)Num<_R>::real(v);
			}
			template<typename _R, typename _I>
			static T make(_R const& v, _I const&)
			{
				return (T)Num<_R>::real(v);
			}
			template<typename _P>
			static T polar(_P const& p)
			{
				return polar<_P,_P>(p, 1);
			}
			template<typename _P, typename _M>
			static T polar(_P const& p, _M const& m = 1)
			{
				return (p == 0) ? m : (p == fft::internal::Pi::get<T>() ? -m : 0);
			}
		};
		template<>
		struct Num<float>
		{
			typedef float type;
			typedef float part_type;
			enum {
				holds_integer = 1,
				holds_real = 1,
				holds_complex = 0
			};
			enum {
				is_integer = 0,
				is_real = 1,
				is_complex = 0
			};
			enum {
				number_type = 1
			};
			typedef void holds_integer_t;
			typedef void holds_real_t;
			typedef void is_real_t;

			static float real(float v)
			{
				return v;
			}
			static float imag(float v)
			{
				return 0;
			}

			static float conj(float v)
			{
				return v;
			}
			template<typename _P>
			static float sin(_P const& p)
			{
				return sinf((float)Num<_P>::real(p));
			}
			template<typename _P>
			static float cos(_P const& p)
			{
				return cosf((float)Num<_P>::real(p));
			}

			template<typename _R>
			static float make(_R const& v)
			{
				return (float)Num<_R>::real(v);
			}
			template<typename _R, typename _I>
			static float make(_R const& v, _I const&)
			{
				return (float)Num<_R>::real(v);
			}
			template<typename _P>
			static float polar(_P const& p)
			{
				return polar<_P,_P>(p, 1);
			}
			template<typename _P, typename _M>
			static float polar(_P const& p, _M const& m = 1)
			{
				return cosf((float)Num<_P>::real(p)) * (float)Num<_M>::real(m);
			}
		};
		template<>
		struct Num<double>
		{
			typedef double type;
			typedef double part_type;
			enum {
				holds_integer = 1,
				holds_real = 1,
				holds_complex = 0
			};
			enum {
				is_integer = 0,
				is_real = 1,
				is_complex = 0
			};
			enum {
				number_type = 1
			};
			typedef void holds_integer_t;
			typedef void holds_real_t;
			typedef void is_real_t;

			static double real(double v)
			{
				return v;
			}
			static double imag(double v)
			{
				return 0;
			}

			static double conj(double v)
			{
				return v;
			}
			template<typename _P>
			static double sin(_P const& p)
			{
				return ::sin((double)Num<_P>::real(p));
			}
			template<typename _P>
			static double cos(_P const& p)
			{
				return ::cos((double)Num<_P>::real(p));
			}

			template<typename _R>
			static double make(_R const& v)
			{
				return (double)Num<_R>::real(v);
			}
			template<typename _R, typename _I>
			static double make(_R const& v, _I const&)
			{
				return (double)Num<_R>::real(v);
			}
			template<typename _P>
			static double polar(_P const& p)
			{
				return polar<_P,_P>(p, 1);
			}
			template<typename _P, typename _M>
			static double polar(_P const& p, _M const& m = 1)
			{
				return ::cos((double)Num<_P>::real(p)) * (double)Num<_M>::real(m);
			}
		};
		template<>
		struct Num<long double>
		{
			typedef long double type;
			typedef long double part_type;
			enum {
				holds_integer = 1,
				holds_real = 1,
				holds_complex = 0
			};
			enum {
				is_integer = 0,
				is_real = 1,
				is_complex = 0
			};
			enum {
				number_type = 1
			};
			typedef void holds_integer_t;
			typedef void holds_real_t;
			typedef void is_real_t;

			static long double real(long double v)
			{
				return v;
			}
			static long double imag(long double v)
			{
				return 0;
			}

			static long double conj(long double v)
			{
				return v;
			}
			template<typename _P>
			static long double sin(_P const& p)
			{
				return sinl((long double)Num<_P>::real(p));
			}
			template<typename _P>
			static long double cos(_P const& p)
			{
				return cosl((long double)Num<_P>::real(p));
			}

			template<typename _R>
			static long double make(_R const& v)
			{
				return (long double)Num<_R>::real(v);
			}
			template<typename _R, typename _I>
			static long double make(_R const& v, _I const&)
			{
				return (long double)Num<_R>::real(v);
			}
			template<typename _P>
			static long double polar(_P const& p)
			{
				return polar<_P,_P>(p, 1);
			}
			template<typename _P, typename _M>
			static long double polar(_P const& p, _M const& m = 1)
			{
				return cosl((long double)Num<_P>::real(p)) * (long double)Num<_M>::real(m);
			}
		};
		#ifdef _FFT_USE_C99_COMPLEX
		template<>
		struct Num<_Complex float>
		{
			typedef _Complex float type;
			typedef float part_type;
			enum {
				holds_integer = 1,
				holds_real = 1,
				holds_complex = 1
			};
			enum {
				is_integer = 0,
				is_real = 0,
				is_complex = 1
			};
			enum {
				number_type = 2
			};
			typedef void holds_integer_t;
			typedef void holds_real_t;
			typedef void holds_complex_t;
			typedef void is_complex_t;

			static float real(_Complex float v)
			{
				return crealf(v);
			}
			static float imag(_Complex float v)
			{
				return cimagf(v);
			}

			static _Complex float conj(_Complex float v)
			{
				return conjf(v);
			}
			template<typename _P>
			static _Complex float sin(_P const& p)
			{
				return (_Complex float)sinf((float)Num<_P>::real(p));
			}
			template<typename _P>
			static _Complex float cos(_P const& p)
			{
				return (_Complex float)cosf((float)Num<_P>::real(p));
			}

			template<typename _R>
			static _Complex float make(_R const& r)
			{
				return (float)Num<_R>::real(r) + (float)Num<_R>::imag(r) * I;
			}
			template<typename _R, typename _I>
			static _Complex float make(_R const& r, _I const& i)
			{
				return (float)Num<_R>::real(r) + (float)Num<_I>::real(i) * I;
			}
			template<typename _P>
			static _Complex float polar(_P const& p)
			{
				return polar<_P,_P>(p, 1);
			}
			template<typename _P, typename _M>
			static _Complex float polar(_P const& p, _M const& m = 1)
			{
				return (cosf((float)Num<_P>::real(p)) + sinf((float)Num<_P>::real(p)) * I) * (float)Num<_M>::real(m);
			}
		};
		template<>
		struct Num<_Complex double>
		{
			typedef _Complex double type;
			typedef double part_type;
			enum {
				holds_integer = 1,
				holds_real = 1,
				holds_complex = 1
			};
			enum {
				is_integer = 0,
				is_real = 0,
				is_complex = 1
			};
			enum {
				number_type = 2
			};
			typedef void holds_integer_t;
			typedef void holds_real_t;
			typedef void holds_complex_t;
			typedef void is_complex_t;

			static double real(_Complex double v)
			{
				return creal(v);
			}
			static double imag(_Complex double v)
			{
				return cimag(v);
			}

			static _Complex double conj(_Complex double v)
			{
				return ::conj(v);
			}
			template<typename _P>
			static _Complex double sin(_P const& p)
			{
				return (_Complex double)::sin((double)Num<_P>::real(p));
			}
			template<typename _P>
			static _Complex double cos(_P const& p)
			{
				return (_Complex double)::cos((double)Num<_P>::real(p));
			}

			template<typename _R>
			static _Complex double make(_R const& r)
			{
				return (double)Num<_R>::real(r) + (double)Num<_R>::imag(r) * I;
			}
			template<typename _R, typename _I>
			static _Complex double make(_R const& r, _I const& i)
			{
				return (double)Num<_R>::real(r) + (double)Num<_I>::real(i) * I;
			}
			template<typename _P>
			static _Complex double polar(_P const& p)
			{
				return polar<_P,_P>(p, 1);
			}
			template<typename _P, typename _M>
			static _Complex double polar(_P const& p, _M const& m = 1)
			{
				return (::cos((double)Num<_P>::real(p)) + ::sin((double)Num<_P>::real(p)) * I) * (double)Num<_M>::real(m);
			}
		};
		template<>
		struct Num<_Complex long double>
		{
			typedef _Complex long double type;
			typedef long double part_type;
			enum {
				holds_integer = 1,
				holds_real = 1,
				holds_complex = 1
			};
			enum {
				is_integer = 0,
				is_real = 0,
				is_complex = 1
			};
			enum {
				number_type = 2
			};
			typedef void holds_integer_t;
			typedef void holds_real_t;
			typedef void holds_complex_t;
			typedef void is_complex_t;

			static long double real(_Complex long double v)
			{
				return creall(v);
			}
			static long double imag(_Complex long double v)
			{
				return cimagl(v);
			}

			static _Complex long double conj(_Complex long double v)
			{
				return conjl(v);
			}
			template<typename _P>
			static _Complex long double sin(_P const& p)
			{
				return (_Complex long double)sinl((long double)Num<_P>::real(p));
			}
			template<typename _P>
			static _Complex long double cos(_P const& p)
			{
				return (_Complex long double)cosl((long double)Num<_P>::real(p));
			}

			template<typename _R>
			static _Complex long double make(_R const& r)
			{
				return (long double)Num<_R>::real(r) + (long double)Num<_R>::imag(r) * I;
			}
			template<typename _R, typename _I>
			static _Complex long double make(_R const& r, _I const& i)
			{
				return (long double)Num<_R>::real(r) + (long double)Num<_I>::real(i) * I;
			}
			template<typename _P>
			static _Complex long double polar(_P const& p)
			{
				return polar<_P,_P>(p, 1);
			}
			template<typename _P, typename _M>
			static _Complex long double polar(_P const& p, _M const& m)
			{
				return (cosl((long double)Num<_P>::real(p)) + sinl((long double)Num<_P>::real(p)) * I) * (long double)Num<_M>::real(m);
			}
		};
		#endif // _FFT_USE_C99_COMPLEX
		template<typename T>
		struct Num<std::complex<T> >
		{
			typedef std::complex<T> type;
			typedef T part_type;
			enum {
				holds_integer = 1,
				holds_real = 1,
				holds_complex = 1
			};
			enum {
				is_integer = 0,
				is_real = 0,
				is_complex = 1
			};
			enum {
				number_type = 3
			};
			typedef void holds_integer_t;
			typedef void holds_real_t;
			typedef void holds_complex_t;
			typedef void is_complex_t;

			static T real(std::complex<T> const& v)
			{
				return v.real();
			}
			static T imag(std::complex<T> const& v)
			{
				return v.imag();
			}

			static std::complex<T> conj(std::complex<T> const& v)
			{
				return std::conj(v);
			}
			template<typename _P>
			static std::complex<T> sin(_P const& p)
			{
				return std::complex<T>(Num<T>::sin((T)Num<_P>::real(p)));
			}
			template<typename _P>
			static std::complex<T> cos(_P const& p)
			{
				return std::complex<T>(Num<T>::cos((T)Num<_P>::real(p)));
			}

			template<typename _R>
			static std::complex<T> make(_R const& r)
			{
				return std::complex<T>((T)Num<_R>::real(r), (T)Num<_R>::imag(r));
			}
			template<typename _R, typename _I>
			static std::complex<T> make(_R const& r, _I const& i)
			{
				return std::complex<T>((T)Num<_R>::real(r), (T)Num<_I>::real(i));
			}
			template<typename _P>
			static std::complex<T> polar(_P const& p)
			{
				return polar<_P,_P>(p, 1);
			}
			template<typename _P, typename _M>
			static std::complex<T> polar(_P const& p, _M const& m)
			{
				return std::polar<T>((T)Num<_M>::real(m), (T)Num<_P>::real(p));
			}
		};
	}
}

#endif // __FFT_INTERNAL_NUM_H__
