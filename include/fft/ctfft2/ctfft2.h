#ifndef __FFT_CTFFT2_CTFFT2_H__
#define __FFT_CTFFT2_CTFFT2_H__

#include <complex>
#include <iterator>
#include <algorithm>
#include <math.h>
#ifdef _FFT_USE_C99_COMPLEX
#include <complex.h>
#endif // _FFT_USE_C99_COMPLEX
#include "../internal/pi.h"
#include "../internal/num.h"
#include "../internal/pow.h"
#include "../internal/log.h"
#include "../internal/unitcircle.h"
#include "../internal/bitreversedcounter.h"

namespace fft
{
	namespace experimental
	{
		namespace ctfft2
		{
			// // N = size, VT = Value Type, CT = Complex Type, RT = Real Type
			// // RT has to be implicitly convertible to CTs value_type, and explicitly or implicitly convertible to VT
			template<size_t N, typename CT>
			class CTFFT2
			{
			private:
				enum { nlog2 = fft::internal::Log<2,N>::value };
			public:
				enum { size = fft::internal::Pow<2,nlog2>::value };
				typedef CT complex_type;
			private:
				typedef fft::internal::Num<CT> _Num_CT;
				fft::internal::UnitCircle<CT,size/2,0,size> unitCircle;
			public:
				CT buffer[size];
			public:
				CTFFT2() {}

				template<typename _InIter, typename _OutIter>
				bool fft(_InIter in, _OutIter out)
				{
                    fft(in, out, typename std::iterator_traits<_InIter>::iterator_category());
                    return true;
				}
			private:
				template<typename _InIter, typename _OutIter>
				bool fft(_InIter in, _OutIter out, std::random_access_iterator_tag)
				{
			        typedef typename std::iterator_traits<_InIter>::value_type IT;
			        //typedef typename std::iterator_traits<_OutIter>::value_type OT;
					typedef fft::internal::Num<IT> _Num_IT;
					//typedef fft::internal::Num<OT> _Num_OT;
					fft::internal::BitReversedCounter<nlog2-2> brc;
					for (size_t i = size / 4; i--; --brc)
					{
						IT i0 = *in;
						IT i1 = *++in;
						IT i2 = *++in;
						IT i3 = *++in;
						++in;

						IT t0 = i0 + i2;
						i0 -= i2;

						IT t1 = i1 + i3;
						i1 -= i3;
						CT t0c = _Num_CT::make(_Num_IT::imag(i1), -_Num_IT::real(i1));

						CT* p = &buffer[brc.value << 2];
						*p = _Num_CT::make(t0 + t1);
						*++p = _Num_CT::make(i0) + t0c;
						*++p = _Num_CT::make(t0 - t1);
						*++p = _Num_CT::make(i0) - t0c;
					}
					std::copy(&buffer[0], &buffer[size], out);
					return false;
				}
				template<typename _InIter, typename _OutIter>
				void fft(_InIter in, _OutIter out, std::input_iterator_tag)
				{
					static typename std::iterator_traits<_InIter>::value_type temp[size];
					typename std::iterator_traits<_InIter>::value_type* p = &temp;
					for (size_t i = size; i--;)
						*p++ = *in++;
					return fft(&temp, out, std::random_access_iterator_tag());
				}
			};
			template<typename CT>
			class CTFFT2<0,CT>
			{
			public:
				enum { size = 0 };
				typedef CT complex_type;
			public:
				CTFFT2() {}

				template<typename _InIter, typename _OutIter>
				bool fft(_InIter in, _OutIter out)
				{
			        return true;
				}
			};
			template<typename CT>
			class CTFFT2<1,CT>
			{
			public:
				enum { size = 1 };
				typedef CT complex_type;
			public:
				CTFFT2() {}

				template<typename _InIter, typename _OutIter>
				bool fft(_InIter in, _OutIter out)
				{
					*out++ = fft::internal::Num<typename std::iterator_traits<_OutIter>::value_type>::make(*in++);
			        return true;
				}
			};
			template<typename CT>
			class CTFFT2<2,CT>
			{
			public:
				enum { size = 2 };
				typedef CT complex_type;
			public:
				CTFFT2() {}

				template<typename _InIter, typename _OutIter>
				bool fft(_InIter in, _OutIter out)
				{
			        typedef typename std::iterator_traits<_InIter>::value_type IT;
					typedef fft::internal::Num<typename std::iterator_traits<_OutIter>::value_type> _Num_OT;
					IT e = *in++;
					IT o = *in++;
					*out++ = _Num_OT::make(e + o);
					*out++ = _Num_OT::make(e - o);
			        return true;
				}
			};
			template<typename CT>
			class CTFFT2<3,CT> : public CTFFT2<2,CT> {};
			template<typename CT>
			class CTFFT2<4,CT>
			{
			public:
				enum { size = 4 };
				typedef CT complex_type;
			public:
				CTFFT2() {}

				template<typename _InIter, typename _OutIter>
				bool fft(_InIter in, _OutIter out)
				{
			        typedef typename std::iterator_traits<_InIter>::value_type IT;
			        typedef typename std::iterator_traits<_OutIter>::value_type OT;
					typedef fft::internal::Num<OT> _Num_OT;
					IT e1 = *in++;
					IT e2 = *in++;
					IT o1 = *in++;
					IT o2 = *in++;

					IT t1 = e1 + o1;
					e1 -= o1;

					IT t2 = e2 + o2;
					OT e2c = _Num_OT::make(e2 - o2);
					e2c = _Num_OT::make(e2c.imag(), -e2c.real());

					*out++ = t1 + t2;
					*out++ = _Num_OT::make(e1 + e2c);
					*out++ = t1 - t2;
					*out++ = _Num_OT::make(e1 - e2c);
					return true;
				}
			};
		}
	}
}

#endif // __FFT_CTFFT2_CTFFT2_H__
