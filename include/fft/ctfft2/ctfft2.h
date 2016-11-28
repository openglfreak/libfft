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
			template<size_t N, typename CT>
			class CTFFT2
			{
			private:
				enum { size_log2 = fft::internal::Log<2,N>::value };
			public:
				enum { size = fft::internal::Pow<2,size_log2>::value };
				typedef CT complex_type;
			private:
				typedef fft::internal::Num<CT> _Num_CT;
				typedef fft::internal::UnitCircle<CT,size/2,0,size> _UnitCircle;

				CT unitCircle[(size / 2) + (size / 2 - 4)];
				CT buffer[size];
				CTFFT2(CTFFT2 const&) {}
			public:
				CTFFT2()
				{
					_UnitCircle::calc(&unitCircle[0]);
					for (size_t i = size / 2; i--;)
						unitCircle[i] = _Num_CT::conj(unitCircle[i]);
					for (CT* p = &unitCircle[0], *p2 = &unitCircle[size/2]; (p2 - p) > 4; p += 2, p2++)
						*p2 = *p;
				}

				template<typename _InIter, typename _OutIter>
				bool fft(_InIter in, _OutIter out)
				{
                    return fft(in, out, typename std::iterator_traits<_InIter>::iterator_category(), typename std::iterator_traits<_OutIter>::iterator_category());
				}
			private:
				template<typename _InIter, typename _OutIter>
				bool fft(_InIter in, _OutIter out, std::random_access_iterator_tag, std::random_access_iterator_tag)
				{
			        typedef typename std::iterator_traits<_InIter>::value_type IT;
			        typedef typename std::iterator_traits<_OutIter>::value_type OT;
					typedef fft::internal::Num<IT> _Num_IT;
					typedef fft::internal::Num<OT> _Num_OT;

					fft::internal::BitReversedCounter<size_log2-2> brc = fft::internal::BitReversedCounter<size_log2-2>();
					for (size_t i = size / 4; i--;)
					{
						IT i0 = in[i];
						IT i1 = in[i + size / 4];
						IT i2 = in[i + size / 2];
						IT i3 = in[i + size / 4 * 3];

						IT t0 = i0 + i2;
						i0 -= i2;

						IT t1 = i1 + i3;
						i1 -= i3;
						CT t0c = _Num_CT::make(_Num_IT::imag(i1), -_Num_IT::real(i1));

						CT* p = &buffer[(--brc).value << 2];
						*p = _Num_CT::make(t0 + t1);
						*++p = _Num_CT::make(i0) + t0c;
						*++p = _Num_CT::make(t0 - t1);
						*++p = _Num_CT::make(i0) - t0c;
					}

					if (size > 8)
					{
						CT* p = &unitCircle[(size /2) + (size / 2 - 4) - 4];
						for (size_t block = 4, shift = size_log2 - 3; block < size / 2; p -= (block <<= 1), shift--)
				            for (size_t j = size; j; j -= block)
				                for (size_t k = block; k--;)
				                {
				                    CT& i0 = buffer[--j - block];
					                CT i1 = p[k] * buffer[j];
				                    buffer[j] = i0 - i1;
				                    buffer[j - block] = i0 + i1;
				                }
					}
					for (size_t k = size / 2; k--;)
			        {
			            CT& i0 = buffer[k];
			            CT i1 = unitCircle[k] * buffer[k + size / 2];
			            out[k + size / 2] = _Num_OT::make(i0 - i1);
			            out[k] = _Num_OT::make(i0 + i1);
			        }
					return true;
				}
				template<typename _InIter, typename _OutIter>
				void fft(_InIter in, _OutIter out, std::random_access_iterator_tag, std::output_iterator_tag)
				{
					typedef typename std::iterator_traits<_OutIter>::value_type OT;
					typedef fft::internal::Num<OT> _Num_OT;

					if (!fft(in, &buffer[0], std::random_access_iterator_tag(), std::random_access_iterator_tag()))
						return false;
					for (CT* p = &buffer[0]; p != &buffer[size];)
					{
						*out = _Num_OT::make(*p);
						*++out = _Num_OT::make(*++p);
						*++out = _Num_OT::make(*++p);
						*++out = _Num_OT::make(*++p);
						++out; ++p;
					}
					return true;
				}
				template<typename _InIter, typename _OutIter>
				void fft(_InIter in, _OutIter out, std::forward_iterator_tag, std::output_iterator_tag)
				{
					static typename std::iterator_traits<_InIter>::value_type temp[size];
					typename std::iterator_traits<_InIter>::value_type* p = &temp;
					for (size_t i = size / 4; i--;)
					{
						*p = *in;
						*++p = *++in;
						*++p = *++in;
						*++p = *++in;
						++p; ++in;
					}
					return fft(&temp, out, std::random_access_iterator_tag(), typename std::iterator_traits<_OutIter>::iterator_category());
				}
			};
			template<typename CT>
			class CTFFT2<0,CT>
			{
			public:
				enum { size = 0 };
				typedef CT complex_type;
			public:
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
