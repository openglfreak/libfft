#pragma once
#ifndef __FFT_EXPERIMENTAL_CTFFT_CTFFT_H__
#define __FFT_EXPERIMENTAL_CTFFT_CTFFT_H__

#include <stddef.h>

#include <iterator>
#include <algorithm>
#include <vector>
#include <utility>

#include <fft/internal/num.h>
#include <fft/internal/primefactorize.h>

namespace fft
{
	namespace experimental
	{
		template<typename CT>
		class CTFFT
		{
		public:
			typedef CT value_type;
			typedef typename fft::internal::Num<CT>::part_type part_type;
		private:
			typedef fft::internal::Num<CT> _Num_CT;
			typedef fft::internal::PrimeFactorize<size_t> _PrimeFactorize;
		public:
			CTFFT() : size(), primeFactors(), twiddleFactors() {}
			CTFFT(size_t size) : size(size), primeFactors(_PrimeFactorize::factorize(size))
			{
				computeTwiddleFactors();
			}
			~CTFFT()
			{
				delete[] twiddleFactors;
			}
			CTFFT(const CTFFT<CT>& other) : size(other.size), primeFactors(other.primeFactors)
			{

			}
			#if __cplusplus >= 201103L
			CTFFT(CTFFT<CT>&& other) : CTFFT()
			{
				swap(*this, other);
			}
			#endif

			friend void swap(CTFFT& first, CTFFT& second)
			{
				using std::swap;
				swap(first.size, second.size);
				swap(first.primeFactors, second.primeFactors);
			}

			CTFFT& operator=(CTFFT other)
			{
				swap(*this, other);
				return *this;
			}

			template<typename InIter, typename OutIter>
			bool fft(InIter in, OutIter out)
			{
				typedef typename std::iterator_traits<InIter>::value_type IT;
		        typedef typename std::iterator_traits<OutIter>::value_type OT;
				typedef fft::internal::Num<IT> _Num_IT;
				typedef fft::internal::Num<OT> _Num_OT;
		        switch (size)
		        {
					case 0:
						break;
		            case 1:
		                *out++ = _Num_OT::make(*in++);
						break;
		            case 2:
		            {
		                IT i0 = *in++;
		                IT i1 = *in++;
		                *out++ = _Num_OT::make(i0 + i1);
		                *out++ = _Num_OT::make(i0 - i1);
						break;
		            }
		            case 4:
		            {
		                IT i0 = *in++;
		                IT i1 = *in++;
		                IT i2 = *in++;
		                IT i3 = *in++;

		                IT t0 = i0 + i2;
		                i0 -= i2;

		                IT t1 = i1 + i3;
		                i1 -= i3;
		                OT i1c = _Num_OT::make(_Num_IT::imag(i1), -_Num_IT::real(i1));

		                *out++ = _Num_OT::make(t0 + t1);
		                *out++ = _Num_OT::make(i0) + i1c;
		                *out++ = _Num_OT::make(t0 - t1);
		                *out++ = _Num_OT::make(i0) - i1c;
						break;
		            }
		            default:
		                return _fft(in, out, typename std::iterator_traits<InIter>::iterator_category(), typename std::iterator_traits<OutIter>::iterator_category());
		        }
				return true;
			}
		private:
			size_t size;
			_PrimeFactorize::list_t primeFactors;
			CT* twiddleFactors;

			template<typename InIter, typename OutIter>
			void _fft(InIter in, OutIter out, std::random_access_iterator_tag, std::random_access_iterator_tag)
			{

			}

			void computeTwiddleFactors()
			{
				twiddleFactors = new CT[size / primeFactors.first];
			}
		};
	}
}

#endif // __FFT_EXPERIMENTAL_CTFFT_CTFFT_H__
