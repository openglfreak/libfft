#pragma once
#ifndef __FFT_INTERNAL_PRIMEFACTORIZE_H__
#define __FFT_INTERNAL_PRIMEFACTORIZE_H__

#include <vector>
#include <utility>

namespace fft
{
	namespace internal
	{
		template<typename T>
		class PrimeFactorize
		{
		private:
			static void countFactor(T& i, T fact, T& out)
			{
				out = 0;
				while ((i % fact) == 0)
				{
					i /= fact;
					out++;
				}
			}

			template<T fact>
			static void countFactor_const(T& i, T& out)
			{
				out = 0;
				while ((i % fact) == 0)
				{
					i /= fact;
					out++;
				}
			}

			static void countFactor_const2(T& i, T& out)
			{
				out = 0;
				while ((i & 0x01) == 0)
				{
					i >>= 1;
					out++;
				}
			}
		public:
			typedef std::pair<T,unsigned char> factor_t;
			typedef std::vector<factor_t> list_t;
			typedef std::vector<T> expanded_list_t;

			static expanded_list_t expand(const list_t& l)
			{
				size_t s = 0;
				for (size_t i = l.size(); i--;)
					s += l[i].second;
				expanded_list_t ret = expanded_list_t();
				ret.reserve(s);
				for (size_t i = 0; i < l.size(); i++)
					ret.insert(ret.end(), l[i].second, l[i].first);
				return ret;
			}

			static list_t factorize(T i)
			{
				static const factor_t v = factor_t(1, 1);
				static const unsigned char steps[8] = { 4, 2, 4, 2, 4, 6, 2, 6 };

				switch (i)
				{
					case 0:
						return list_t();
					case 1:
						return list_t(&v, &v + 1);
					default:
					{
						list_t ret = list_t();
						T count = 0;

						countFactor_const2(i, count);
						if (count)
							ret.push_back(factor_t(2, count));
						if (i == 1)
							return ret;

						countFactor_const<3>(i, count);
						if (count)
							ret.push_back(factor_t(3, count));
						if (i == 1)
							return ret;

						countFactor_const<5>(i, count);
						if (count)
							ret.push_back(factor_t(5, count));
						if (i == 1)
							return ret;

						T fact = 7;
						do {
							for (const unsigned char* ptr = &steps[0]; ptr != &steps[8] && i > 1; fact += *ptr++)
							{
								countFactor(i, fact, count);
								if (count)
									ret.push_back(factor_t(fact, count));

								if (ret.empty() && ((fact * fact) > i))
								{
									ret.push_back(factor_t(i, 1));
									return ret;
								}
							}
						} while (i > 1);
						return ret;
					}
				}
			}
		};
	}
}

#endif // __FFT_INTERNAL_PRIMEFACTORIZE_H__
