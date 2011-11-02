//============================================================================
// Name        : cph
// Author      : Andrew Tausz <atausz@stanford.edu>
// Version     : 1.0
// Copyright   : Copyright © 2011 Andrew Tausz
// Description : A basic package for persistent homology in C++
//============================================================================

#ifndef RIGHT_OPEN_INTERVAL_HPP_
#define RIGHT_OPEN_INTERVAL_HPP_

namespace cph
{

	template<class T>
	class right_open_interval
	{
	private:
		T _start, _finish;
		bool _is_infinite;

	public:

		static const right_open_interval<T> make_interval(const T start, const T finish)
		{
			return right_open_interval(start, finish);
		}

		static const right_open_interval<T> make_interval(const T start)
		{
			return right_open_interval(start);
		}

		right_open_interval(const T start, const T finish) :
			_start(start), _finish(finish), _is_infinite(false)
		{
		}

		right_open_interval(const T start) :
			_start(start), _finish(0), _is_infinite(true)
		{
		}
		virtual ~right_open_interval()
		{
		}

		const T start() const
		{
			return _start;
		}

		const T finish() const
		{
			return _finish;
		}

		const bool is_infinite() const
		{
			return this->_is_infinite;
		}

		int compare(const right_open_interval<T> & other) const
		{
			if (_start < other._start)
			{
				return -1;
			}

			if (_start > other._start)
			{
				return 1;
			}

			if (this->_is_infinite)
			{
				if (other._is_infinite)
				{
					// both infinite
					return 0;
				}
				else
				{
					return 1;
				}
			}
			else
			{
				if (other._is_infinite)
				{
					return -1;
				}
				else
				{
					// both finite
					if (_finish < other._finish)
					{
						return -1;
					}

					if (_finish > other._finish)
					{
						return 1;
					}
				}
			}

			return 0;
		}

		bool operator <(const right_open_interval<T> & other) const
		{
			return (this->compare(other) < 0);
		}

		bool operator >(const right_open_interval<T> & other) const
		{
			return (this->compare(other) > 0);
		}

		bool operator <=(const right_open_interval<T> & other) const
		{
			return (this->compare(other) <= 0);
		}

		bool operator >=(const right_open_interval<T> & other) const
		{
			return (this->compare(other) >= 0);
		}

		bool operator ==(const right_open_interval<T> & other) const
		{
			return (this->compare(other) == 0);
		}

		friend std::ostream & operator <<(std::ostream & s, const right_open_interval<T> & value)
		{
			s << '[';

			s << value._start;
			s << ',';

			if (value._is_infinite)
			{
				s << "infinity";
			}
			else
			{
				s << value._finish;
			}

			s << ')';

			return s;
		}
	};

} /* namespace cph */
#endif /* RIGHT_OPEN_INTERVAL_HPP_ */
