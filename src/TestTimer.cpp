#define WIN32_LEAN_AND_MEAN
#include <iostream>
#include <Windows.h>
#include <stdint.h> // portable: uint64_t   MSVC: __int64
#include "TestTimer.h"

	TestTimer::TestTimer() = default;
	TestTimer::~TestTimer() = default;

	struct TestTimer::timeval {
		long tv_sec;
		long tv_usec;
	};

	struct TestTimer::timeval m_time_start
	{
		0, 0
	};

	struct TestTimer::timeval m_time_end
	{
		0, 0
	};

	int TestTimer::gettimeofday(struct timeval* tp, struct timezone* tzp)
	{
		// Note: some broken versions only have 8 trailing zero's;
		// the correct epoch has 9 trailing zero's
		// This magic number is the number of 100 nanosecond intervals
		// since January 1, 1601 (UTC) until 00:00:00 January 1, 1970
		static const uint64_t EPOCH = ((uint64_t)116444736000000000ULL);

		SYSTEMTIME  system_time;
		FILETIME	file_time;
		uint64_t	time;

		GetSystemTime(&system_time);
		SystemTimeToFileTime(&system_time, &file_time);
		time = ((uint64_t)file_time.dwLowDateTime);
		time += ((uint64_t)file_time.dwHighDateTime) << 32;

		tp->tv_sec = (long)((time - EPOCH) / 10000000L);
		tp->tv_usec = (long)(system_time.wMilliseconds * 1000);
		return 0;
	}

	void TestTimer::tic()
	{
		m_running = true;
		gettimeofday(&m_time_start, NULL);
	}

	double TestTimer::toc(const char *printInfo = nullptr)
	{
		double res = 0.0;
		if (m_running)
		{
			gettimeofday(&m_time_end, NULL);
			m_running = false;
			res = static_cast<double>(m_time_end.tv_sec - m_time_start.tv_sec)
				+ static_cast<double>(m_time_end.tv_usec - m_time_start.tv_usec) / 1000000.0;
			if (nullptr != printInfo)
			{
				std::cout << "Time elapsed - " << printInfo << ": ";
				if (res > 1.0)
					std::cout << res << " s\n";
				else if (res > 0.001)
					std::cout << res * 1000.0 << " ms\n";
				else
					std::cout << res * 1000000.0 << " us\n";
			}
		}
		return res;
	}
