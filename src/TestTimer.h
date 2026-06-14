#pragma once

class TestTimer
{
public:
	TestTimer();
	~TestTimer();

	struct timeval;	// MSVC defines this in winsock2.h!? 
	typedef struct timeval timeval;

	void tic();
	double toc(const char *printInfo);

private:
	bool m_running;
	int gettimeofday(struct timeval* tp, struct timezone* tzp);
};
