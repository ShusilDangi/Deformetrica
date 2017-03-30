/***************************************************************************************
*                                                                                      *
*                                     Deformetrica                                     *
*                                                                                      *
*    Copyright Inria and the University of Utah.  All rights reserved. This file is    *
*    distributed under the terms of the Inria Non-Commercial License Agreement.        *
*                                                                                      *
*                                                                                      *
****************************************************************************************/

#ifndef _SimpleTimer_h
#define _SimpleTimer_h

#include <ctime>

/**
 *	\brief 		A simple timer.
 *
 *	\copyright		Inria and the University of Utah
 *	\version		Deformetrica 2.0
 *
 *	\details	The SimpleTimer class handles simple timings.
 */
class SimpleTimer
{

public:

	////////////////////////////////////////////////////////////////////////////////////////////////////
	// Constructor(s) / Destructor :
	////////////////////////////////////////////////////////////////////////////////////////////////////

	SimpleTimer();

	~SimpleTimer() {}



	////////////////////////////////////////////////////////////////////////////////////////////////////
	// Encapsulation method(s) :
	////////////////////////////////////////////////////////////////////////////////////////////////////

	/// Returns the number of elapsed time in seconds.
	inline double GetElapsedTimeInSecondsOnly() { this->Stop(); return m_ElapsedTimeInSecondsOnly; }

	/// Returns the number of elapsed CPU time in seconds.
	inline double GetElapsedCPUTimeInSecondsOnly() { this->Stop(); return m_ElapsedCPUTimeInSecondsOnly; }

	/// Returns the number of elapsed hours.
	inline unsigned int GetElapsedHours() { this->Stop(); return m_ElapsedHours; }
	/// Returns the number of elapsed minutes.
	inline unsigned int GetElapsedMinutes() { this->Stop(); return m_ElapsedMinutes; }
	/// Returns the number of elapsed seconds.
	inline unsigned int GetElapsedSeconds() { this->Stop(); return m_ElapsedSeconds; }


	////////////////////////////////////////////////////////////////////////////////////////////////////
	// Other method(s) :
	////////////////////////////////////////////////////////////////////////////////////////////////////

	/// Starts the timer.
	void Start();

	/// Stops the timer.
	void Stop();



private:

	////////////////////////////////////////////////////////////////////////////////////////////////////
	// Attribute(s)
	////////////////////////////////////////////////////////////////////////////////////////////////////

	clock_t m_StartClock;
	time_t m_StartTime;

	unsigned int m_ElapsedHours;
	unsigned int m_ElapsedMinutes;
	unsigned int m_ElapsedSeconds;

	double m_ElapsedTimeInSecondsOnly;
	double m_ElapsedCPUTimeInSecondsOnly;

	bool m_Started;


}; /* class SimpleTimer */


#endif /* _SimpleTimer_h */
