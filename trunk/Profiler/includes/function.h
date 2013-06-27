//-------------------------------------------------------------------------------------------------------------

#pragma once
//-------------------------------------------------------------------------------------------------------------

//! Class defining the function object of the target application.

//! Class defining the function object of the target application.
class Function
{
public:
	
	//! Constructor.
	//! @param name the name of the function.
	Function( const TCHAR* name);

	//! Function to add the specified time to child function time
	//! to child function execution time.
	//! @param time the child function execution time.
	void addChildTime( float time);

	//! Function to add the specified time to total time.
	//! Total time includes self execution time and child function execution time.
	//! @param time the total time of execution for the function.
	void addTime( float time);

	//! Function to get total execution time of the function.
	//! @returns a total time of the function execution.
	float getTotalTime() const;

	//! Function to get execution time of the child function(s).
	//! @returns a total time of the function execution.
	float getChildTime() const;

	//! Function to get the function name.
	//! @returns the function name.
	const TCHAR* getName() const;
	
	//! Function to get the number of times this function is called.
	//! @returns number of times this function is called.
	int	 getCalls() const;

	//! Function to start the timer for logging the function execution time.
	void startTimer();
	
	//! Function to stop the timer for logging the function execution time.
	//! @param total time elapsed since timer is started.
	void stopTimer( float& time);

protected:
	std::wstring	mFuncName;
	int				mNumCalls;
	float			mTotalTime;
	float			mChildTime;
	clock_t			mStartTime;	//Timer start time
};
//-------------------------------------------------------------------------------------------------------------

