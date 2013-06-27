//-------------------------------------------------------------------------------------------------------------

#pragma once
//-------------------------------------------------------------------------------------------------------------

//! Function to provide _penter function hook as a result of /Gh compiler option.
//! This function will be called at the start of every method or function.
extern "C" void PROFILER_API _cdecl _penter( void );

//! Function to provide _pexit function hook as a result of /GH compiler option.
//! This function will be called at the end of every method or function.
extern "C" void PROFILER_API _cdecl _pexit( void );
//-------------------------------------------------------------------------------------------------------------

class Function;
class CallStack;
class Symbol;
//-------------------------------------------------------------------------------------------------------------

//! Class defining the profiler to log function execution time.

//! Class defining the profiler to log function execution time.
class Profiler
{
public:

	//! Function to get the instance of a profiler object.
	//! @returns the Profiler object.
	static Profiler &getProfiler();
	
	//! Destructor.
	~Profiler();

	//! Function to start profiling after entering the target function.
	//! @param va the virtual memory address from the target function.
	void enterFunc( DWORD va);
	
	//! Function to stop profiling before exiting the target function.
	//! @param va the virtual memory address from the target function.
	void exitFunc( DWORD va);
	
	//! Function to clear all the profiler counters.
	void clear();
	
	//! Function to serialize the profiler data to a log file.
	//! Log file is defined using "PROFILER_LOG" environment variable.
	void serialize();

protected:

	//! Constructor.
	Profiler();

	//! Function to get the function object from a virtual address.
	//! @param va the virtual memory address from the target function.
	Function *getFunc( DWORD va);

	std::map<DWORD,Function *> mVAFuncMap;
	std::map<std::pair<DWORD,DWORD>,Symbol*> mhModSymbolMap;
	CallStack mCallStack;
};
//-------------------------------------------------------------------------------------------------------------
