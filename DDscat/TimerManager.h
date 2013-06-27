#ifndef __TIMEMANAGER_H__
#define __TIMEMANAGER_H__

#include "Definitions.h"
#define NTIMERS 12

class TimerManager
{
protected:
	static TimerManager *item;
	int ntimers;
	real *timers;
	TimerManager(void);
	virtual ~TimerManager(void);

public:
	static TimerManager *GetInstance();
	static void Kill();

public:
	void SetValue(int index, real value);
	real GetValue(int index);
	int GetTimerNumber() { return ntimers; }
};

#endif // __TIMEMANAGER_H__

//
// Timers are for:
// (0): total GPU time
// (1): # of iterations
// (2): not used
// (3-6): Scat, Evale, Alphadiag, Evalq
// (7-10): Scat, Evale, Alphadiag, Evalq again
// (11): Scat

