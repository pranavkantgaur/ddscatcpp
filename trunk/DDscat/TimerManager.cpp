#include "StdAfx.h"

#include "TimerManager.h"
#include "DDscatMain.h"

TimerManager * TimerManager::item = NULL;
TimerManager::TimerManager(void)
{
	ntimers = NTIMERS;
	timers = new real[ntimers];
}

TimerManager::~TimerManager(void)
{
	CleanDelete2(timers);
}

TimerManager *TimerManager::GetInstance()
{
	if (!item)
		item = new TimerManager;

	return item;
}

void TimerManager::Kill()
{
	CleanDelete(item);
}

void TimerManager::SetValue(int index, real value)
{
	if (index >= 0 && index < ntimers)
	{
		timers[index] = value;
	}
}

real TimerManager::GetValue(int index)
{
	if ((index >= 0) && (index < ntimers))
	{
		return timers[index];		
	}
	else
		return -(real)1.;
}
