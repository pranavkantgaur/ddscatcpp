//-------------------------------------------------------------------------------------------------------------

#include "includes.h"
#include "function.h"
#include "callstack.h"
#include "symbol.h"
#include "profiler.h"
//-------------------------------------------------------------------------------------------------------------

void _cdecl enterFunc( long retAddress)
{
	Profiler::getProfiler().enterFunc( retAddress);
}
//-------------------------------------------------------------------------------------------------------------

void _cdecl exitFunc( long retAddress)
{
	Profiler::getProfiler().exitFunc( retAddress);
}
//-------------------------------------------------------------------------------------------------------------

extern "C" void __declspec(naked) _cdecl _penter( void ) 
{
	_asm 
	{
		//Prolog instructions
		pushad
		//calculate the pointer to the return address by adding 4*8 bytes 
		//(8 register values are pushed onto stack which must be removed)
        mov  eax, esp
        add  eax, 32
         // retrieve return address from stack
        mov  eax, dword ptr[eax]
        // subtract 5 bytes as instruction for call _penter is 5 bytes long on 32-bit machines, e.g. E8 <00 00 00 00>
        sub  eax, 5
        // provide return address to recordFunctionCall
		push eax
		call enterFunc
		pop eax
		
		//Epilog instructions
		popad
		ret
    }
}
//-------------------------------------------------------------------------------------------------------------

extern "C" void __declspec(naked) _cdecl _pexit( void ) 
{
   _asm 
   {
		//Prolog instructions
		pushad
	
		//calculate the pointer to the return address by adding 4*7 bytes 
		//(7 register values are pushed onto stack which must be removed)
        mov  eax, esp
        add  eax, 28
 
        // retrieve return address from stack
        mov  eax, dword ptr[eax]
 
        // subtract 5 bytes as instruction for call _penter is 5 bytes long on 32-bit machines, e.g. E8 <00 00 00 00>
        sub  eax, 5
 
        // provide return address to recordFunctionCall
		push eax
		call exitFunc
		pop eax
		
		//Epilog instructions
		popad
		ret
    }
}
//-------------------------------------------------------------------------------------------------------------

Profiler &Profiler::getProfiler()
{
	static Profiler profiler;
	return profiler;
}
//-------------------------------------------------------------------------------------------------------------

Profiler::Profiler() 
:	mVAFuncMap(), 
	mhModSymbolMap(),
	mCallStack() 
{
	::CoInitialize(NULL);
}
//-------------------------------------------------------------------------------------------------------------

Profiler::~Profiler() 
{ 
	serialize();
	clear();
	::CoUninitialize();
}
//-------------------------------------------------------------------------------------------------------------

void Profiler::enterFunc( DWORD va)
{
	Function *func = getFunc( va);
	mCallStack.push( func);
	func->startTimer();
}
//-------------------------------------------------------------------------------------------------------------

void Profiler::exitFunc( DWORD va)
{
	Function *func = mCallStack.top();
	mCallStack.pop();
	
	float time = 0.0;
	func->stopTimer( time);
	func->addTime( time);

	func = mCallStack.top();
	if ( NULL != func)
	{
		func->addChildTime( time);
	}
}
//-------------------------------------------------------------------------------------------------------------

void Profiler::clear()
{
	for ( std::map<DWORD, Function*>::iterator vafIter = mVAFuncMap.begin(); 
					vafIter != mVAFuncMap.end(); vafIter++)
	{
		delete vafIter->second;
	}
	mVAFuncMap.clear();

	for ( std::map<std::pair<DWORD,DWORD>,Symbol*>::iterator hIter = mhModSymbolMap.begin(); 
					hIter != mhModSymbolMap.end(); hIter++)
	{
		delete hIter->second;
	}
	mhModSymbolMap.clear();
}
//-------------------------------------------------------------------------------------------------------------

struct FuncSorter
{
	bool operator()( const Function* func1, const Function* func2) const
	{
		return ( 
				(func1->getTotalTime() - func1->getChildTime())  
							> 
				(func2->getTotalTime() - func2->getChildTime())
				);
	}
};
//-------------------------------------------------------------------------------------------------------------

void Profiler::serialize()
{
//	TCHAR logfile[MAX_PATH];
//	if( 0 != GetEnvironmentVariable( _T("PROFILER_LOG"), logfile, MAX_PATH))
//	{
//		std::wstring strLogFile( logfile);
		
		std::wofstream ofLogFile( "d:\\Profiler.log", std::ios::out|std::ios::trunc);
		if( ofLogFile)
		{
			std::set<Function*> entrySet; 
			for (std::map<DWORD, Function*>::iterator iter = mVAFuncMap.begin(); iter != mVAFuncMap.end(); iter++)
				entrySet.insert( iter->second);
			
			ofLogFile << _T("Function Name,Calls,Total Time,Child Time,Self Time\n") << std::endl;
			for (std::set<Function*>::const_iterator funcItr = entrySet.begin(); funcItr != entrySet.end(); funcItr++)
			{
				Function* func = *funcItr;
				ofLogFile << func->getName() << (",") 
						  << func->getCalls() << (",")
						  << func->getTotalTime() << (",")	
						  << func->getChildTime() << (",")
						  << (func->getTotalTime() - func->getChildTime()) << std::endl;
				ofLogFile.flush();							
			}
			ofLogFile.close();
		}
//	}
}
//-------------------------------------------------------------------------------------------------------------

Function*	Profiler::getFunc( DWORD va)
{
	std::map<DWORD,Function*>::iterator vaItr = mVAFuncMap.find( va);
	//If function exists in the map of Virtual Address vs Function object then use it
	if( vaItr != mVAFuncMap.end())
	{
		return vaItr->second;
	}
	else
	{
		Symbol* pSymbol = NULL;
		std::map<std::pair<DWORD,DWORD>,Symbol*>::iterator sItr = mhModSymbolMap.begin();
		std::map<std::pair<DWORD,DWORD>,Symbol*>::iterator eItr = mhModSymbolMap.end();	
		for( ; sItr != eItr; sItr++)
		{
			DWORD baseAddress = sItr->first.first;
			DWORD endAddress = sItr->first.second;
			if( baseAddress <= va && va <= endAddress)//Create the symbol table for this module
			{
				pSymbol = sItr->second;
				break;
			}
		}
		if( NULL == pSymbol)
		{
			HANDLE hProc = GetCurrentProcess();
			HMODULE hMods[1024];
			DWORD cbNeeded;
			if( EnumProcessModules( hProc, hMods, sizeof(hMods), &cbNeeded))
			{
				for ( unsigned int i = 0; i < (cbNeeded / sizeof(HMODULE)); i++ )
				{
					MODULEINFO modinfo;
					ZeroMemory( &modinfo, sizeof(MODULEINFO));
					if( GetModuleInformation( hProc, hMods[i], &modinfo, sizeof(MODULEINFO)))
					{
						DWORD loadAddress = (DWORD)(modinfo.lpBaseOfDll);
						DWORD endAddress = loadAddress + modinfo.SizeOfImage;
						if( loadAddress <= va && va <= endAddress)//Create the symbol table for this module
						{
							pSymbol = new Symbol( hProc, hMods[i]);
							mhModSymbolMap.insert( std::make_pair( std::make_pair(loadAddress,endAddress),pSymbol));
							break;
						}
					}
				}
			}
		}
		if( NULL != pSymbol)
		{
			std::wstring strFuncName = pSymbol->getFuncName( va);
			Function* pFunc = new Function( strFuncName.c_str());
			mVAFuncMap.insert( std::make_pair(va,pFunc));
			return pFunc;
		}
		else
			return NULL;
	}
	return NULL;
}
//-------------------------------------------------------------------------------------------------------------
