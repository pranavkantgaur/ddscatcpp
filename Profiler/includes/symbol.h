//-------------------------------------------------------------------------------------------------------------

#pragma once
//-------------------------------------------------------------------------------------------------------------

//! Class defining the symbol table for a loaded module of an application.

//! Class defining the symbol table for a loaded module of an application.
class Symbol
{
public:

	//! Constructor.
	//! @param hProcess the handle to the target application process.
	//! @param hModule the handle to the loaded application module.
	Symbol( HANDLE hProcess, HMODULE hModule);
	
	//! Destructor.
	~Symbol();

	//! Function to get the function name from the virtual memory address.
	//! The virtual memory address should be from the function body either start or end of function.
	//! @returns the function name.
	std::wstring getFuncName( long va);
	
protected:
	CComPtr<IDiaDataSource> mPtrSource;
	CComPtr<IDiaSession>	mPtrSession;
	bool					mIsValid;
};
//-------------------------------------------------------------------------------------------------------------