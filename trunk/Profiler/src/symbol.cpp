//-------------------------------------------------------------------------------------------------------------

#include "includes.h"
#include "symbol.h"
//-------------------------------------------------------------------------------------------------------------

Symbol::Symbol( HANDLE hProcess, HMODULE hModule)
:	mIsValid ( false),
	mPtrSession ( NULL), 
	mPtrSource ( NULL)
{
	HRESULT hr = S_OK;
	TCHAR szModName[MAX_PATH];
    // Get the full path to the module's file.
    if ( GetModuleFileNameEx( hProcess, hModule, szModName, sizeof( szModName)/ sizeof( TCHAR)))
    {
		hr = mPtrSource.CoCreateInstance( CLSID_DiaSource, NULL, CLSCTX_INPROC_SERVER);
		if (SUCCEEDED(hr))
		{
			if ( SUCCEEDED( mPtrSource->loadDataForExe( szModName, NULL, NULL ) ) )
			{
				if ( SUCCEEDED( mPtrSource->openSession( &mPtrSession ) ) ) 
				{
					hr = mPtrSession->put_loadAddress( (ULONGLONG)hModule);
					if( SUCCEEDED(hr))
						mIsValid = true;
				}
			}
		}    
    }
}
//-------------------------------------------------------------------------------------------------------------

Symbol::~Symbol()
{
	mPtrSession.Release();
	mPtrSource.Release();
}
//-------------------------------------------------------------------------------------------------------------

std::wstring Symbol::getFuncName( long va)
{
	std::wstring strFuncName;
	if( mIsValid)
	{
		IDiaSymbol* pFunc;
		if( SUCCEEDED( mPtrSession->findSymbolByVA( va, SymTagFunction, &pFunc)))
		{
			_bstr_t funcName;
			if( NULL != pFunc && SUCCEEDED( pFunc->get_name(funcName.GetAddress())))
			{
				strFuncName.append((TCHAR*)funcName);
			}
		}
	}
	return strFuncName;
}
//-------------------------------------------------------------------------------------------------------------

