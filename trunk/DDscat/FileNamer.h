#pragma once

class FileNamer
{
protected:
	static FileNamer *item;
	string cflpol1, cflpol2, cflsca, cflavg, cfle1, cfle2, cfleb1, cfleb2, cflfml;
	int norichar, cashedIwav, cashedIrad, cashedIdir;

protected:
	FileNamer(void);
	virtual ~FileNamer(void);

public:
	static FileNamer *GetInstance();
	static void Kill(void);

public:
	void Namer(int iwav, int irad, int idir);
	void Init(int nr);

public:
	inline const char *GetCflpol1() { return cflpol1.c_str(); }
	inline const char *GetCflpol2() { return cflpol2.c_str(); }
	inline const char *GetCflsca() { return cflsca.c_str(); }
	inline const char *GetCflavg() { return cflavg.c_str(); }
	inline const char *GetCfle1() { return cfle1.c_str(); }
	inline const char *GetCfle2() { return cfle2.c_str(); }
	inline const char *GetCfleb1() { return cfleb1.c_str(); }
	inline const char *GetCfleb2() { return cfleb2.c_str(); }
	inline const char *GetCflfml() { return cflfml.c_str(); }
};
