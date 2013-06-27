#ifndef CleanDelete
	#define CleanDelete(x) {if(x)delete(x);(x)=NULL;}
	#define CleanDelete2(x) {if(x)delete[](x);(x)=NULL;}
#endif
