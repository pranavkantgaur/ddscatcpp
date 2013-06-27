#include "StdAfx.h"
#include "AbstractSolver.h"

AbstractSolver *AbstractSolver::solver = NULL;
SolMethod AbstractSolver::method = SolMethod_End;
map<SolMethod, SolverCreator> *AbstractSolver::factory = NULL;
real AbstractSolver::errscal = (real)0.;

const Complex coner((real)1., (real)0.);
const Complex czero((real)0., (real)0.);

AbstractSolver::AbstractSolver(void)
{
	method = SolMethod_End;
	solver = NULL;
	wk = NULL;
}

AbstractSolver::~AbstractSolver(void)
{
	method = SolMethod_End;
	solver = NULL;
	CleanDelete2(wk);
}

AbstractSolver *AbstractSolver::GetInstance(SolMethod mmethod)
{
	if (method == mmethod)
		return solver;
	else
	{
		if (solver)
			delete solver;
		solver = NULL;
		map<SolMethod, SolverCreator>::iterator ia = factory->find(mmethod);
		if (ia != factory->end())
		{
			solver = ia->second();
			method = mmethod;
		}
		return solver;
	}
}

void AbstractSolver::Kill(void)
{
	if (solver)
	{
		delete solver;
		if (factory)
			delete factory;
	}	
}

bool AbstractSolver::RegisterSolver(SolMethod method, SolverCreator creator)
{
	if (!factory)
		factory = new map<SolMethod, SolverCreator>;
	map<SolMethod, SolverCreator>::iterator ia = factory->find(method);
	if (ia != factory->end())
		return false;
	else
	{
		factory->insert(pair<SolMethod, SolverCreator>(method, creator));
		return true;
	}
}

void AbstractSolver::DeleteSolver()
{
	if (solver)
	{
		delete solver;
		solver = NULL;
		method = SolMethod_End;
	}	
}

void AbstractSolver::SetMatvecFunction(void (*Matvec)(Complex *, Complex *, int *))
{
	this->Matvec = Matvec;
}

Complex AbstractSolver::DotProduct(Complex *cx, Complex *cy, int n)
{
	Complex res;
	for(int i=0; i<n; ++i)
	{
		res += cx[i].conjg() * cy[i];
	}
	return res;
}