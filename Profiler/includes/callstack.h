//-------------------------------------------------------------------------------------------------------------

#pragma once
//-------------------------------------------------------------------------------------------------------------

class Function;
//-------------------------------------------------------------------------------------------------------------

//! Class defining the function call stack of the target application.

//! Class defining the function call stack of the target application.
class CallStack
{
public:

	//! Constructor.
	CallStack();

	//! Function to push the function object at the start of it's execution.
	//! @param ent the pointer to the function object.
	void		push( Function *ent);
	
	//! Function to pop the function object at the end of it's execution.
	void		pop();
	
	//! Function to retrieve the first function on the call stack.
	//! @returns a Function pointer at the top of the call stack.
	Function*	top();

	//! Function to test if the call stack is empty.
	//! @returns a true if the call stack is not empty else false.
	bool		hasChild();

protected:
	std::stack<Function*>	mEntryStack;
};
//-------------------------------------------------------------------------------------------------------------