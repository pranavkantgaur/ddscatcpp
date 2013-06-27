#ifndef __ARRAYF_H__
#define __ARRAYF_H__

#include <cstdlib>
#ifdef _WIN32
	#include <malloc.h>
#endif	

// The model of Fortran array
template <typename T>
class AbstractArrayF
{
protected:
	T *data;
	unsigned int *sizes;
	unsigned int dimension;
	unsigned int totalSize;		// in items

public:
	AbstractArrayF(void) : data(NULL), sizes(NULL), dimension(0), totalSize(0)  { }
	~AbstractArrayF(void) 
	{ 
		Deallocate();
	}

public:
	inline unsigned int GetTotalSize() { return totalSize; }
	inline T *GetData() { return data; }
	inline unsigned int GetDimension() { return dimension; }
	unsigned int GetSize(unsigned int index)
	{
		return (index < dimension) ? sizes[index] : 0;
	}
	void Clear()
	{
		if (totalSize > 0)
		{
			memset(data, 0, totalSize * sizeof(T));
		}
	}
	void Deallocate()
	{
		if (data) { free(data); data = NULL; } 
		if (sizes) { free(sizes); sizes = NULL; }
	}
};

// Dimensions and accessors to the array s.b. given in Fortran way
// you don't need to change the fortran indexes when access.
// Anyway, the data is places in C way: the last index changes first.
template <typename T>
class Array2F : public AbstractArrayF<T>
{
public:
	Array2F(void) { this->dimension = 2; }
	~Array2F() { }
	T *Dimension(unsigned int xSize, unsigned int ySize)
	{
		this->sizes = (unsigned int *)calloc(this->dimension, sizeof(unsigned int));
		this->sizes[0] = xSize;
		this->sizes[1] = ySize;
		this->totalSize = xSize * ySize;
		this->data = (T *)calloc(this->totalSize, sizeof(T));
		return this->data;
	}
	T *Extend(unsigned int delta)
	{
		this->sizes[0] += delta;
		this->totalSize = this->sizes[0] * this->sizes[1];
		this->data = (T *)realloc(this->data, this->totalSize * sizeof(T));
		return this->data;
	}
	T *Close(unsigned int size)
	{
		if (size != this->sizes[0])
		{
			this->sizes[0] = size;
			this->totalSize = this->sizes[0] * this->sizes[1];
			this->data = (T *)realloc(this->data, this->totalSize * sizeof(T));
		}
		return this->data;
	}
	T *CopyItself(unsigned int toto, unsigned int from, unsigned int cnt)
	{
		for(unsigned int i=0; i<cnt; ++i)
		{
            for(unsigned int j=0; j<this->sizes[1]; ++j)
			{
                this->data[(toto + i) * this->sizes[1] + j] = this->data[(from + i) * this->sizes[1] + j];
			}
		}
		return this->data;
	}
	T *MoveItself(unsigned int toto, unsigned int from, unsigned int cnt)			// ChB: experiment
	{
		unsigned int rowLength = this->sizes[1] * sizeof(T);
		for(unsigned int i=0; i<cnt; ++i)
		{
			memmove(this->data + (toto + i) * rowLength, this->data + (from + i) * rowLength, rowLength);
		}
		return this->data;
	}
	T &Value(unsigned int x, unsigned int y)
	{
		return this->data[x * this->sizes[1] + y];
	}
	T *Row(unsigned int x)
	{
		return this->data + x * this->sizes[1];
	}
	bool Fill3(int indX, T a, T b, T c)
	{
		if (this->sizes[1] != 3)
			return false;
		unsigned int pos = indX * this->sizes[1];
		this->data[pos] = a;
		this->data[pos + 1] = b;
		this->data[pos + 2] = c;
		return true;
	}
	bool Add3(int indX, T a, T b, T c)
	{
		if (this->sizes[1] != 3)
			return false;
		unsigned int pos = indX * this->sizes[1];
		this->data[pos] += a;
		this->data[pos + 1] += b;
		this->data[pos + 2] += c;
		return true;
	}
	bool Fill3(int indX, T a)
	{
		if (this->sizes[1] != 3)
			return false;
		unsigned int pos = indX * this->sizes[1];
        	this->data[pos] = this->data[pos + 1] = this->data[pos + 2] = a;
		return true;
	}
	bool Add3(int indX, T a)
	{
		if (this->sizes[1] != 3)
			return false;
		unsigned int pos = indX * this->sizes[1];
		this->data[pos] += a;
		this->data[pos + 1] += a;
		this->data[pos + 2] += a;
		return true;
	}
	void PrintLine(FILE *file, const char *Format, int line)
	{
		unsigned int pos = line * this->sizes[1];
		for(unsigned int i=0; i<this->sizes[1]; ++i)
		{
			fprintf(file, Format, this->data[pos+i]);
		}
	}
   	void ReadLine(int file, int line)
	{
		read(file, Row(line), this->sizes[1] * sizeof(T));
	}
	void WriteLine(int file, int line)
	{
		write(file, Row(line), this->sizes[1] * sizeof(T));
	}
};

template <typename T>
class Array3F : public AbstractArrayF<T>
{
public:
	Array3F(void) { this->dimension = 3; }
	~Array3F(void) { }
	T *Dimension(unsigned int xSize, unsigned int ySize, unsigned int zSize)
	{
		this->sizes = (unsigned int *)calloc(this->dimension, sizeof(unsigned int));
		this->sizes[0] = xSize;
		this->sizes[1] = ySize;
		this->sizes[2] = zSize;
		this->totalSize = xSize * ySize * zSize;
		this->data = (T *)calloc(this->totalSize, sizeof(T));
		return this->data;
	}
	T &Value(unsigned int x, unsigned int y, unsigned int z)
	{
		return this->data[z + (y + x * this->sizes[1]) * this->sizes[2]];
	}
};

// C-like array with 4 dimensions, organized as stack of Array3F, namely
// Array3F datas are placed in single data array, but the last index 
// points to data of correspondent subarray
template <typename T>
class Array4Stacked : public AbstractArrayF<T>
{
protected:
	unsigned int lastDimStep;

public:
	Array4Stacked(void) { this->dimension = 4; }
	~Array4Stacked(void) { }
	T *Dimension(unsigned int xSize, unsigned int ySize, unsigned int zSize, unsigned int aSize)
	{
		this->sizes = (unsigned int *)calloc(this->dimension, sizeof(unsigned int));
		this->sizes[0] = xSize;
		this->sizes[1] = ySize;
		this->sizes[2] = zSize;
		this->sizes[3] = aSize;
		this->totalSize = xSize * ySize * zSize * aSize;
		this->data = (T *)calloc(this->totalSize, sizeof(T));
		this->lastDimStep = xSize * ySize * zSize;
		return this->data;
	}
	T &Value(unsigned int x, unsigned int y, unsigned int z, unsigned int a)
	{
//		return (this->data + a * this->lastDimStep)[z + (y + x * this->sizes[1]) * this->sizes[2]];			// C-like access
		return (this->data + a * this->lastDimStep)[x + (y + z * this->sizes[1]) * this->sizes[0]];			// Fortran accessor	
	}
	T *SubarrayData(unsigned int index)
	{
		return this->data + this->lastDimStep * index;
	}
//		For F-version of f array
	void InsertLastDim(unsigned int ax, unsigned int ay, unsigned int az, T *f)
	{
		unsigned int pos = ax + this->sizes[0] * (ay + this->sizes[1] * az);
		for(unsigned int i=0; i<this->sizes[3]; ++i)
		{
			this->data[pos] = f[i];
			pos += this->lastDimStep;
		}
	}
//		For F-version of f array
	void ClearIntoPlanes(unsigned int ax, unsigned int ay, unsigned int az)
	{
		unsigned int i, j, k, m, mm, index;
		for(k=0; k<this->sizes[2]; ++k)
		{
			for(j=0; j<this->sizes[1]; ++j)
			{
				index = ax + (j + k * this->sizes[1]) * this->sizes[0];
				for(m=mm=0; m<this->sizes[3]; ++m, mm+=this->lastDimStep)
				{
					this->data[index + mm].clear();
				}	
			}
		}
		for(k=0; k<this->sizes[2]; ++k)
		{
			for(i=0; i<this->sizes[0]; ++i)
			{
				index = i + (ay + k * this->sizes[1]) * this->sizes[0];
				for(m=mm=0; m<this->sizes[3]; ++m, mm+=this->lastDimStep)
				{
					this->data[index + mm].clear();
				}	
			}
		}
		for(j=0; j<this->sizes[1]; ++j)
		{
			for(i=0; i<this->sizes[0]; ++i)
			{
				index = i + (j + az * this->sizes[1]) * this->sizes[0];
				for(m=mm=0; m<this->sizes[3]; ++m, mm+=this->lastDimStep)
				{
					this->data[index + mm].clear();
				}	
			}
		}
	}
	/* **		For C-version of f array
	void InsertLastDim(unsigned int ax, unsigned int ay, unsigned int az, T *f)
	{
		unsigned int pos = az + sizes[2] * (ay + sizes[1] * ax);
		for(unsigned int i=0; i<sizes[3]; ++i)
		{
			data[pos] = f[i];
			pos += lastDimStep;
		}
	}
	** */
	/* **		For C-version of this array
	void ClearIntoPlanes(unsigned int ax, unsigned int ay, unsigned int az)
	{
		unsigned int i, j, k, m, mm, index;
		for(j=0; j<sizes[1]; ++j)
		{
			for(k=0; k<sizes[2]; ++k)
			{
				index = k + (j + ax * sizes[1]) * sizes[2];
				for(m=mm=0; m<sizes[3]; ++m, mm+=lastDimStep)
				{
					data[index + mm].clear();
				}	
			}
		}
		for(i=0; i<sizes[0]; ++i)
		{
			for(k=0; k<sizes[2]; ++k)
			{
				index = k + (ay + i * sizes[1]) * sizes[2];
				for(m=mm=0; m<sizes[3]; ++m, mm+=lastDimStep)
				{
					data[index + mm].clear();
				}	
			}
		}
		for(i=0; i<sizes[0]; ++i)
		{
			for(j=0; j<sizes[1]; ++j)	
			{
				index = az + (j + i * sizes[1]) * sizes[2];
				for(m=mm=0; m<sizes[3]; ++m, mm+=lastDimStep)
				{
					data[index + mm].clear();
				}	
			}
		}
	}
	** */
};

#endif

/* **
template <typename T>
class __declspec(dllexport) Array4Fft : public AbstractArrayF<T>
{
protected:
	unsigned int lastDimStep;

public:
	Array4Fft(void) { dimension = 4; }
	~Array4Fft(void) { }
	T *Dimension(unsigned int xSize, unsigned int ySize, unsigned int zSize, unsigned int aSize)
	{
		sizes = (unsigned int *)malloc(dimension * sizeof(unsigned int));
		sizes[0] = xSize;
		sizes[1] = ySize;
		sizes[2] = zSize;
		sizes[3] = aSize;
		totalSize = xSize * ySize * zSize * aSize;
		data = (T *)malloc(totalSize * sizeof(T));
		lastDimStep = xSize * ySize * zSize;
		return data;
	}
	T &Value(unsigned int x, unsigned int y, unsigned int z, unsigned int a)
	{
		return data[x + (y + (z + a*sizes[2])*sizes[1])*sizes[0]];
	}
	T *SubarrayData(unsigned int index)
	{
		return data + lastDimStep*index;
	}
	void ExtractLastDim(unsigned int ax, unsigned int ay, unsigned int az, T *f, unsigned int num)
	{
		unsigned int pos = ax + sizes[0]*(ay + sizes[1]*az);
		for(unsigned int i=0; i<num; ++i)
		{
			f[i] = data[pos];
			pos += lastDimStep;
		}
	}
	void Reorganize(unsigned int ax, unsigned int ay, unsigned int az, T *tmp)
	{
		unsigned int pos1 = ax + sizes[0]*(ay + sizes[1]*az);
		unsigned int pos2 = pos1 + lastDimStep;
		unsigned int pos3 = pos2 + lastDimStep;
		Complex xx = tmp[0] * data[pos1] + tmp[1] * data[pos2] + tmp[2] * data[pos3];
		Complex yy = tmp[1] * data[pos1] + tmp[3] * data[pos2] + tmp[4] * data[pos3];
		Complex zz = tmp[2] * data[pos1] + tmp[4] * data[pos2] + tmp[5] * data[pos3];
		data[pos1] = xx;
		data[pos2] = yy;
		data[pos3] = zz;
	}
	void SetIntoPlanes(unsigned int ax, unsigned int ay, unsigned int az, const T &tmp)
	{
		unsigned int i, j, k, m, mm, index;
		for(k=0; k<sizes[2]; ++k)
		{
			for(j=0; j<sizes[1]; ++j)
			{
				index = ax + (j + k * sizes[1]) * sizes[0];
				for(m=mm=0; m<sizes[3]; ++m, mm+=lastDimStep)
				{
					data[index + mm] = tmp;
				}	
			}
		}
		for(k=0; k<sizes[2]; ++k)
		{
			for(i=0; i<sizes[0]; ++i)
			{
				index = i + (ay + k * sizes[1]) * sizes[0];
				for(m=mm=0; m<sizes[3]; ++m, mm+=lastDimStep)
				{
					data[index + mm] = tmp;
				}	
			}
		}
		for(j=0; j<sizes[1]; ++j)
		{
			for(i=0; i<sizes[0]; ++i)
			{
				index = i + (j + az * sizes[1]) * sizes[0];
				for(m=mm=0; m<sizes[3]; ++m, mm+=lastDimStep)
				{
					data[index + mm] = tmp;
				}	
			}
		}
	}
	void InsertLastDim(unsigned int ax, unsigned int ay, unsigned int az, T *f, unsigned int num)
	{
		unsigned int pos = ax + sizes[0]*(ay + sizes[1]*az);
		for(unsigned int i=0; i<num; ++i)
		{
			data[pos] = f[i];
			pos += lastDimStep;
		}
	}
};
** */
