#pragma once

class Line
{
public:
	static int nx, ny, nz;
	int ix, iy, iz, icc1, icc2, icc3;

public:
	Line();
	~Line();
	int Index() const
	{
		return ((ix * ny + iy) * nz + iz);
	}
	bool operator<(const Line &op) const
	{
		return Index() < op.Index();
	}
};
