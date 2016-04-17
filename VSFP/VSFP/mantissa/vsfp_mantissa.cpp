#include <stdexcept>

#include "vsfp_mantissa.h"

VSFP_MANTISSA::VSFP_MANTISSA(U16 size)
{
	if (!size)
		throw std::logic_error("Zero size value.");
		
	mantissa = new U64[size];
	mantissa_size = size;
	zero();
}
VSFP_MANTISSA::VSFP_MANTISSA(U16 size, float arg)
{
	if (!size)
		throw std::logic_error("Zero size value.");

	mantissa = new U64[size];
	mantissa_size = size;
	zero();
	mantissa[0] = static_cast<U64>((*reinterpret_cast<U32*>(&arg)) << 9) << 8*sizeof(float);
}
VSFP_MANTISSA::VSFP_MANTISSA(U16 size, double arg)
{
	if (!size)
		throw std::logic_error("Zero size value.");

	mantissa = new U64[size];
	mantissa_size = size;
	zero();
	mantissa[0] = (*reinterpret_cast<U64*>(&arg)) << 12;
}
VSFP_MANTISSA::~VSFP_MANTISSA()
{
	delete[] mantissa;
}

void VSFP_MANTISSA::zero()
{
	for (U16 i = 0; i < mantissa_size; i++)
		mantissa[i] = 0;
	return;
}
void VSFP_MANTISSA::change_size(U16 new_size)
{
	if (new_size == mantissa_size)
		throw std::logic_error("New size is the same as the old size.");

	U64* new_mantissa = new U64[new_size];

	U16 i = 0;
	if (new_size > mantissa_size)
	{
		for (; i < mantissa_size; i++)
			new_mantissa[i] = mantissa[i];
		for (; i < new_size; i++)
			new_mantissa[i] = 0;
	}
	else for (; i < new_size; i++)
			new_mantissa[i] = mantissa[i];

	mantissa_size = new_size;
	return;
}

VSFP_MANTISSA& VSFP_MANTISSA::operator<< (U16 shift)
{
	if (!shift)
		throw std::logic_error("Zero shift value.");

	U16 i = 0;
	if (shift/64)
	{
		for (; i+shift/64 < mantissa_size; i++)
			mantissa[i] = mantissa[i+shift/64];
		for (; i < mantissa_size; i++)
			mantissa[i] = 0;
		for (i = 0; i+1+shift/64 < mantissa_size; i++)
			mantissa[i] = (mantissa[i] << shift%64) | (mantissa[i+1] >> (64-shift%64));
		mantissa[i] <<= shift%64;
	}
	else
	{
		for (; i+1 < mantissa_size; i++)
			mantissa[i] = (mantissa[i] << shift) | (mantissa[i+1] >> (64-shift));
		mantissa[i] <<= shift;
	}

	return *this;
}
VSFP_MANTISSA& VSFP_MANTISSA::operator>> (U16 shift)
{
	if (!shift)
		throw std::logic_error("Zero shift value.");

	U16 i = mantissa_size-1;
	if (shift/64)
	{
		for (; ;i--)
		{	mantissa[i] = mantissa[i-shift/64];
			if (!(i-shift/64))
				break;	}
		for (i--; ;i--)
		{	mantissa[i] = 0;
			if (!i)
				break;	}
		for (i = mantissa_size-1; i > 0; i--)
			mantissa[i] = (mantissa[i] >> shift) | (mantissa[i-1] << (64-shift));
		mantissa[i] >>= shift;
	}
	else
	{
		for (; i > 0; i--)
			mantissa[i] = (mantissa[i] >> shift) | (mantissa[i-1] << (64-shift));
		mantissa[i] >>= shift;
	}

	return *this;
}

U8 VSFP_MANTISSA::operator+= (const VSFP_MANTISSA& arg)
{
//	for (U16 i = min(mantissa_size, arg.mantissa_size)-1; ; i--)
	{

	}
	return 0;
}