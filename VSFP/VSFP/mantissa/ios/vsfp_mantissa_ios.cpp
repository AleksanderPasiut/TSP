#include "vsfp_mantissa_ios.h"

ostream& operator<<(ostream& ostr, const VSFP_MANTISSA& vm)
{
	ostr << hex;
	for (U16 i = 0; i < vm.size(); i++)
			ostr << setfill('0') << setw(16) << vm.mantissa[i] << " ";

	ostr << endl;
	return ostr;
}