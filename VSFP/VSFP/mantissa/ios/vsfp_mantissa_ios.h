// This file provides additional functions for VSFP_MANTISSA class to be
// comfortably used with standard input/output library (iostream).

// It is used only for debugging purposes and is to be excluded from release version.

#pragma once

#include <iostream>
#include <iomanip>
using namespace std;

#include "..\vsfp_mantissa.h"

ostream& operator<<(ostream& ostr, const VSFP_MANTISSA& vm);
