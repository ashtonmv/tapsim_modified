#ifndef TETGEN_PREDICATES_CXX
#define TETGEN_PREDICATES_CXX

#define TETLIBARY

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#ifdef CPU86
#include <float.h>
#endif

#ifdef LINUX
#include <fpu_control.h>
#endif

#include "dist_tetgen.h"

namespace TetGen
{
	#include "./dist_1.4.3/predicates.cxx"
}

#endif

