#ifndef GUARD_bd_h
#define GUARD_bd_h

#include "rng.h"
#include "info.h"
#include "tree.h"
#include "logging.h"

#ifdef MPIBART
bool bd(tree& x, xinfo& xi, pinfo& pi, RNG& gen, size_t numslaves);
#else
bool bd(tree& x, xinfo& xi, dinfo& di, double* phi, pinfo& pi, RNG& gen, Logger logger);
#endif

#endif
