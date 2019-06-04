#ifndef GUARD_bd_h
#define GUARD_bd_h

#include "rng.h"
#include "info.h"
#include "tree.h"

#ifdef MPIBART
bool bd(tree& x, xinfo& xi, pinfo& pi, RNG& gen, size_t numslaves);
#else
bool bd(tree& x, xinfo& xi, dinfo& di, pinfo& pi, RNG& gen);
bool bdprec(tree& x, xinfo& xi, dinfo& di, pinfo& pi, RNG& gen);
bool bdhet(tree& x, xinfo& xi, dinfo& di, double* phi, pinfo& pi, RNG& gen);
bool bd_rj(tree& x, xinfo& xi, dinfo& di, pinfo& pi, RNG& gen);
#endif

#endif
