#ifndef ISET_H
#define ISET_H

#include <deque>
#include <algorithm>
#include <stdexcept>
#include <set>
#include <cstring>
#include <sstream>

#include "R.h"
#include "Rinternals.h"
#include "R_ext/Rdynload.h"
#include "R_ext/Visibility.h"

extern "C" {

SEXP linear_olaps(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);

SEXP paired_olaps(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);

SEXP expand_pair_links(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);

SEXP get_box_bounds(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);

}

#endif
