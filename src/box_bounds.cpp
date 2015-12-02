#include "iset.h"

extern "C" {

SEXP get_box_bounds (SEXP ids, SEXP reflevels, SEXP adex, SEXP chrs, SEXP starts, SEXP ends) try {

    if (!isInteger(ids)) { throw std::runtime_error("grouping vector should be integer"); }
    if (!isInteger(adex)) { throw std::runtime_error("anchor index vector should be integer"); }
    const int npts=LENGTH(ids);
    if (LENGTH(adex)!=npts) { throw std::runtime_error("anchor index and grouping vectors are not of same length"); }
    const int * iptr=INTEGER(ids);
    const int * aptr=INTEGER(adex);

    if (!isString(reflevels)) { throw std::runtime_error("grouping reference levels should be character"); }
    const int nlevels=LENGTH(reflevels);

    if (!isInteger(chrs)) { throw std::runtime_error("chromosome index vector should be integer"); }
    if (!isInteger(starts)) { throw std::runtime_error("start vector should be integer"); }
    if (!isInteger(ends)) { throw std::runtime_error("end vector should be integer"); }
    const int nregs=LENGTH(chrs);
    if (nregs!=LENGTH(starts) || nregs!=LENGTH(ends)) { throw std::runtime_error("chromosome/start/end vectors should have the same length"); }
    const int * cptr=INTEGER(chrs);
    const int * sptr=INTEGER(starts);
    const int * eptr=INTEGER(ends);

    int ngroups=(npts > 0 ? 1 : 0);
    for (int i=1; i<npts; ++i) {
        if (iptr[i]!=iptr[i-1]) { ++ngroups; }
    }

    SEXP output=PROTECT(allocVector(VECSXP, 4));
try {
    SET_VECTOR_ELT(output, 0, allocVector(INTSXP, ngroups));
    int* fac_ptr=INTEGER(VECTOR_ELT(output, 0));
    SET_VECTOR_ELT(output, 1, allocVector(INTSXP, ngroups));
    int* chr_ptr=INTEGER(VECTOR_ELT(output, 1));
    SET_VECTOR_ELT(output, 2, allocVector(INTSXP, ngroups));
    int* start_ptr=INTEGER(VECTOR_ELT(output, 2));
    SET_VECTOR_ELT(output, 3, allocVector(INTSXP, ngroups));
    int* end_ptr=INTEGER(VECTOR_ELT(output, 3));

    if (ngroups > 0) {
        fac_ptr[0]=iptr[0];
        chr_ptr[0]=cptr[aptr[0]];
        start_ptr[0]=sptr[aptr[0]];
        end_ptr[0]=eptr[aptr[0]];
    }

    int curgroup=0;
    for (int i=1; i<npts; ++i) {
        const int& curdex=aptr[i];
        if (curdex < 0 || curdex >= nregs) { throw std::runtime_error("anchor index is out of bounds"); }

        if (iptr[i]!=iptr[i-1]) {
            ++curgroup;
            if (curgroup >= ngroups) { throw std::runtime_error("internal indexing error for groups"); }
            fac_ptr[curgroup]=iptr[i];
            chr_ptr[curgroup]=cptr[curdex];
            start_ptr[curgroup]=sptr[curdex];
            end_ptr[curgroup]=eptr[curdex];
        } else {
            if (chr_ptr[curgroup]!=cptr[curdex]) {
                if (iptr[i] >= nlevels) {
                    throw std::runtime_error("grouping index outside range of reference levels");
                }
                std::stringstream err_out;
                err_out << "multiple chromosomes for group '" << CHAR(STRING_ELT(reflevels, iptr[i])) << "'";
                throw std::runtime_error(err_out.str());
            }
            if (start_ptr[curgroup] > sptr[curdex]) { start_ptr[curgroup] = sptr[curdex]; }
            if (end_ptr[curgroup] < eptr[curdex]) { end_ptr[curgroup] = eptr[curdex]; }
        }
    }
} catch (std::exception& e) {
    UNPROTECT(1);
    throw;
}
    UNPROTECT(1);
    return output;
} catch (std::exception &e) {
    return mkString(e.what());
}


}
