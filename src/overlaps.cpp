#include "R.h"
#include "Rinternals.h"
#include <deque>
#include <algorithm>
#include <stdexcept>

extern "C" {

/* Checking validity of inputs */

void check_indices (const int* qsptr, const int* qeptr, const int Nq, const int* sjptr, const int Ns, const int Ns_all) {
    int curqs, curqe;
    for (int checkdex=0; checkdex < Nq; ++checkdex) { 
        curqs = qsptr[checkdex];
        curqe = qeptr[checkdex];
        if (curqs==NA_INTEGER || curqe==NA_INTEGER) { throw std::runtime_error("query indices must be finite integers"); }
        if (curqs >= curqe)  { continue; } // Empty range, no need to check actual values.
        if (curqs >= Ns || curqs < 0) { throw std::runtime_error("query start index out of bounds"); }
        if (curqe > Ns || curqe < 0) { throw std::runtime_error("query end index out of bounds"); }
    }
    
    if (Ns_all < 0) { throw std::runtime_error("total number of subjects must be non-negative"); }
    int curs;
    for (int checkdex=0; checkdex < Ns; ++checkdex) { 
        curs = sjptr[checkdex];
        if (curs >= Ns_all || curs < 0 || curs==NA_INTEGER) { throw std::runtime_error("subject index out of bounds"); }
    }
    return;
}

/* Detects all overlaps between linear ranges and either interacting loci in a pair */

SEXP expand_olaps (SEXP anchor1, SEXP anchor2, SEXP querystarts, SEXP queryends, SEXP subject, SEXP nsubjects) try {
    if (!isInteger(anchor1) || !isInteger(anchor2)) { throw std::runtime_error("anchors must be integer vectors"); }
    const int Npairs = LENGTH(anchor1);
    if (Npairs != LENGTH(anchor2)) { throw std::runtime_error("anchor vectors must be of equal length"); } 
    const int* a1ptr=INTEGER(anchor1), *a2ptr=INTEGER(anchor2);
    
    if (!isInteger(querystarts) || !isInteger(queryends)) { throw std::runtime_error("query indices must be integer vectors"); }
    const int Nq = LENGTH(querystarts);
    if (Nq != LENGTH(queryends)) { throw std::runtime_error("query indices must be of equal length"); }
    const int* qsptr=INTEGER(querystarts), *qeptr=INTEGER(queryends);
    
    if (!isInteger(subject)) { throw std::runtime_error("subject indices must be integer"); }
    const int Ns = LENGTH(subject);
    const int *sjptr=INTEGER(subject);
    if (!isInteger(nsubjects) || LENGTH(nsubjects)!=1) { throw std::runtime_error("total number of subjects must be an integer scalar"); }
    const int Ns_all = asInteger(nsubjects);
   
    // Checking indices. 
    check_indices(qsptr, qeptr, Nq, sjptr, Ns, Ns_all);

    /* Constructing an output deque */
    std::deque<int> new_query, new_subject;
    int* latest_pair=(int*)R_alloc(Ns_all, sizeof(int));
    for (int checkdex=0; checkdex < Ns_all; ++checkdex) { latest_pair[checkdex] = -1; }

    int curpair=0, mode=0, curq=0, curindex=0, curs=0, just_added=0;
    for (curpair=0; curpair<Npairs; ++curpair) {
        just_added=0;

        for (mode=0; mode<2;  ++mode) { 
            if (mode == 0) {
                curq = a1ptr[curpair];
            } else if (curq != a2ptr[curpair]) {
                curq = a2ptr[curpair];
            } else { break; }

            if (curq >= Nq || curq < 0 || curq==NA_INTEGER) { throw std::runtime_error("region index out of bounds"); }
            for (curindex=qsptr[curq]; curindex<qeptr[curq]; ++curindex) {
                curs=sjptr[curindex];
                if (latest_pair[curs] < curpair) { 
                    new_query.push_back(curpair);
                    new_subject.push_back(curs);
                    latest_pair[curs] = curpair;
                    ++just_added;
                }
            }
        }

        // Sorting for output.
        std::sort(new_subject.end()-just_added, new_subject.end());
    }

    // Transplanting the object.
    SEXP output=PROTECT(allocVector(VECSXP, 2));
    try { 
        SET_VECTOR_ELT(output, 0, allocVector(INTSXP, new_query.size()));
        SET_VECTOR_ELT(output, 1, allocVector(INTSXP, new_subject.size()));
        int * new_qptr=INTEGER(VECTOR_ELT(output, 0));
        std::copy(new_query.begin(), new_query.end(), new_qptr);
        int * new_sptr=INTEGER(VECTOR_ELT(output, 1));
        std::copy(new_subject.begin(), new_subject.end(), new_sptr);
    } catch (std::exception& e) {
        UNPROTECT(1);
        throw;
    }
    
    UNPROTECT(1);
    return output;
} catch (std::exception& e) { 
    return mkString(e.what());
}

/* Useful classes to get different outputs from the same code. */

class output_store {
public:
    output_store() {};
    virtual ~output_store() {};
    virtual void prime(int, int)=0;
    virtual void acknowledge(int, int)=0;
    virtual void postprocess()=0;
    virtual SEXP generate() const=0;
};

class expanded_overlap : public output_store { // Stores all new query:subject relations.
public:
    expanded_overlap() : just_added(0) {};
    ~expanded_overlap() {};
    void prime(int nq, int ns) {};
    void acknowledge(int q, int s) {
        new_query.push_back(q);
        new_subject.push_back(s);
        ++just_added;
    }
    void postprocess() {
        std::sort(new_subject.end()-just_added, new_subject.end());
        just_added=0;
    }
    SEXP generate () const {
        SEXP output=PROTECT(allocVector(VECSXP, 2));
        try { 
            SET_VECTOR_ELT(output, 0, allocVector(INTSXP, new_query.size()));
            SET_VECTOR_ELT(output, 1, allocVector(INTSXP, new_subject.size()));
            int * new_qptr=INTEGER(VECTOR_ELT(output, 0));
            std::copy(new_query.begin(), new_query.end(), new_qptr);
            int * new_sptr=INTEGER(VECTOR_ELT(output, 1));
            std::copy(new_subject.begin(), new_subject.end(), new_sptr);
        } catch (std::exception& e) {
            UNPROTECT(1);
            throw;
        }     
        UNPROTECT(1);
        return output;
    } 
private:
    std::deque<int> new_query, new_subject;
    int just_added;
};

class query_overlap : public output_store { // Stores number of times 'query' was hit.
public:
    query_overlap() : nquery(0) {};
    ~query_overlap() {};
    void prime(int nq, int ns) {
        nquery=nq;
        query_hit.resize(nq);    
    };
    void acknowledge(int q, int s) {
        if (q >= nquery) { throw std::runtime_error("requested query index out of range"); }
        ++query_hit[q];
    }
    void postprocess() {}
    SEXP generate () const {
        SEXP output=PROTECT(allocVector(INTSXP, nquery));
        try { 
            int * new_optr=INTEGER(output);
            std::copy(query_hit.begin(), query_hit.end(), new_optr);
        } catch (std::exception& e) {
            UNPROTECT(1);
            throw;
        }     
        UNPROTECT(1);
        return output;
    } 
private:
    int nquery;
    std::deque<int> query_hit;
};

class subject_overlap : public output_store { // Stores number of times 'subject' was hit.
public:
    subject_overlap() : nsubject(0) {};
    ~subject_overlap() {};
    void prime(int nq, int ns) {
        nsubject=ns;
        subject_hit.resize(ns);    
    };
    void acknowledge(int q, int s) {
        if (s >= nsubject) { throw std::runtime_error("requested subject index out of range"); }
        ++subject_hit[s];
    }
    void postprocess() {}
    SEXP generate () const {
        SEXP output=PROTECT(allocVector(INTSXP, nsubject));
        try { 
            int * new_optr=INTEGER(output);
            std::copy(subject_hit.begin(), subject_hit.end(), new_optr);
        } catch (std::exception& e) {
            UNPROTECT(1);
            throw;
        }     
        UNPROTECT(1);
        return output;
    } 
private:
    int nsubject;
    std::deque<int> subject_hit;
};

/* Base function, to detect all overlaps between paired ranges and both interacting loci in a pair */

void detect_paired_olaps(output_store* output, SEXP anchor1, SEXP anchor2, SEXP querystarts1, SEXP queryends1, SEXP subject1, 
        SEXP querystarts2, SEXP queryends2, SEXP subject2, SEXP nsubjects) {

    if (!isInteger(anchor1) || !isInteger(anchor2)) { throw std::runtime_error("anchors must be integer vectors"); }
    const int Npairs = LENGTH(anchor1);
    if (Npairs != LENGTH(anchor2)) { throw std::runtime_error("anchor vectors must be of equal length"); } 
    const int* a1ptr=INTEGER(anchor1), *a2ptr=INTEGER(anchor2);
    
    if (!isInteger(querystarts1) || !isInteger(queryends1)) { throw std::runtime_error("query indices (1) must be integer vectors"); }
    const int Nq = LENGTH(querystarts1);
    if (Nq != LENGTH(queryends1)) { throw std::runtime_error("query indices (1) must be of equal length"); }
    const int* qsptr1=INTEGER(querystarts1), *qeptr1=INTEGER(queryends1);
    
    if (!isInteger(subject1)) { throw std::runtime_error("subject indices (1) must be integer"); }
    const int Ns1 = LENGTH(subject1);
    const int *sjptr1=INTEGER(subject1);

    if (!isInteger(querystarts2) || !isInteger(queryends2)) { throw std::runtime_error("query indices (2) must be integer vectors"); }
    if (Nq != LENGTH(querystarts2) || Nq != LENGTH(queryends2)) { throw std::runtime_error("query indices (2) must be of equal length"); }
    const int* qsptr2=INTEGER(querystarts2), *qeptr2=INTEGER(queryends2);
    
    if (!isInteger(subject2)) { throw std::runtime_error("subject indices (2) must be integer"); }
    const int Ns2 = LENGTH(subject2);
    const int *sjptr2=INTEGER(subject2);

    if (!isInteger(nsubjects) || LENGTH(nsubjects)!=1) { throw std::runtime_error("total number of subjects must be an integer scalar"); }
    const int Ns_all = asInteger(nsubjects);
   
    // Check indices.
    check_indices(qsptr1, qeptr1, Nq, sjptr1, Ns1, Ns_all);
    check_indices(qsptr2, qeptr2, Nq, sjptr2, Ns2, Ns_all);

    // Setting up logging arrays. 
    output->prime(Npairs, Ns_all);
    int* latest_pair_A=(int*)R_alloc(Ns_all, sizeof(int));
    int* latest_pair_B=(int*)R_alloc(Ns_all, sizeof(int));
    bool* is_stored_A=(bool*)R_alloc(Ns_all, sizeof(bool));
    bool* is_stored_B=(bool*)R_alloc(Ns_all, sizeof(bool));
    for (int checkdex=0; checkdex < Ns_all; ++checkdex) { 
        latest_pair_A[checkdex] = latest_pair_B[checkdex] = -1; 
        is_stored_A[checkdex] = is_stored_B[checkdex] = true;
    }

    int curpair=0, mode=0, maxmode, curq1=0, curq2=0, curindex=0, curs=0;
    int * latest_pair;
    bool * is_stored;
    for (curpair=0; curpair<Npairs; ++curpair) {
        maxmode = (a1ptr[curpair] == a2ptr[curpair] ? 1 : 2);

        /* Checking whether the first and second anchor overlaps anything in the opposing query sets.
         * Doing this twice; first and second anchors to the first and second query sets (A), then
         * the first and second anchors to the second and first query sets (B).
         */
        for (mode=0; mode<maxmode; ++mode) { 
            if (mode==0) { 
                curq1 = a1ptr[curpair];
                curq2 = a2ptr[curpair];
                if (curq1 >= Nq || curq1 < 0 || curq1==NA_INTEGER) { throw std::runtime_error("region index (1) out of bounds"); }
                if (curq2 >= Nq || curq2 < 0 || curq2==NA_INTEGER) { throw std::runtime_error("region index (2) out of bounds"); }
                latest_pair = latest_pair_A;
                is_stored = is_stored_A;
            } else {
                curq2 = a1ptr[curpair];
                curq1 = a2ptr[curpair];
                latest_pair = latest_pair_B;
                is_stored = is_stored_B;
            }

            for (curindex=qsptr1[curq1]; curindex<qeptr1[curq1]; ++curindex) {
                curs=sjptr1[curindex];
                if (mode && latest_pair_A[curs] == curpair && is_stored_A[curs]) { continue; } // Already added in first cycle.
                if (latest_pair[curs] < curpair) { 
                    latest_pair[curs] = curpair;
                    is_stored[curs] = false;
                }
            }

            for (curindex=qsptr2[curq2]; curindex<qeptr2[curq2]; ++curindex) {
                curs=sjptr2[curindex];
                if (mode && latest_pair_A[curs] == curpair && is_stored[curs]) { continue; }
                if (latest_pair[curs] == curpair && !is_stored[curs]) {
                    output->acknowledge(curpair, curs);
                    is_stored[curs] = true;
                }
            }
        }

        output->postprocess();
    }

    return;
}

/* Functions based on derived classes. */

SEXP expand_paired_olaps(SEXP anchor1, SEXP anchor2, SEXP querystarts1, SEXP queryends1, SEXP subject1, 
        SEXP querystarts2, SEXP queryends2, SEXP subject2, SEXP nsubjects) try {
    expanded_overlap x;
    detect_paired_olaps(&x, anchor1, anchor2, querystarts1, queryends1, subject1,
            querystarts2, queryends2, subject2, nsubjects);
    return x.generate();
} catch (std::exception& e) { 
    return mkString(e.what());
}

SEXP queryhit_paired_olaps(SEXP anchor1, SEXP anchor2, SEXP querystarts1, SEXP queryends1, SEXP subject1, 
        SEXP querystarts2, SEXP queryends2, SEXP subject2, SEXP nsubjects) try {
    query_overlap x;
    detect_paired_olaps(&x, anchor1, anchor2, querystarts1, queryends1, subject1,
            querystarts2, queryends2, subject2, nsubjects);
    return x.generate();
} catch (std::exception& e) { 
    return mkString(e.what());
}

SEXP subjecthit_paired_olaps(SEXP anchor1, SEXP anchor2, SEXP querystarts1, SEXP queryends1, SEXP subject1, 
        SEXP querystarts2, SEXP queryends2, SEXP subject2, SEXP nsubjects) try {
    subject_overlap x;
    detect_paired_olaps(&x, anchor1, anchor2, querystarts1, queryends1, subject1,
            querystarts2, queryends2, subject2, nsubjects);
    return x.generate();
} catch (std::exception& e) { 
    return mkString(e.what());
}

}
