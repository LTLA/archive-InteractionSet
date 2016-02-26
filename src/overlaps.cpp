#include "iset.h"

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

void set_mode_values (SEXP use_both, int& start, int& end) {
    if (!isInteger(use_both) || LENGTH(use_both)!=1) { throw std::runtime_error("'use_both' specifier should be an integer scalar"); }
    switch (asInteger(use_both)) {
        case 1:
            start=0; end=2; break;
        case 2:
            start=0; end=1; break;
        case 3:
            start=1; end=2; break;
        default:
            throw std::runtime_error("invalid specification for 'use_both'");
    }
    return;
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
    void prime(int nq, int ns) {
        (void)nq; (void)ns;
    };
    void acknowledge(int q, int s) {
        new_query.push_back(q);
        new_subject.push_back(s);
        ++just_added;
    }
    void postprocess() { // Sorting for output.
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
        (void)ns;
        nquery=nq;
        query_hit.resize(nq);    
    };
    void acknowledge(int q, int s) {
        (void)s;
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
        (void)nq;
        nsubject=ns;
        subject_hit.resize(ns);    
    };
    void acknowledge(int q, int s) {
        (void)q;
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

/* Base function, to detect all overlaps between linear ranges and either interacting loci in a pair */

void detect_olaps(output_store* output, SEXP anchor1, SEXP anchor2, SEXP querystarts, SEXP queryends, SEXP subject, SEXP nsubjects, SEXP use_both) {
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

    int true_mode_start, true_mode_end;
    set_mode_values(use_both, true_mode_start, true_mode_end);

    // Checking indices. 
    check_indices(qsptr, qeptr, Nq, sjptr, Ns, Ns_all);

    /* Constructing an output deque */
    output->prime(Npairs, Ns_all);
    int* latest_pair=(int*)R_alloc(Ns_all, sizeof(int));
    for (int checkdex=0; checkdex < Ns_all; ++checkdex) { latest_pair[checkdex] = -1; }

    int curpair=0, mode=0, maxmode, curq=0, curindex=0, curs=0;
    for (curpair=0; curpair<Npairs; ++curpair) {
        maxmode=(a1ptr[curpair]==a2ptr[curpair] ? true_mode_start+1 : true_mode_end); // just check one or the other, if they're the same.

        for (mode=true_mode_start; mode<maxmode; ++mode) { 
            if (mode == 0) {
                curq = a1ptr[curpair];
            } else {
                curq = a2ptr[curpair];
            }

            if (curq >= Nq || curq < 0 || curq==NA_INTEGER) { throw std::runtime_error("region index out of bounds"); }
            for (curindex=qsptr[curq]; curindex<qeptr[curq]; ++curindex) {
                curs=sjptr[curindex];
                if (latest_pair[curs] < curpair) { 
                    output->acknowledge(curpair, curs);
                    latest_pair[curs] = curpair;
                }
            }
        }

        output->postprocess();
    }

    return;
}

/* Base function, to detect all 2D overlaps between two GInteractions objects.
 * This is pretty complicated:
 *
 *   anchor1, anchor2 hold the anchor indices of the query GInteractions.
 *   querystarts, queryend hold the first and one-after-last row of the Hits object obtained by findOverlaps(regions(query), regions(subject))
 *   subject holds the queryHits of the aforementioned Hits object.
 *   next_anchor_start1, next_anchor_end1 holds the first and one-after-last position of the sorted subject@anchor1.
 *   next_id1 holds the order of subject@anchor1.
 *   next_anchor_start2, next_anchor_end2 holds the first and one-after-last position of the sorted subject@anchor2.
 *   next_id2 holds the order of subject@anchor2.
 *   next_num_pairs holds the total number of pairs in the subject.
 * 
 * The idea is to:
 *   1) go through each query pair;
 *   2) retrieve the regions(subject) overlapping each anchor query;
 *   3) identify all subject anchors matching each retrieved regions(subject)
 *   4) cross-reference all identified subject anchors 1 and 2 to identify 2D overlaps.
 */

void detect_paired_olaps(output_store* output, SEXP anchor1, SEXP anchor2, 
        SEXP querystarts, SEXP queryends, SEXP subject, 
        SEXP next_anchor_start1, SEXP next_anchor_end1, SEXP next_id1,
        SEXP next_anchor_start2, SEXP next_anchor_end2, SEXP next_id2,
        SEXP num_next_pairs, SEXP use_both) {

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

    if (!isInteger(next_anchor_start1) || !isInteger(next_anchor_end1)) { throw std::runtime_error("next indices (1) must be integer vectors"); }
    const int Nas=LENGTH(next_anchor_start1);
    if (Nas != LENGTH(next_anchor_end1)) { throw std::runtime_error("next indices (1) must be of equal length"); }
    const int* nasptr1=INTEGER(next_anchor_start1), *naeptr1=INTEGER(next_anchor_end1);  
    if (!isInteger(next_id1)) { throw std::runtime_error("next ID indices (1) must be integer"); }
    const int *niptr1=INTEGER(next_id1);

    if (!isInteger(next_anchor_start2) || !isInteger(next_anchor_end2)) { throw std::runtime_error("next indices (2) must be integer vectors"); }
    if (Nas != LENGTH(next_anchor_start2) || Nas != LENGTH(next_anchor_end2)) { throw std::runtime_error("next indices (2) must be of equal length"); }
    const int* nasptr2=INTEGER(next_anchor_start2), *naeptr2=INTEGER(next_anchor_end2);  
    if (!isInteger(next_id2)) { throw std::runtime_error("next ID indices (2) must be integer"); }
    const int *niptr2=INTEGER(next_id2);

    if (!isInteger(num_next_pairs) || LENGTH(num_next_pairs)!=1) { throw std::runtime_error("total number of next pairs must be an integer scalar"); }
    const int Nnp = asInteger(num_next_pairs);
    if (LENGTH(next_id1)!=Nnp || LENGTH(next_id2)!=Nnp) { throw std::runtime_error("number of next IDs is not equal to specified number of pairs"); }

    int true_mode_start, true_mode_end;
    set_mode_values(use_both, true_mode_start, true_mode_end);
 
    // Check indices.
    check_indices(qsptr, qeptr, Nq, sjptr, Ns, Nas);
    check_indices(nasptr1, naeptr1, Nas, niptr1, Nnp, Nnp);
    check_indices(nasptr2, naeptr2, Nas, niptr2, Nnp, Nnp);

    // Setting up logging arrays. 
    output->prime(Npairs, Nnp);
    int* latest_pair_A=(int*)R_alloc(Nnp, sizeof(int));
    int* latest_pair_B=(int*)R_alloc(Nnp, sizeof(int));
    bool* is_stored_A=(bool*)R_alloc(Nnp, sizeof(bool));
    bool* is_stored_B=(bool*)R_alloc(Nnp, sizeof(bool));
    for (int checkdex=0; checkdex < Nnp; ++checkdex) { 
        latest_pair_A[checkdex] = latest_pair_B[checkdex] = -1; 
        is_stored_A[checkdex] = is_stored_B[checkdex] = true;
    }

    int curpair=0, mode=0, maxmode, curq1=0, curq2=0, 
        curindex=0, cur_subreg=0, cur_nextanch=0, cur_nextid;
    int * latest_pair;
    bool * is_stored;
    for (curpair=0; curpair<Npairs; ++curpair) {
        maxmode = (a1ptr[curpair] == a2ptr[curpair] ? true_mode_start+1 : true_mode_end);

        /* Checking whether the first and second anchor overlaps anything in the opposing query sets.
         * Doing this twice; first and second anchors to the first and second query sets (A), then
         * the first and second anchors to the second and first query sets (B).
         */
        for (mode=true_mode_start; mode<maxmode; ++mode) { 
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

            for (curindex=qsptr[curq1]; curindex<qeptr[curq1]; ++curindex) {
                cur_subreg=sjptr[curindex];
                for (cur_nextanch=nasptr1[cur_subreg]; cur_nextanch<naeptr1[cur_subreg]; ++cur_nextanch) {
                    cur_nextid=niptr1[cur_nextanch];
                    if (mode && latest_pair_A[cur_nextid] == curpair && is_stored_A[cur_nextid]) { continue; } // Already added in first cycle.
                    if (latest_pair[cur_nextid] < curpair) { 
                        latest_pair[cur_nextid] = curpair;
                        is_stored[cur_nextid] = false;
                    }
                }
            }

            for (curindex=qsptr[curq2]; curindex<qeptr[curq2]; ++curindex) {
                cur_subreg=sjptr[curindex];
                for (cur_nextanch=nasptr2[cur_subreg]; cur_nextanch<naeptr2[cur_subreg]; ++cur_nextanch) {
                    cur_nextid=niptr2[cur_nextanch];
                    if (mode && latest_pair_A[cur_nextid] == curpair && is_stored_A[cur_nextid]) { continue; }
                    if (latest_pair[cur_nextid] == curpair && !is_stored[cur_nextid]) {
                        output->acknowledge(curpair, cur_nextid);
                        is_stored[cur_nextid] = true;
                    }
                }
            }
        }

        output->postprocess();
    }

    return;
}

}

/* Functions based on derived classes. */

SEXP expand_olaps(SEXP anchor1, SEXP anchor2, SEXP querystarts, SEXP queryends, SEXP subject, SEXP nsubjects, SEXP use_both) try {
    expanded_overlap x;
    detect_olaps(&x, anchor1, anchor2, querystarts, queryends, subject, nsubjects, use_both);
    return x.generate();
} catch (std::exception& e) {
    return mkString(e.what());
}

SEXP queryhit_olaps(SEXP anchor1, SEXP anchor2, SEXP querystarts, SEXP queryends, SEXP subject, SEXP nsubjects, SEXP use_both) try {
    query_overlap x;
    detect_olaps(&x, anchor1, anchor2, querystarts, queryends, subject, nsubjects, use_both);
    return x.generate();
} catch (std::exception& e) {
    return mkString(e.what());
}

SEXP subjecthit_olaps(SEXP anchor1, SEXP anchor2, SEXP querystarts, SEXP queryends, SEXP subject, SEXP nsubjects, SEXP use_both) try {
    subject_overlap x;
    detect_olaps(&x, anchor1, anchor2, querystarts, queryends, subject, nsubjects, use_both);
    return x.generate();
} catch (std::exception& e) {
    return mkString(e.what());
}

SEXP expand_paired_olaps(SEXP anchor1, SEXP anchor2, 
        SEXP querystarts, SEXP queryends, SEXP subject,
        SEXP next_anchor_start1, SEXP next_anchor_end1, SEXP next_id1,
        SEXP next_anchor_start2, SEXP next_anchor_end2, SEXP next_id2,
        SEXP num_next_pairs, SEXP use_both) try {
    expanded_overlap x;
    detect_paired_olaps(&x, anchor1, anchor2, querystarts, queryends, subject,
            next_anchor_start1, next_anchor_end1, next_id1,
            next_anchor_start2, next_anchor_end2, next_id2,
            num_next_pairs, use_both);
    return x.generate();
} catch (std::exception& e) { 
    return mkString(e.what());
}

SEXP queryhit_paired_olaps(SEXP anchor1, SEXP anchor2, 
        SEXP querystarts, SEXP queryends, SEXP subject,
        SEXP next_anchor_start1, SEXP next_anchor_end1, SEXP next_id1,
        SEXP next_anchor_start2, SEXP next_anchor_end2, SEXP next_id2,
        SEXP num_next_pairs, SEXP use_both) try {
    query_overlap x;
    detect_paired_olaps(&x, anchor1, anchor2, querystarts, queryends, subject,
            next_anchor_start1, next_anchor_end1, next_id1,
            next_anchor_start2, next_anchor_end2, next_id2,
            num_next_pairs, use_both);
    return x.generate();
} catch (std::exception& e) { 
    return mkString(e.what());
}

SEXP subjecthit_paired_olaps(SEXP anchor1, SEXP anchor2, 
        SEXP querystarts, SEXP queryends, SEXP subject,
        SEXP next_anchor_start1, SEXP next_anchor_end1, SEXP next_id1,
        SEXP next_anchor_start2, SEXP next_anchor_end2, SEXP next_id2,
        SEXP num_next_pairs, SEXP use_both) try {
    subject_overlap x;
    detect_paired_olaps(&x, anchor1, anchor2, querystarts, queryends, subject,
            next_anchor_start1, next_anchor_end1, next_id1,
            next_anchor_start2, next_anchor_end2, next_id2,
            num_next_pairs, use_both);
    return x.generate();
} catch (std::exception& e) { 
    return mkString(e.what());
}

/* Another function that links two regions together. This does something quite similar
 * to detect_paired_olaps, but whereas that function requires both anchor regions to
 * overlap the regions from the same subject pair, here we only require that the
 * anchor regions overlap one region from either set. We then store all possible 
 * combinations of regions that are mediated by our interactions.
 */

SEXP expand_pair_links(SEXP anchor1, SEXP anchor2, SEXP querystarts1, SEXP queryends1, SEXP subject1, SEXP nsubjects1, 
        SEXP querystarts2, SEXP queryends2, SEXP subject2, SEXP nsubjects2, SEXP sameness, SEXP use_both) try {

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

    if (!isInteger(nsubjects1) || LENGTH(nsubjects1)!=1) { throw std::runtime_error("total number of subjects (1) must be an integer scalar"); }
    const int Ns_all1 = asInteger(nsubjects1);
    if (!isInteger(nsubjects2) || LENGTH(nsubjects2)!=1) { throw std::runtime_error("total number of subjects (2) must be an integer scalar"); }
    const int Ns_all2 = asInteger(nsubjects2);

    int true_mode_start, true_mode_end;
    set_mode_values(use_both, true_mode_start, true_mode_end);           
    if (!isLogical(sameness) || LENGTH(sameness)!=1) { throw std::runtime_error("same region indicator should be a logical scalar"); }
    const bool is_same=asLogical(sameness);
    if (is_same) { true_mode_end=true_mode_start+1; } // No point examining the flipped ones.
   
    // Check indices.
    check_indices(qsptr1, qeptr1, Nq, sjptr1, Ns1, Ns_all1);
    check_indices(qsptr2, qeptr2, Nq, sjptr2, Ns2, Ns_all2);

    // Setting up the set. 
    typedef std::pair<int, int> link;
    std::set<link> currently_active;
    std::set<link>::const_iterator itca;
    std::deque<link> stored_inactive;
    std::deque<int> interactions;

    int curpair=0, mode=0, maxmode, curq1=0, curq2=0, curindex1=0, curindex2=0;
    for (curpair=0; curpair<Npairs; ++curpair) {
        maxmode = (a1ptr[curpair] == a2ptr[curpair] ? true_mode_start+1 : true_mode_end);

        // Repeating with switched anchors, if the query sets are not the same.
        for (mode=true_mode_start; mode<maxmode; ++mode) { 
            if (mode==0) { 
                curq1 = a1ptr[curpair];
                curq2 = a2ptr[curpair];
                if (curq1 >= Nq || curq1 < 0 || curq1==NA_INTEGER) { throw std::runtime_error("region index (1) out of bounds"); }
                if (curq2 >= Nq || curq2 < 0 || curq2==NA_INTEGER) { throw std::runtime_error("region index (2) out of bounds"); }
            } else {
                curq2 = a1ptr[curpair];
                curq1 = a2ptr[curpair];
            }

            /* Storing all combinations associated with this pair. Avoiding redundant combinations
             * for self-linking (we can't use a simple rule like curindex2 < curindex1 to restrict 
             * the loop, as the second anchor can overlap regions above the first anchor if the 
             * latter is nested within the former).
             */
            for (curindex1=qsptr1[curq1]; curindex1<qeptr1[curq1]; ++curindex1) {
                for (curindex2=qsptr2[curq2]; curindex2<qeptr2[curq2]; ++curindex2) {
                    if (is_same && sjptr1[curindex1] < sjptr2[curindex2]) {
                        currently_active.insert(link(sjptr2[curindex2], sjptr1[curindex1]));
                    } else {
                        currently_active.insert(link(sjptr1[curindex1], sjptr2[curindex2]));
                    }
                }
            }
        }

        // Relieving the set by storing all inactive entries (automatically sorted as well).
        for (itca=currently_active.begin(); itca!=currently_active.end(); ++itca) {
            stored_inactive.push_back(*itca);
        }
        interactions.resize(interactions.size() + currently_active.size(), curpair);
        currently_active.clear();
    }

    // Popping back a list of information.
    SEXP output=PROTECT(allocVector(VECSXP, 3));
    try {
        const int total_entries=interactions.size();
        SET_VECTOR_ELT(output, 0, allocVector(INTSXP, total_entries));
        int * oiptr=INTEGER(VECTOR_ELT(output, 0));
        SET_VECTOR_ELT(output, 1, allocVector(INTSXP, total_entries));
        int * os1ptr=INTEGER(VECTOR_ELT(output, 1));
        SET_VECTOR_ELT(output, 2, allocVector(INTSXP, total_entries));
        int * os2ptr=INTEGER(VECTOR_ELT(output, 2));
    
        std::copy(interactions.begin(), interactions.end(), oiptr);
        for (int curdex=0; curdex < total_entries; ++curdex) {
            os1ptr[curdex]=stored_inactive[curdex].first;
            os2ptr[curdex]=stored_inactive[curdex].second;
        }
    } catch (std::exception& e) {
        UNPROTECT(1);
        throw;
    }

    UNPROTECT(1);
    return output;
} catch (std::exception& e) {
    return mkString(e.what());
}

