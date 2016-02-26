#include "iset.h"
#define REGISTER(x, i) {#x, (DL_FUNC) &x, i}

extern "C" {

static const R_CallMethodDef all_call_entries[] = {
    REGISTER(expand_olaps, 7),
    REGISTER(queryhit_olaps, 7),
    REGISTER(subjecthit_olaps, 7),
    REGISTER(expand_paired_olaps, 13),
    REGISTER(queryhit_paired_olaps, 13),
    REGISTER(subjecthit_paired_olaps, 13),
    REGISTER(expand_pair_links, 12),
    REGISTER(get_box_bounds, 6),
    {NULL, NULL, 0}
};

void attribute_visible R_init_InteractionSet(DllInfo *dll) {
    R_registerRoutines(dll, NULL, all_call_entries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
    R_forceSymbols(dll, TRUE);
}

}

