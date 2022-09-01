#include <stdlib.h>
#include "vectormap.h"
#include "memory.h"

void InitializeIndMap(Vectmap *v, const int max, const int dim)
{
    int i0, arrlen = 0;

    v->dim = dim;
    v->shift = max;
    v->ne = 2 * max;

    arrlen = v->ne;
    for (i0 = 1; i0 < dim; i0++)
    {
        arrlen *= v->ne;
    }
    v->maplen = arrlen;
    v->indmap = SafeCalloc(v->maplen, sizeof(int));
}

void DeallocIndMap(Vectmap *v)
{
    free(v->indmap);
}