#ifndef _VECTORMAP_H
#define _VECTORMAP_H
//
// This struct is a map between a vector in an array and it's index in that array
typedef struct vectmap
{
    int dim;     // Dimension
    int ne;      // # elements per dim
    int shift;   // Center index of each dim
    int *indmap; // Vector <-> index map
    int maplen;
} Vectmap;
// Shift coord so >0, then store x, y, z in
// first, second, and third digit of number, base ne.
#define VECMAPIND3(v, x, y, z) (v->indmap[(x + v->shift) +           \
                                          ((y + v->shift) * v->ne) + \
                                          ((z + v->shift) * v->ne * v->ne)])
#define VECMAPGETIND3(v, x, y, z) (x + v->shift) +               \
                                      ((y + v->shift) * v->ne) + \
                                      ((z + v->shift) * v->ne * v->ne)

// Properly initializes Vectmap, and allocates spaces for index mapping vector.
void InitializeIndMap(Vectmap *v, const int max, const int dim);
void DeallocIndMap(Vectmap *v);

#endif
