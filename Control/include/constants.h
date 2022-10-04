#ifndef _CONSTANTS_H
#define _CONSTANTS_H
// FUNDAMENTAL
#define PI          3.14159265358979323846e+0       // -
#define HPLANCK     6.62607004e-34                  // [Js]
#define Q0          1.6021766208e-19                // [C]
#define M0          9.10938356e-31                  // [Kg]
#define KB          1.3806488e-23                   // [J/K]
#define EPS0        8.854187817620e-12              // [F/m]
#define C0          0.299792458e+9                  // [m/s]
#define SR2         1.41421356237309504880e+0       // -
#define SR3         1.73205080756887729352e+0       // -
#define ONEBYTHREE  3.33333333333333333333E-1       // -
#define TWOBYTHREE  6.66666666666666666666E-1       // -
#define RYD	    	13.6056981	                    //ev
#define BOHR		5.291772e-11	                //m
// DERIVED
#define TWOPI       (2.0e0*PI)                      // [1]
#define PIHALF      (PI/2.0e0)                      // [1]
#define HBAR        (HPLANCK/TWOPI)                 // [Js]
#define H2M0        ((HBAR*HBAR)/(2.0*M0))          // [Jm^2]
#define M0H2        ((2.0*M0)/(HBAR*HBAR))          // [J^-1m^-2]
#define H2M0Q0      ((HBAR*HBAR)/(2.0*M0*Q0))       // [eVm^2]
#define IH2M0Q0     ((2.0*Q0*M0)/(HBAR*HBAR))       // [eV^-1m^-2]
#define HM0         (HBAR/M0)                       // [m^2/s]
#define M0H         (M0/HBAR)                       // [s/m^2]
#define HBEV        (HBAR/Q0)                       // [eVs]
#define IHBEV       (Q0/HBAR)                       // [1/eVs]
#define EKQQ        ((EPS0*KB)/(Q0*Q0))             // [m]
#define KBEV        (KB/Q0)                         // [eV/K]
// TOLERANCE
#define TOL         1.0e-11
#define MINXRES     1.0e-3        // Minimum compositional resolution for material interpolation
//
// G-VECTOR RELATED
#define GMAX    9   // Max value for G-vec coordinate
#define GDIM    3   // Dimensionality of G-vectors
//
#if (__SIZEOF_DOUBLE__ == 8)
typedef union 
{
    double d;
    struct {
        unsigned long mantissa : 52;
        unsigned int exponent  : 11;
        unsigned int sign      : 1;
    } parts;
} double_cast;
#endif
#if (__SIZEOF_DOUBLE__ == 4)
typedef union 
{
    double d;
    struct {
        unsigned int mantissa : 23;
        unsigned int exponent : 8;
        unsigned int sign     : 1;
    } parts;
} float_cast;
#endif
//
#endif
