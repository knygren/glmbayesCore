// @source_type: h
// @source_origin: Rinternals.h
// @includes: Rconfig.h, Boolean.h, Complex.h
// @depends: Rconfig
// @provides: R_INTERNALS_H_, Rbyte, R_len_t, R_xlen_t, SEXP, PROTECT_INDEX

#ifndef R_INTERNALS_H_
#define R_INTERNALS_H_

// openclport: include directives disabled for OpenCL C compilation.
// openclport: preload equivalent ported headers/shims in program assembly.
// openclport-disabled-include: #include <Rconfig.h>
// openclport-disabled-include: #include <R_ext/Boolean.h>
// openclport-disabled-include: #include <R_ext/Complex.h>

/*
 * Minimal R internals shim surface for parsing selected R_ext headers.
 * Opaque SEXP plus length/type aliases only.
 */
typedef unsigned char Rbyte;
typedef int R_len_t;

#if (SIZEOF_SIZE_T > 4)
typedef long long R_xlen_t;
#else
typedef int R_xlen_t;
#endif

typedef struct SEXPREC *SEXP;
typedef int PROTECT_INDEX;

#endif /* R_INTERNALS_H_ */
