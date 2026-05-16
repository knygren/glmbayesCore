// @source_type: h
// @source_origin: Rdefines.h
// @includes: Rinternals.h
// @depends: Rinternals
// @provides: R_DEFINES_H_, NA_STRING

#ifndef R_DEFINES_H_
#define R_DEFINES_H_

// openclport: include directives disabled for OpenCL C compilation.
// openclport: preload equivalent ported headers/shims in program assembly.
// openclport-disabled-include: #include <Rinternals.h>

/* Minimal compatibility aliases used by some downstream headers. */
#ifndef NA_STRING
#define NA_STRING ((SEXP)0)
#endif

#endif /* R_DEFINES_H_ */
