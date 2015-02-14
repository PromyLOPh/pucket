#pragma once

#if defined(__clang__)
#if __has_extension(attribute_ext_vector_type)
/* LLVM/clang */
typedef double double2 __attribute__ ((ext_vector_type (2)));
typedef double double4 __attribute__ ((ext_vector_type (4)));
#endif
#else
/* GCC */
typedef double double2 __attribute__ ((vector_size (sizeof (double)*2)));
typedef double double4 __attribute__ ((vector_size (sizeof (double)*4)));
#endif

