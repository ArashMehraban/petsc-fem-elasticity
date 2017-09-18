#ifndef _fe_align_h
#define _fe_align_h

#if defined(__GNUC__)
#  define _align __attribute__((aligned(16))) /* SSE instructions need 16-byte alignment */
#else
#  define _align
#endif

#endif
