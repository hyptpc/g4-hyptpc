// -*- C++ -*-

#ifndef GET_NUMBER_FROM_KERNEL_ENTROPY_POOL_HH
#define GET_NUMBER_FROM_KERNEL_ENTROPY_POOL_HH

#ifdef __cplusplus
extern "C" {
#endif

  int   GetIntFromKernelEntropyPool();
  short GetShortFromKernelEntropyPool();
  long  GetLongFromKernelEntropyPool();

#ifdef __cplusplus
}
#endif

#endif
