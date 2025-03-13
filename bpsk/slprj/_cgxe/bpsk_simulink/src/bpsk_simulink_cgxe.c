/* Include files */

#include "bpsk_simulink_cgxe.h"
#include "m_GvkxvAbVM3TqdJK85FbpIG.h"

unsigned int cgxe_bpsk_simulink_method_dispatcher(SimStruct* S, int_T method,
  void* data)
{
  if (ssGetChecksum0(S) == 4201392949 &&
      ssGetChecksum1(S) == 3411258648 &&
      ssGetChecksum2(S) == 1294117970 &&
      ssGetChecksum3(S) == 73053349) {
    method_dispatcher_GvkxvAbVM3TqdJK85FbpIG(S, method, data);
    return 1;
  }

  return 0;
}
