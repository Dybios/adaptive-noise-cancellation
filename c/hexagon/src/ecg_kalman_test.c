#include <stdio.h>
#include <dlfcn.h>

//#include "test_main.h"

#include "HAP_perf.h"
#include "HAP_farf.h"

#include "ecg_kalman.h"

#define TEST_SUCCESS 0
#define TEST_FAILURE 1

#define TRY(exception, func) \
   if (TEST_SUCCESS != (exception = func)) {\
      goto exception##bail; \
   }

#define THROW(exception, errno) \
   exception = errno; \
   goto exception##bail;

#define CATCH(exception) exception##bail: if (exception != TEST_SUCCESS)


int test_main_start(int argc, char *argv[])
{
   long long int cyclesStart;
   long long int cyclesEnd;
   int nErr = TEST_SUCCESS;
//   FILE *input = NULL;
//   FILE *output = NULL;

   // Use defaults
//   input = fopen("/home/dybios/Qualcomm/Hexagon_SDK/3.5.4/examples/common/ecg_kalman/data/data_synthesized_0dB_64k.csv", "r");
//   output = fopen("/home/dybios/Qualcomm/Hexagon_SDK/3.5.4/examples/common/ecg_kalman/data/output/clean_output.csv", "w");

   if (argc < 2) {
      argv[1] = "/home/dybios/Qualcomm/Hexagon_SDK/3.5.4/examples/common/ecg_kalman/data/data_synthesized_0dB_64k.csv";
      argv[2] = "/home/dybios/Qualcomm/Hexagon_SDK/3.5.4/examples/common/ecg_kalman/data/output/clean_output.csv";
   }
   FARF(HIGH, "-- start lib test --                                                ");

   FARF(HIGH, "Calling template_so(%10d)                                       ", (int)&nErr);
   cyclesStart = HAP_perf_get_pcycles();
   nErr = ecg_kalman_main((int)&nErr, argv[1], argv[2]);
   cyclesEnd = HAP_perf_get_pcycles();

   FARF(ALWAYS, "Calling template_so() took %10d cycles                          ", (int)(cyclesEnd - cyclesStart));

   if (nErr == (int)&nErr) {
     nErr = TEST_SUCCESS;
   } else {
     FARF(ERROR, "template_so returned %10d instead of %10d", nErr, (int)&nErr);
     THROW(nErr, TEST_FAILURE);
   }


   FARF(HIGH, "Test Passed                                                         ");

   CATCH(nErr){};

   FARF(HIGH, "-- end lib test test --                                             ");

   return nErr;
}