#include <stdio.h>
#include <dlfcn.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>

#include "test_main.h"

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

int test_main_start(int argc, char *argv[]);

int main(int argc, char* argv[])
{
    return test_main_start(argc, argv);
}

int test_main_start(int argc, char *argv[])
{
   long long int cyclesStart;
   long long int cyclesEnd;
   int nErr = TEST_SUCCESS;

   if (argc < 2) {
      argv[1] = "/home/dybios/Qualcomm/Hexagon_SDK/3.5.4/examples/common/ecg_kalman/data/data_synthesized_0dB_64k.csv";
      argv[2] = "/home/dybios/Qualcomm/Hexagon_SDK/3.5.4/examples/common/ecg_kalman/data/output/clean_output.csv";
   }
   FARF(HIGH, "argv[1] = %s", argv[1]);
   FARF(HIGH, "argv[2] = %s", argv[2]);

   int rows = 0, cols = 0;
   char line[100];

   /** This section reads the files into the program memory.
    *  Input can be provided by either cmdline arguments or if run as is, defaults are taken.
    */
   FILE *input = NULL;
   FILE *output_preprocess = NULL;
   FILE *output = NULL;

   // Use defaults
   input = fopen(argv[1], "r");
   output_preprocess = fopen("/home/dybios/Qualcomm/Hexagon_SDK/3.5.4/examples/common/ecg_kalman/data/output/preprocessed.csv", "w");
   output = fopen(argv[2], "w");

   // Count rows and columns
   while (fgets(line, sizeof(line), input) != NULL) {
       cols = 0;
       char *token = strtok(line, ",");
       while (token != NULL) {
           cols++;
           token = strtok(NULL, ",");
       }
       rows++; // data in each channel length
   }

   FARF(HIGH, " fopen done                                               ");
   fseek(input, 0, SEEK_SET); // Move back to beginning of file
   
   // Allocate memory for the 2D array
   double **data = (double **)malloc(rows * sizeof(double *));
   for (int i = 0; i < rows; i++) {
       data[i] = (double *)malloc(cols * sizeof(double));
   }

   double *output_data = (double *)malloc(rows * sizeof(double));

   FARF(HIGH, " data and output_data allocated                                               ");

   // Read the CSV data into the array
   int i = 0, j = 0;
   while (fgets(line, sizeof(line), input) != NULL) {
       j = 0;
       char *token = strtok(line, ",");
       while (token != NULL) {
           data[i][j++] = atof(token);
           token = strtok(NULL, ",");
       }
       i++;
   }
   
   FARF(HIGH, " put into the array                                              ");
   fclose(input); 

   FARF(HIGH, "-- start lib test --                                                ");

   FARF(HIGH, "Calling ecg_kalman(%10d)                                       ", (int)&nErr);
   cyclesStart = HAP_perf_get_pcycles();
   nErr = ecg_kalman_main((int)&nErr, data, rows, cols, output_data);
   cyclesEnd = HAP_perf_get_pcycles();

   FARF(ALWAYS, "Calling ecg_kalman() took %10d cycles                          ", (int)(cyclesEnd - cyclesStart));

   if (nErr == (int)&nErr) {
     nErr = TEST_SUCCESS;
   } else {
     FARF(ERROR, "ecg_kalman returned %10d instead of %10d", nErr, (int)&nErr);
     THROW(nErr, TEST_FAILURE);
   }
   
   // Append the output buffer to output file
   for (int i = 0; i < rows; i++) {
       fprintf(output, "%f%c", output_data[i], '\n');
   }

   // Free the allocated memory
   for (int i = 0; i < rows; i++) {
       free(data[i]);
   }
   free(data);
   free(output);

   FARF(HIGH, "Test Passed                                                         ");

   CATCH(nErr){};

   FARF(HIGH, "-- end lib test test --                                             ");

   return nErr;
}
