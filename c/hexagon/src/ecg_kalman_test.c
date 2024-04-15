#include <stdio.h>
#include <dlfcn.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>

#include "test_main.h"

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
   int ecg_complex_length = 0;
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
//   FILE *output_preprocess = NULL;
   FILE *output = NULL;

   // Use defaults
   input = fopen(argv[1], "r");
//   output_preprocess = fopen("/home/dybios/Qualcomm/Hexagon_SDK/3.5.4/examples/common/ecg_kalman/data/output/preprocessed.csv", "w");
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

   // Preprocess the ECG data to get the value of T
   FARF(HIGH, "ECG preprocessing starts...");
   double **preprocessed_data = (double **)malloc(rows * sizeof(double *));
   for (int i = 0; i < rows; i++) {
       preprocessed_data[i] = (double *)malloc(cols * sizeof(double));
   }
   nErr = preprocess_ecg_data(data, rows, preprocessed_data, &ecg_complex_length);
   FARF(HIGH, "ecg_complex_length (T) = %d", ecg_complex_length);
   FARF(HIGH, "ECG preprocessing done.");

   // Output for preprocessing stage
//   for (int i = 0; i < rows; i++) {
//       for (int j = 0; j < cols; j++) {
//           if (j != (cols - 1)) {
//               fprintf(output_preprocess, "%f%c", preprocessed_data[i][j], ',');
//           } else {
//               fprintf(output_preprocess, "%f", preprocessed_data[i][j]);
//           }
//       }
//       fprintf(output_preprocess, "%c", '\n'); // New line after every row
//   }

   FARF(HIGH, "Calling process_kalman(%10d)                                       ", (int)&nErr);
   int rem = rows % ecg_complex_length;
   int output_bufsize = rows + ecg_complex_length - rem;
   double *output_data = (double *)malloc(output_bufsize * sizeof(double));
   nErr = process_kalman(preprocessed_data, rows, cols, ecg_complex_length, output_data);
   FARF(ALWAYS, "process_kalman() complete.");

   if (nErr == (int)&nErr) {
     nErr = TEST_SUCCESS;
   } else {
     FARF(ERROR, "process_kalman returned %10d instead of %10d", nErr, (int)&nErr);
     THROW(nErr, TEST_FAILURE);
   }

   // Append the output buffer to output file
   for (int i = 0; i < rows; i++) {
       fprintf(output, "%f%c", output_data[i], '\n');
   }

   // Free the allocated memory
   for (int i = rows-1; i >= 0; i++) {
       free(preprocessed_data[i]);
       free(data[i]);
   }
   free(preprocessed_data);
   free(data);
   free(output_data);
   fclose(output);

   FARF(HIGH, "Test Passed                                                         ");

   CATCH(nErr){};

   FARF(HIGH, "-- end lib test test --                                             ");

   return nErr;
}
