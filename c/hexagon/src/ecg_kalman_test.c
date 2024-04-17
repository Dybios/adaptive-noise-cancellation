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
   FARF(ALWAYS, "argv[1] = %s", argv[1]);
   FARF(ALWAYS, "argv[2] = %s", argv[2]);

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

   fseek(input, 0, SEEK_SET); // Move back to beginning of file

   // Allocate memory for the 2D array
   float **data = (float **)malloc(rows * sizeof(float *));
   for (int i = 0; i < rows; i++) {
       data[i] = (float *)malloc(cols * sizeof(float));
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

   FARF(ALWAYS, " data in memory                                              ");
   fclose(input);

   // Preprocess the ECG data to get the value of T
   FARF(ALWAYS, "ECG preprocessing starts...");
   float **preprocessed_data = (float **)malloc(rows * sizeof(float *));
   for (int i = 0; i < rows; i++) {
       preprocessed_data[i] = (float *)malloc(cols * sizeof(float));
   }
   nErr = preprocess_ecg_data(data, rows, preprocessed_data, &ecg_complex_length);
   FARF(ALWAYS, "ecg_complex_length (T) = %d", ecg_complex_length);
   FARF(ALWAYS, "ECG preprocessing done.");

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

   FARF(ALWAYS, "Calling process_kalman(%10d)                                       ", (int)&nErr);
   int rem = rows % ecg_complex_length;
   int output_bufsize = rows + ecg_complex_length - rem;
   float *output_data = (float *)malloc(output_bufsize * sizeof(float));
   nErr = process_kalman(preprocessed_data, rows, cols, ecg_complex_length, output_data);
   FARF(ALWAYS, "process_kalman() complete.");

   if (nErr == (int)&nErr) {
     nErr = TEST_SUCCESS;
   } else {
     FARF(ERROR, "process_kalman returned %10d instead of %10d", nErr, (int)&nErr);
     THROW(nErr, TEST_FAILURE);
   }

   // Append the output buffer to output file
//   for (int i = 0; i < rows; i++) {
//       fprintf(output, "%f%c", output_data[i], '\n');
//   }

   // Free the allocated memory
   for (int i = 0; i < rows; i++) {
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
