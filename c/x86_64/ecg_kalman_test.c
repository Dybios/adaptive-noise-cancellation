#include "ecg_kalman.h"

int main(int argc, char *argv[])
{
   int nErr = 0;
   int ecg_complex_length = 0;
   int rows = 0, cols = 0;
   char line[100];

   /** This section reads the files into the program memory.
    *  Input can be provided by either cmdline arguments or if run as is, defaults are taken.
    */
    FILE *input = NULL;
    FILE *output_preprocess = NULL;
    FILE *output = NULL;
    if (argc == 1) {
       // Use defaults
       input = fopen("../../data/data_synthesized_8dB.csv", "r");
       output_preprocess = fopen("../../data/output/preprocessed.csv", "w");
       output = fopen("../../data/output/clean_output.csv", "w");
    }
    else if (argc > 1 && argc < 3) {
       printf("Not enough arguments.\n");
       exit(1);
    }
    else if (argc > 3) {
       printf("Too many arguments.\n");
       exit(1);
    }
    else {
       input = fopen(argv[1], "r");
       output = fopen(argv[2], "w");
       output_preprocess = fopen("../../data/output/preprocessed.csv", "w");
    }

    if (input == NULL || output == NULL) {
           printf("Cannot open file.\n");
    }

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
   fclose(input);

   // Preprocess the ECG data to get the value of T
   printf("ECG preprocessing starts...\n");
   double **preprocessed_data = (double **)malloc(rows * sizeof(double *));
   for (int i = 0; i < rows; i++) {
       preprocessed_data[i] = (double *)malloc(cols * sizeof(double));
   }
   nErr = preprocess_ecg_data(data, rows, preprocessed_data, &ecg_complex_length);
   printf("ecg_complex_length (T) = %d\n", ecg_complex_length);
   printf("ECG preprocessing done. Writing to file...\n");

   // Output for preprocessing stage
   for (int i = 0; i < rows; i++) {
       for (int j = 0; j < cols; j++) {
           if (j != (cols - 1)) {
               fprintf(output_preprocess, "%f%c", preprocessed_data[i][j], ',');
           } else {
               fprintf(output_preprocess, "%f", preprocessed_data[i][j]);
           }

       }
       fprintf(output_preprocess, "%c", '\n'); // New line after every row
   }

   printf("process_kalman main process starts...\n");
   int rem = (int)rows % ecg_complex_length;
   int output_bufsize = rows + ecg_complex_length - rem;
   double *output_data = (double *)malloc(output_bufsize * sizeof(double));
   nErr = process_kalman(preprocessed_data, rows, cols, ecg_complex_length, output_data);
   printf("process_kalman completed. Writing to file...\n");

   // Append the output buffer to output file
   for (int i = 0; i < output_bufsize; i++) {
       fprintf(output, "%f%c", output_data[i], '\n');
   }

   // Free the allocated memory
   printf("freeing memory...\n");
   for (int i = rows-1; i >= 0; i++) {
       free(data[i]);
   }
   free(data);
   free(output_data);
   fclose(output);

   return nErr;
}
