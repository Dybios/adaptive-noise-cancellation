CFLAGS = -Wall -Werror -O2 -g #-fsanitize=address
LDLIBS = -lm

ecg_kalman_test: ecg_kalman_test.o ecg_kalman.o matrix.o filter.o qrs_detection.o qrs_filter.o
	gcc -o ecg_kalman_test ecg_kalman_test.o ecg_kalman.o matrix.o filter-c/filter.o qrs_detector/c89/qrs_detection.o qrs_detector/c89/filter.o $(LDLIBS) $(CFLAGS)

ecg_kalman.o: ecg_kalman.c ecg_kalman.h
	gcc -o ecg_kalman.o -c ecg_kalman.c $(LDLIBS) $(CFLAGS)

matrix.o: matrix.c matrix.h
	gcc -o matrix.o -c matrix.c $(LDLIBS) $(CFLAGS)

filter.o: filter-c/filter.c filter-c/filter.h
	gcc -o filter-c/filter.o -c filter-c/filter.c $(LDLIBS) $(CFLAGS)

qrs_detection.o: qrs_detector/c89/qrs_detection.c qrs_detector/c89/qrs_detection.h
	gcc -o qrs_detector/c89/qrs_detection.o -c qrs_detector/c89/qrs_detection.c  $(LDLIBS) $(CFLAGS)

qrs_filter.o: qrs_detector/c89/filter.c qrs_detector/c89/filter.h
	gcc -o qrs_detector/c89/filter.o -c qrs_detector/c89/filter.c $(LDLIBS) $(CFLAGS)

clean:
	$(RM) ecg_kalman_test
	$(RM) *.o filter-c/filter.o qrs_detector/c89/qrs_detection.o qrs_detector/c89/qrs_filter.o
	$(RM) *.gch filter-c/*.gch qrs_detector/c89/*.gch
