CC = g++
CFLAGS = -Wall -std=c++11

OBJ = main.o hw7_1.o hw7_2.o hw7_3.o hw7_4.o hw7_5.o hw7_6.o equations.o

%.o: %.cpp hw7.h
	$(CC) $(CFLAGS) -c -o $@ $< 

hw7: $(OBJ)
	$(CC) $(CFLAGS) $^ -o $@

.PHONY: clean

clean:
	rm -f *.o