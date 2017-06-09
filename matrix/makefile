CC = g++
CFLAGS = -Wall -std=c++11

OBJ = main.o

%.o: %.cpp matrix.h
	$(CC) $(CFLAGS) -c -o $@ $< 

hw8: $(OBJ)
	$(CC) $(CFLAGS) $^ -o $@

.PHONY: clean

clean:
	rm -f *.o