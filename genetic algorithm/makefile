CC = g++
CFLAGS = -std=c++11

OBJ = main.o

%.o: %.cpp GA.h selection.h recombination.h mutation.h
	$(CC) $(CFLAGS) -c -o $@ $< 

hw9: $(OBJ)
	$(CC) $(CFLAGS) $^ -o $@

.PHONY: clean

clean:
	rm -f *.o