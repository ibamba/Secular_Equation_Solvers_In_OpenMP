CC= gcc
CFLAGS= -g -O3 -Wall
EXEC= secular
OBJ= gragg.o hybride.o
LIB= -lm -fopenmp

all: $(OBJ)
	$(CC) $(CFLAGS) $(LIB) $(OBJ) -o $(EXEC)

CMD_EXEC= ./$(EXEC) 10000 1

gragg.o: gragg.c
	$(CC) -c $(CFLAGS) $< $(LIB)

hybride.o: hybride.c
%.o: %.c
	$(CC) -c $(CFLAGS) gragg.c hybride.c $(LIB)

.PHONY: clean

clean :
	rm -rf $(EXEC) *~ *.o
