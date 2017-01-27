CFLAGS = -O3 -funroll-all-loops

all: program.x


program.x: program.o
	gcc program.o -o program.x -lstdc++ -lm $(CFLAGS)

program.o: somik.cpp
	gcc -c somik.cpp -o program.o $(CFLAGS)

clean:
	rm -fr program.o program.x
