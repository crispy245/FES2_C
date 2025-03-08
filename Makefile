N = 5
M = 5

run: fes2
	./fes2

malloc: fes2_malloc
	./fes2_malloc

build: fes2 fes2_malloc

fes2: fes2.c
	gcc -DN=$(N) -DM=$(M) fes2.c -o fes2

fes2_malloc: fes2_malloc.c
	gcc -DN=$(N) -DM=$(M) fes2_malloc.c -o fes2_malloc

clean:
	rm -f fes2 fes2_malloc