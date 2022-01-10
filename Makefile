MVBIN = -mv *.o gaus/*.o libge.a bin/
FLAGS = -Wall -pedantic -Wextra
CCC = cc -c
CCO = cc -o
DEL = -rm bin/*.o bin/libge.a aprox intrp prosta gen

all: aprox gen intrp prosta

aprox: main.o splines.o points.o aproksymator_na_bazie.o libge.a
	$(CCO) $@  $^ -L gaus

intrp: main.o splines.o points.o interpolator.o libge.a
	$(CCO) $@  $^ -L gaus

prosta: main.o splines.o points.o prosta.o
	$(CCO) $@  $^
	$(MVBIN)

aproksymator_na_bazie.o: aprx/aproksymator_na_bazie.c aprx/makespl.h aprx/points.h gaus/piv_ge_solver.h
	$(CCC) -I gaus $< $(FLAGS)

interpolator.o: aprx/interpolator.c aprx/makespl.h aprx/points.h gaus/piv_ge_solver.h
	$(CCC) -I gaus $< $(FLAGS)

main.o: aprx/main.c aprx/points.h aprx/splines.h aprx/makespl.h
	$(CCC) $< $(FLAGS)
	-mkdir bin

prosta.o: aprx/prosta.c aprx/makespl.h
	$(CCC) $< $(FLAGS)

splines.o: aprx/splines.c aprx/splines.h
	$(CCC) $< $(FLAGS)

points.o: aprx/points.c aprx/points.h
	$(CCC) $< $(FLAGS)

libge.a: gaus/matrix.o gaus/pivot.o gaus/piv_ge_solver.o
	ar rvs $@ $^

matrix.o: gaus/matrix.c gaus/matrix.h
	$(CCC) $< $(FLAGS)

pivot.o: gaus/pivot.c gaus/matrix.h
	$(CCC) $< $(FLAGS)

piv_ge_solver.o: gaus/piv_ge_solver.c gaus/piv_ge_solver.h gaus/matrix.h
	$(CCC) $< $(FLAGS)

gen: test/gen.o
	$(CCO) $@ $<
	-mv test/*.o bin/

gen.o: test/gen.c
	$(CCC) $< $(FLAGS)

.PHONY: clean

clean:
	$(DEL)
