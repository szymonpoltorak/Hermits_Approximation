MVBIN = -mv *.o gaus/*.o libge.a bin/
FLAGS = -Wall -pedantic -Wextra
CCC = cc -c
CCO = cc -o
DEL = -rm bin/*.o bin/libge.a aprox intrp prosta gen hermit

all: aprox gen intrp hermit prosta

aprox: main.o splines.o points.o aproksymator_na_bazie.o libge.a
	$(CCO) $@  $^ -L gaus

intrp: main.o splines.o points.o interpolator.o libge.a
	$(CCO) $@  $^ -L gaus

hermit: main.o splines.o points.o aprx_hermit.o libge.a
	$(CCO) $@ $^ -L gaus

test_hermit: hermit
	echo
	echo Testuje uzycie punktow z aproksymacja
	echo
	./hermit -s spl -p test/dane.1 -g myplot -f 5.1 -t 5.7 -n 300
	gnuplot -p -e "plot 'test/dane.1', 'myplot'"
	echo
	echo Testuje uzycie funkcji sklejanych bez uzycia punktow
	echo
	./hermit -s spl  -g innyplot -f 5.6 -t 5.997 -n 300
	gnuplot -p -e "plot 'myplot', 'innyplot'"
	echo
	echo Testuje uzycie funkcji sklejach z uzyciem punktow
	echo
	./hermit -s spl -p test/dane.1 -g 3spl
	gnuplot -p -e "plot 'test/dane.1', '3spl'"
	$(MVBIN)
	-rm myplot innyplot 3spl

test_aprox: aprox
	echo
	echo Testuje uzycie punktow z aproksymacja
	echo
	./aprox -s spl -p test/dane.1 -g myplot -f 5.1 -t 5.7 -n 300
	gnuplot -p -e "plot 'test/dane.1', 'myplot'"
	echo
	echo Testuje uzycie funkcji sklejanych bez uzycia punktow
	echo
	./aprox -s spl  -g innyplot -f 5.6 -t 5.997 -n 300
	gnuplot -p -e "plot 'myplot', 'innyplot'"
	echo
	echo Testuje uzycie funkcji sklejach z uzyciem punktow
	echo
	./aprox -s spl -p test/dane.1 -g 3spl
	gnuplot -p -e "plot 'test/dane.1', '3spl'"
	$(MVBIN)
	-rm myplot innyplot 3spl

test_intrp: intrp
	echo
	echo Testuje uzycie punktow z aproksymacja
	echo
	./intrp -s spl -p test/dane.1 -g myplot -f 5.1 -t 5.7 -n 300
	gnuplot -p -e "plot 'test/dane.1', 'myplot'"
	echo
	echo Testuje uzycie funkcji sklejanych bez uzycia punktow
	echo
	./intrp -s spl  -g innyplot -f 5.6 -t 5.997 -n 300
	gnuplot -p -e "plot 'myplot', 'innyplot'"
	echo
	echo Testuje uzycie funkcji sklejach z uzyciem punktow
	echo
	./intrp -s spl -p test/dane.1 -g 3spl
	gnuplot -p -e "plot 'test/dane.1', '3spl'"
	$(MVBIN)
	-rm myplot innyplot 3spl

prosta: main.o splines.o points.o prosta.o
	$(CCO) $@  $^
	$(MVBIN)

aprx_hermit.o: aprx/aprx_hermit.c aprx/makespl.h gaus/piv_ge_solver.h aprx/points.h
	$(CCC) -I gaus $< $(FLAGS)

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
