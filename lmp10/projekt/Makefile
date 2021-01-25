aprox: main.o splines.o points.o aproksymator_na_bazie.o gaus/libge.a
	$(CC) -o aprox  main.o splines.o points.o aproksymator_na_bazie.o -L gaus -l ge

intrp: main.o splines.o points.o interpolator.o gaus/libge.a
	$(CC) -o intrp  main.o splines.o points.o interpolator.o -L gaus -l ge

prosta: main.o splines.o points.o prosta.o
	$(CC) -o prosta  main.o splines.o points.o prosta.o	

aproksymator_na_bazie.o: makespl.h points.h gaus/piv_ge_solver.h
	$(CC) -I gaus -c aproksymator_na_bazie.c

interpolator.o: makespl.h points.h gaus/piv_ge_solver.h
	$(CC) -I gaus -c interpolator.c

aproksymacja_czebyszewa.o: programy/makespl.h programy/points.h gaus/piv_ge_solver.h
	$(CC) -I gaus -I programy -c programy/aproksymacja_czebyszewa.c

main.o: programy/main.c programy/points.h programy/splines.h programy/makespl.h
	$(CC) -I programy -c programy/main.c 

splines.o: programy/splines.c programy/splines.h
	$(CC) -I programy -c programy/splines.c

points.o: programy/points.c programy/points.h
	$(CC) -I programy -c programy/points.c

ap: main.o splines.o points.o aproksymacja_czebyszewa.o gaus/libge.a
	$(CC) -o aprox  main.o splines.o points.o aproksymacja_czebyszewa.o -lm -L gaus -l ge 

.PHONY: clean

clean:
	-rm *.o aprox
