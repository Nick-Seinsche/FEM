clean:
	rm -f *.o
	rm -f *.out

run:
	@chmod a+x main.out
	@./main.out

compile:
	g++ -I../Libraries/mfem-4.5/mfem/include -L../Libraries/mfem-4.5/mfem/lib -lmfem main.cpp -o main.out


