clean:
	rm -f *.o
	rm -f *.out
	rm -f *.gf

run:
	g++ -std=c++11 main.cpp -o main.out \
	-I"../Libraries/mfem-4.5/mfem/include" \
	-L../Libraries/mfem-4.5/mfem/lib -lmfem -lrt
	@chmod a+x main.out
	@./main.out

vis-sol:
	/bin/sh -c '../Libraries/glvis-4.2/glvis -m mesh.mesh -g sol_eps.gf sol_hom.gf'

vis-mesh:
	../Libraries/glvis-4.2/glvis -m mesh.mesh




