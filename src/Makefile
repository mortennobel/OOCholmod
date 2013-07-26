include make.inc

lib:
	rm -rf *.o liboochol.a
	$(CXX) -c $(INC) $(FLAGS) config_singleton.cpp dense_matrix.cpp factor.cpp sparse_matrix.cpp oo_blas.cpp oo_lapack.cpp 
	ar cr liboochol.a *.o
	rm -rf *.o

nullptr:
	rm -rf nullptr
	g++ -std=c++0x nullptr.cpp -o nullptr

unique:
	rm -rf unique
	g++ -std=c++0x unique_ptr.cpp -o unique

clean: 
	rm -rf *.o *.a

