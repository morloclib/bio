all:
	morloc make foo.loc

.PHONY:
c:
	g++ --std=c++17 -g -Wall -o pool-cpp.out pool.cpp -I/home/z/.morloc/include -I/home/z/.morloc/src/morloc/plain/morloclib/matrix/eigen -I/home/z/.morloc/src/morloc/plain/morloclib/bio/tree -I/home/z/.morloc/src/morloc/plain/morloclib/cppbase -I/home/z/.morloc/src/morloc/plain/morloclib/bio/algo

.PHONY:
run:
	./nexus.py upgma '[[0,17,21,31,23],[17,0,30,34,21],[21,30,0,28,39],[31,34,28,0,43],[23,21,39,43,0]]'

# (((a,b),e),(c,d))
# (((0,1),4),(2,3))
