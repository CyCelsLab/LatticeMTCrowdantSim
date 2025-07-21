CC=g++
CFLAGS=-c -Wall -std=c++20
SRC=./src/
mmc: ./build/ ./build/main.o ./build/diffuse.o ./build/depolymerise.o ./build/polymerise.o ./build/test.o ./build/print.o ./build/constructor.o ./build/populating.o ./build/kinetics.o ./build/preseed.o
	mkdir -p build bin 
#	g++ -L/build/ main.o diffuse.o depolymerise.o polymerise.o test.o print.o constructor.o populating.o -o ./bin/mmc
	g++ ./build/main.o ./build/diffuse.o ./build/depolymerise.o ./build/polymerise.o ./build/test.o ./build/print.o ./build/constructor.o ./build/populating.o ./build/kinetics.o ./build/preseed.o -o ./bin/mmc -lyaml-cpp

./build/:
	mkdir build
./build/main.o:
	$(CC) ${CFLAGS} ${SRC}main.cpp -o ./build/main.o -lyaml-cpp
./build/diffuse.o: ${SRC}class.hpp ${SRC}random.hpp
	$(CC) ${CFLAGS} ${SRC}diffuse.cpp -o ./build/diffuse.o
./build/polymerise.o: ${SRC}class.hpp ${SRC}random.hpp
	$(CC) ${CFLAGS} ${SRC}polymerise.cpp -o ./build/polymerise.o
./build/depolymerise.o: ${SRC}class.hpp ${SRC}random.hpp
	$(CC) ${CFLAGS} ${SRC}depolymerise.cpp -o ./build/depolymerise.o
./build/test.o: ${SRC}class.hpp 
	$(CC) ${CFLAGS} ${SRC}test.cpp -o ./build/test.o
./build/print.o:${SRC}class.hpp
	$(CC) ${CFLAGS} ${SRC}print.cpp -o ./build/print.o
./build/populating.o: ${SRC}class.hpp ${SRC}random.hpp 
	${CC} ${CFLAGS} ${SRC}populating.cpp -o ./build/populating.o
./build/constructor.o:
	${CC} ${CFLAGS} ${SRC}constructor.cpp -o ./build/constructor.o
./build/kinetics.o:
	${CC} ${CFLAGS} ${SRC}kinetics.cpp -o ./build/kinetics.o
./build/preseed.o:
	${CC} ${CFLAGS} ${SRC}preseed.cpp -o ./build/preseed.o
clean:
	rm -rf ./bin/*
	rm -rf ./build/*
done:
	make clean
	make -j 10
	
install:
	git clone https://github.com/jbeder/yaml-cpp.git
	cd yaml-cpp && mkdir build && cd build && cmake .. && sudo make install
	make done
	
