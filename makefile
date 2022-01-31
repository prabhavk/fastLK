fastLK: fastLK.o
	g++ -o fastLK fastLK.o

fastLK.o: fastLK.cpp 
	g++ -c fastLK.cpp -o fastLK.o -Ofast -I eigen3 -I . -std=c++17

clean:
	rm fastLK fastLK.o
