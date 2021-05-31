fastLK: fastLK.o
	g++ -o fastLK fastLK.o

fastLK.o: fastLK.cpp 
	g++ -c fastLK.cpp -o fastLK.o -Ofast -I eigen3

clean:
	rm fastLK fastLK.o
