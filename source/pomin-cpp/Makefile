CXX = g++
CXXFLAGS = -Wall -g -std=c++11
TARGETS = main.o pomin.o particle.o RK4.o

main: $(TARGETS)
	$(CXX) $(CXXFLAGS) -o main $(TARGETS)

clean:
	rm main $(TARGETS)
