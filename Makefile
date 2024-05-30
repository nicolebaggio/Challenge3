CC = mpicxx
CXXFLAGS = -std=c++20
CPPFLAGS= -Wall -03 -I include

DOXYFILE = Doxyfile

SRCS = main.cpp 
OBJS = $(SRCS:.cpp=.o)

all: main 

main: $(OBJS)
	$(CXX) $(CXXFLAGS) $(OBJS) -o $@ -lmpi

%.o: %.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@


clean:
	rm *.o main 

doc:
	doxygen $(DOXYFILE)


