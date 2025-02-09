
UNAME_S := $(shell uname -s)
ifeq ($(UNAME_S),Darwin)
  CXX = g++-14
  CXXFLAGS = -O3 -Wall -Wextra -std=c++17 -I ./toofus
  GLLINK = `pkg-config --libs glut` -framework OpenGL
	GLFLAGS = `pkg-config --cflags glut`
else
  CXX = g++
  CXXFLAGS = -O3 -Wall -std=c++17 -I ./toofus
  GLLINK = -lGLU -lGL -L/usr/X11R6/lib -lglut -lXmu -lXext -lX11 -lXi
endif

# The list of source files
SOURCES = Particle.cpp Pipe.cpp Interaction.cpp InteractionPipe.cpp Loading.cpp PeriodicCell.cpp BuriedPipe.cpp

# Each cpp file listed below corresponds to an object file
OBJECTS = $(SOURCES:%.cpp=%.o)

.PHONY: all clean clone_toofus

all: BuriedPipe see

clean:
	echo "\033[0;32m-> Remove object files\033[0m"
	rm -f *.o
	echo "\033[0;32m-> Remove BuriedPipe and see\033[0m"
	rm -f BuriedPipe see
	echo "\033[0;32m-> Remove libBuriedPipe.a\033[0m"
	rm -f libBuriedPipe.a

clean+: clean
	echo "\033[0;32m-> Remove local folder toofus\033[0m"
	rm -rf ./toofus
	
clone_toofus:
	@if [ ! -d "./toofus" ]; then \
		echo "\033[0;32m-> CLONING ToOfUs (Tools often used)\033[0m"; \
		git clone https://github.com/richefeu/toofus.git; \
	fi

%.o: %.cpp clone_toofus
	@echo "\033[0;32m-> COMPILING OBJECT" $@ "\033[0m"
	$(CXX) $(CXXFLAGS) -c $< -o $@

libBuriedPipe.a: $(OBJECTS)
	@echo "\033[0;32m-> BUILDING LIBRARY" $@ "\033[0m"
	ar rcs $@ $^
	
BuriedPipe: run.cpp libBuriedPipe.a
	@echo "\033[0;32m-> BUILDING APPLICATION" $@ "\033[0m"
	$(CXX) $(CXXFLAGS) -c $< -o run.o
	$(CXX) -o $@ run.o libBuriedPipe.a
	
see: see.cpp libBuriedPipe.a
	@echo "\033[0;32m-> BUILDING APPLICATION" $@ "\033[0m"
	$(CXX) $(CXXFLAGS) -c $< -o see.o $(GLFLAGS)
	$(CXX) -o $@ see.o libBuriedPipe.a $(GLLINK)
	
	
	
	