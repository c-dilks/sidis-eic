include config.mk


# local dependencies
DEPS += -Isrc
#LIBS += -L. -l$(LARGEX)


# assume each .cpp file has main and build corresponding .exe executable
SOURCES := $(basename $(wildcard *.cpp))
EXES := $(addsuffix .exe, $(SOURCES))


#--------------------------------------------

export
all: 
	@cd mstwpdf; make
	@cd src; make
	make exe

asan: SANFLAGS = -g -fno-omit-frame-pointer -fsanitize=address -static-libasan
asan: all

exe: $(EXES)

%.exe: %.o
	@echo "--- make executable $@"
	$(CXX) -o $@ $< ./$(LARGEXOBJ) $(LIBS)

%.o: %.cpp
	@echo "----- build $@ -----"
	$(CXX) -c $^ -o $@ $(FLAGS) $(DEPS)

clean:
	@cd mstwpdf; make clean
	@cd src; make clean
#$(RM) $(EXES)
