# rootcling dictionary
DICT := $(LARGEX)Dict.cxx
PCM := $(LARGEX)Dict_rdict.pcm
LINKDEF = LinkDef.h


# source code (with $(DICT) and $(LINKDEF) moved to end of lists for rootcling)
SOURCES := $(filter-out $(DICT), $(wildcard *.cxx) $(wildcard sfset/*.cxx) $(wildcard interp/*.cxx)) $(DICT)
HEADERS := $(filter-out $(LINKDEF), $(wildcard *.h) $(wildcard sfset/*.h) $(wildcard interp/*.h) $(wildcard interp/*.ipp)) $(LINKDEF)


#-----------------------------------------------

$(LARGEXOBJ): $(DICT) $(HEADERS) $(SOURCES)
	@echo "----- build $@ -----"
	$(CXX) $(SOURCES) -shared -o ../$@ $(FLAGS) $(SANFLAGS) $(DEPS) $(LIBS)

$(DICT): $(HEADERS)
	@echo "----- generate dictionary $@ -----"
	@rootcling -f $@ $(DEPS) $^
	mv $(PCM) ../

clean:
	$(RM) ../$(LARGEXOBJ) $(DICT) ../$(PCM)

