CXX := g++
STD := -std=c++14
CPPFLAGS := $(STD) -Iinclude
CXXFLAGS := $(STD) -Wall -O3 -Iinclude -fmax-errors=3
# CXXFLAGS := $(STD) -Wall -g -Iinclude -fmax-errors=3
LDFLAGS :=
LDLIBS :=

SRC := src
BIN := bin
BLD := .build
EXT := .cc

.PHONY: all clean

ifeq (0, $(words $(findstring $(MAKECMDGOALS), clean)))

LDFLAGS += $(shell sed -r 's/([^:]+)(:|$$)/ -L\1/g' <<< "$$LIBRARY_PATH")

ROOT_CXXFLAGS := $(shell root-config --cflags | sed 's/ -std=c++[^ ]\+ / /')
ROOT_LDFLAGS  := $(shell root-config --ldflags)
ROOT_LDLIBS   := $(shell root-config --libs)

# RPATH
rpath_script := ldd $(shell root-config --libdir)/libTreePlayer.so \
  | sed -nr 's|.*=> (.+)/.+\.so[.0-9]* \(.*|\1|p' \
  | sort -u \
  | sed -nr '/^(\/usr)?\/lib/!s/^/-Wl,-rpath=/p'
ROOT_LDLIBS += $(shell $(rpath_script))

C_angles1 := $(ROOT_CXXFLAGS)
L_angles1 := $(ROOT_LDLIBS)

C_angles := $(ROOT_CXXFLAGS)
L_angles := $(ROOT_LDLIBS)

C_fit1 := $(ROOT_CXXFLAGS)
L_fit1 := $(ROOT_LDLIBS) -lMinuit

C_fit := -fopenmp $(ROOT_CXXFLAGS)
L_fit := -fopenmp $(ROOT_LDLIBS) -lMinuit

C_mc_test := $(ROOT_CXXFLAGS)
L_mc_test := $(ROOT_LDLIBS) -lMinuit

C_draw1 := $(ROOT_CXXFLAGS)
L_draw1 := $(ROOT_LDLIBS)

C_draw := $(ROOT_CXXFLAGS)
L_draw := $(ROOT_LDLIBS)

C_pars := $(ROOT_CXXFLAGS)
L_pars := $(ROOT_LDLIBS)

C_llr := $(ROOT_CXXFLAGS)
L_llr := $(ROOT_LDLIBS)

SRCS := $(shell find $(SRC) -type f -name '*$(EXT)')
DEPS := $(patsubst $(SRC)/%$(EXT),$(BLD)/%.d,$(SRCS))

GREP_EXES := grep -rl '^ *int \+main *(' $(SRC) --include='*$(EXT)'
EXES := $(patsubst $(SRC)%$(EXT),$(BIN)%,$(shell $(GREP_EXES)))

all: $(EXES)

$(BIN)/angles1 $(BIN)/angles \
$(BIN)/fit1 $(BIN)/fit $(BIN)/mc_test \
$(BIN)/draw1 $(BIN)/draw $(BIN)/pars $(BIN)/llr \
: $(BLD)/program_options.o

-include $(DEPS)

.SECONDEXPANSION:

$(DEPS): $(BLD)/%.d: $(SRC)/%$(EXT) | $(BLD)/$$(dir %)
	$(CXX) $(CPPFLAGS) $(C_$*) -MM -MT '$(@:.d=.o)' $< -MF $@

$(BLD)/%.o: | $(BLD)
	$(CXX) $(CXXFLAGS) $(C_$*) -c $(filter %$(EXT),$^) -o $@

$(BIN)/%: $(BLD)/%.o | $(BIN)
	$(CXX) $(LDFLAGS) $(filter %.o,$^) -o $@ $(LDLIBS) $(L_$*)

$(BIN):
	mkdir -p $@

$(BLD)/%/:
	mkdir -p $@

endif

clean:
	@rm -rfv $(BLD) $(BIN)
