PROJECT := ./

CXX = g++
CXXFLAGS = -O2 \
           -std=c++11 \
           -Wall \
           -Wno-sign-compare \
           -fno-omit-frame-pointer

LD_FLAGS = -static -lgsl -lgslcblas -lpthread

PROG_INC = $(PROJECT)/include
INC_FLAGS = -I$(PROG_INC)
  	  	
BASE_SRC = $(shell find $(PROJECT) -type f -name "*.cpp" -type f ! -name "stat_main.cpp")
BASE_OBJ = $(BASE_SRC:.cpp=.o)

PROG_HEADERS = $(shell find $(PROJECT) -type f -name "*.h")
PROG_SRC     = $(shell find $(PROJECT) -type f -name "*.cpp")
PROG_OBJ = $(PROG_SRC:.cpp=.o)

BIN_DIR = $(PROJECT)/bin
PROG = $(BIN_DIR)/snp_summary_stat

all: path \
	 snp_summary_stat

path: $(BIN_DIR)

$(BIN_DIR):
	mkdir -p $@

$(PROG): $(PROG_OBJ)
	$(CXX) $(PROG_OBJ) $(CXXFLAGS) $(INC_FLAGS) $(LD_FLAGS) -o $@

$(PROG_OBJ): %.o: %.cpp $(PROG_HEADERS)
	$(CXX) $(CXXFLAGS) $(INC_FLAGS) -c $< -o $@


snp_summary_stat: path $(PROG)

clean:
	rm -rf $(BIN_DIR) $(PROG_OBJ)

.PHONY: all path snp_summary_stat clean
