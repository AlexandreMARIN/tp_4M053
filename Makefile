INCLDIR	:= include
OBJDIR	:= obj
SRCDIR	:= src
BINDIR	:= bin

CXX     := g++
VPATH	:=
LDFLAGS :=
LIBRARY :=
CXXFLAGS  := -std=c++11 -Wall -I $(INCLDIR)

#Source and object files (automatic)
SRCS = $(wildcard $(SRCDIR)/*.cpp)
OBJS = $(subst $(SRCDIR)/,$(OBJDIR)/, $(subst .cpp,.o, $(SRCS)))

# Define here your main source files separated by spaces (without suffix!)
EXEC = test_Vect test_Matrix test_iter_solv test_iter_solv2 test_iter_solv3 test_LU test_Cholesky test_LU2 test_opt test_bad_cond test_iter_solv4 test_conj_grad test_COO test_iter_solv5

#Phony = do not represent a file
#.PHONY: all
all : makedir $(EXEC)

# For multiple binaries
$(EXEC) : %: %.cpp $(OBJS)
	$(CXX) $(CXXFLAGS) -o $(BINDIR)/$@ $^

$(OBJDIR)/%.o : $(SRCDIR)/%.cpp
	$(CXX) $(CXXFLAGS) -c -o $@ $<

#Clean: delete every binaries and object files
.PHONY: clean
clean :
	rm -rf $(OBJDIR)/*
	rm -rf $(BINDIR)/*
#Building folders (-p : no error if folder do not exist)
.PHONY: makedir
makedir :
	mkdir -p $(BINDIR)
	mkdir -p $(OBJDIR)

#For some debug
.PHONY: print
print :
	@echo $(SRCS)
	@echo $(OBJS)

#Remarks:
# $@ : filename representing the target
# $< : filename of the first prerequisite
# $^ : filenames of all prerequisites, separated by spaces. Dupplicated are removed.
# $? : names of all prerequisites that are newer than the target, separated by spaces
# $+ : similar to $^ but include dupplicates
# $* : stem of target file (filename without suffix)
