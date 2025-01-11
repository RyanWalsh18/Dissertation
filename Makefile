# Compiler
CXX = clang++

# Compiler flags
CXXFLAGS = -std=c++14 -Wall -g -Xpreprocessor -fopenmp

# Include and library directories
INCLUDE_DIR = -I include -I/usr/local/opt/libomp/include
LIB_DIR = -L lib -L/usr/local/opt/libomp/lib

# Libraries
LIBS = -lSDL2 -lomp

# Source files
SRC_FILES = main.cpp simulation.cpp body.cpp quadTree.cpp

# Object files derived from source files
OBJS = $(SRC_FILES:.cpp=.o)
OBJS_WITHOUT_MAIN = simulation.o body.o quadTree.o

# Header files
HEADERS = body.h simulation.h quadTree.h accuracyTests.h

# Executable name for the main application
EXEC = simulation

# Test files 
# TEST_SRC = $(wildcard test_*.cpp)
TEST_SRC = test_integrationTests.cpp

# Main test file
MAIN_TEST_FILE = main_tests.cpp

# Test executable name
TEST_EXEC = tests

# Accuracy Test Executable Name
ACCURACY_EXEC = generate_accuracy_data

# All targets
all: $(EXEC) $(TEST_EXEC) $(ACCURACY_EXEC)

# Main application executable
$(EXEC): $(OBJS)
	$(CXX) $(CXXFLAGS) $(OBJS) -o $(EXEC) $(INCLUDE_DIR) $(LIB_DIR) $(LIBS)

# Test executable
$(TEST_EXEC): $(MAIN_TEST_FILE) $(TEST_SRC) $(OBJS_WITHOUT_MAIN)
	$(CXX) $(CXXFLAGS) $(MAIN_TEST_FILE) $(TEST_SRC) $(OBJS_WITHOUT_MAIN) -o $@ $(INCLUDE_DIR) $(LIB_DIR) $(LIBS)

# Accuracy Test Executable
$(ACCURACY_EXEC): accuracyTests.o $(OBJS_WITHOUT_MAIN)
	$(CXX) $(CXXFLAGS) accuracyTests.o $(OBJS_WITHOUT_MAIN) -o $@ $(INCLUDE_DIR) $(LIB_DIR) $(LIBS)


# Compile source files to object files
%.o: %.cpp $(HEADERS)
	$(CXX) $(CXXFLAGS) -c $< -o $@ $(INCLUDE_DIR)

# Test compilation 
test: $(TEST_EXEC)

accuracy_data: $(ACCURACY_EXEC)

# Clean up
clean:
	rm -f $(OBJS) $(EXEC) $(TEST_EXEC) *.o

.PHONY: all clean test accuracy_data