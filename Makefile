# Compiler and flags
CXX := g++
CXXFLAGS := -std=c++20 -Wall -Wextra -O3 -Iinclude -I/usr/local/include -DTEST_DATA_DIR=\"$(CURDIR)/tests/data\"
LDFLAGS := -L/usr/local/lib -lgtest -lgtest_main -pthread

# Source and output directories
SRC_DIR := src
INC_DIR := include
TEST_DIR := tests
OBJ_DIR := obj
SRC := $(wildcard $(SRC_DIR)/*.cpp)
OBJ := $(patsubst $(SRC_DIR)/%.cpp,$(OBJ_DIR)/%.o,$(SRC))
LIB := hjbsolver.a
TEST_SRC := $(wildcard $(TEST_DIR)/*.cpp)
TEST_OBJ := $(patsubst $(TEST_DIR)/%.cpp,$(OBJ_DIR)/%.o,$(TEST_SRC))
TEST_BIN := test_runner

# Default target
all: $(LIB)

# Build static library
$(LIB): $(OBJ)
	ar rcs $@ $^

# Build object files for the library
$(OBJ_DIR)/%.o: $(SRC_DIR)/%.cpp | $(OBJ_DIR)
	$(CXX) $(CXXFLAGS) -c $< -o $@

# Build object files for tests
$(OBJ_DIR)/%.o: $(TEST_DIR)/%.cpp | $(OBJ_DIR)
	$(CXX) $(CXXFLAGS) -c $< -o $@

# Create object directory if it doesn't exist
$(OBJ_DIR):
	mkdir -p $@

# Build test executable
$(TEST_BIN): $(LIB) $(TEST_OBJ)
	$(CXX) -o $@ $(TEST_OBJ) $(LIB) $(LDFLAGS)

# Run tests
test: $(TEST_BIN)
	./$(TEST_BIN)

# Clean target
clean:
	rm -rf $(OBJ_DIR) $(LIB) $(TEST_BIN)

# Phony targets
.PHONY: all clean test