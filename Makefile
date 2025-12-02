# Compiler
CXX = g++
CXXFLAGS = -std=c++20 -Wall -Iinclude
LDFLAGS = -lsndfile

# Directories
SRC_DIR = src
BUILD_DIR = build

# Find all .cpp files
SRC = $(wildcard $(SRC_DIR)/*.cpp)

# .o files go into build/
OBJ = $(patsubst $(SRC_DIR)/%.cpp, $(BUILD_DIR)/%.o, $(SRC))

# Output executable
OUT = $(BUILD_DIR)/program

# Default rule
all: $(OUT)

# Link final binary
$(OUT): $(OBJ)
	$(CXX) $(OBJ) -o $(OUT) $(LDFLAGS)

# Compile .cpp → .o
$(BUILD_DIR)/%.o: $(SRC_DIR)/%.cpp
	@mkdir -p $(BUILD_DIR)
	$(CXX) $(CXXFLAGS) -c $< -o $@

clean:
	rm -rf $(BUILD_DIR)

.PHONY: all clean

