include ./Makefile.inc

TARGET_EXEC ?= main

BUILD_DIR ?= ./build_${COMPILER}
SRC_DIR ?= ./src
INC_DIR ?= ./include
TEST_CASE_INC_DIR ?= ./include/test_cases

SRCS := $(shell find $(SRC_DIR) -name "*.cpp" -or -name *.c -or -name *.s)
OBJS := $(SRCS:%=$(BUILD_DIR)/%.o)
DEPS := $(OBJS:.o=.d)

INCS ?= -I$(INC_DIR) -I$(TEST_CASE_INC_DIR) $(MT_FLAGS) $(C_OTHERS)
LIBS ?= 
INC_FLAGS := -I$(INC_DIR) -I$(TEST_CASE_INC_DIR)

CPPFLAGS ?= -O2 -Wall -MMD -MP -Wwrite-strings $(INCS) -MMD -MP
LDFLAGS ?= -MMD -MP $(LSCALAPACK) $(LBLAS) $(MT_FLAGS) $(L_OTHERS)
CXXFLAGS ?=
CFLAGS ?=

$(BUILD_DIR)/$(TARGET_EXEC): $(OBJS)
	$(CXX) $(OBJS) -o $@ $(LDFLAGS)

# c source
$(BUILD_DIR)/%.c.o: %.c
	mkdir -p $(dir $@)
	$(CC) $(CPPFLAGS) $(CFLAGS) -c $< -o $@

# c++ source
$(BUILD_DIR)/%.cpp.o: %.cpp
	mkdir -p $(dir $@)
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -c $< -o $@


.PHONY: clean

clean:
	$(RM) -r $(BUILD_DIR)

MKDIR_P ?= mkdir -p
