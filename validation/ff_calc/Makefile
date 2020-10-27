VPATH := build/lib:inc:src:$(VPATH)

OS := $(shell uname)
INC := -Iinc

# Define compiler command
ifeq ($(OS),Darwin)
	CC_CMD = ${CC} -dynamiclib -flat_namespace -I$(INC)
	LIB_EXT = .dylib
else
	CC_CMD = ${CC} -fPIC -shared $(INC)
	LIB_EXT = .so
endif

# Source files and targets
CPP_FILES := $(notdir $(wildcard ./src/*.cpp))
LIB_FILES := $(addprefix lib, $(CPP_FILES:.cpp=$(LIB_EXT)))


install: prepare $(LIB_FILES)

# General rule
lib%.so: %.hpp %.cpp
	$(CC_CMD) -o build/lib/$@


# Helper targets
.PHONY: prepare
prepare:
	@mkdir -p build/lib

.PHONY: clean
clean:
	@rm -rf build

.PHONY: debug
debug:
	echo $(CC_CMD)
	echo $(CPP_FILES)
	echo $(LIB_FILES)
