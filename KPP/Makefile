# Set here the desired C compiler and its optimization options
CC   = gcc #g++
COPT = -O3 -Wall  

ifndef ROOT
	# Set root directory
	ROOT := ..
endif

# Set src folder
SRC_DIR := ./src

CURR_DIR := KPP

# Define any libraries to link into executable
LDFLAGS := -lm

ifndef BUILD
	BUILD := $(ROOT)/build
endif
ifndef APP_DIR
	APP_DIR := $(BUILD)/apps
endif
OBJ_DIR := $(BUILD)/objects/$(CURR_DIR)

# Define output library:
LIB := libKpp.a

# Define lib folder:
LIB_DIR := $(ROOT)/lib

# Define the executable file:
TARGET := KPP

# Define any directories containing header files
INCLUDE := -I./include

HEADERS = $(wildcard ./include/*.h)

# To create Matlab gateway routines
# Note: use $(CC) as the mex C compiler
MEX  = mex

SPSRC = $(SRC_DIR)/KPP_JacobianSP.c \
		$(SRC_DIR)/KPP_HessianSP.c  \
		$(SRC_DIR)/KPP_StoichiomSP.c

SPOBJ = ${SPSRC:$(SRC_DIR)/%.c=$(OBJ_DIR)/%.o}

# Define the C source files
SRC =   $(SRC_DIR)/KPP_Main.c      $(SRC_DIR)/KPP_Integrator.c   \
		$(SRC_DIR)/KPP_Function.c  $(SRC_DIR)/KPP_Initialize.c   \
		$(SRC_DIR)/KPP_Jacobian.c  $(SRC_DIR)/KPP_LinearAlgebra.c\
		$(SRC_DIR)/KPP_Rates.c     $(SRC_DIR)/KPP_Hessian.c      \
		$(SRC_DIR)/KPP_Stoichiom.c $(SRC_DIR)/KPP_Util.c         \
		$(SRC_DIR)/KPP_Monitor.c   $(SRC_DIR)/KPP_Photol.c

# Define the C object files
OBJ = ${SRC:$(SRC_DIR)/%.c=$(OBJ_DIR)/%.o}

.PHONY: all build clean mex

all:    build $(APP_DIR)/$(TARGET) $(LIB_DIR)/$(LIB)

$(OBJ_DIR)/%.o: $(SRC_DIR)/%.c $(HEADERS)
	@mkdir -p $(@D)
	$(CC) $(COPT) $(INCLUDE) -o $@ -c $<

$(APP_DIR)/$(TARGET): $(HEADERS) $(SPOBJ) $(OBJ)
	@mkdir -p $(@D)

#	$(CC) $(COPT) $(INCLUDE) -o $(APP_DIR)/$(TARGET) $(SPOBJ) $(OBJ) $(LDFLAGS)
#	@echo "---> Executable file $(TARGET) has been compiled! <---"

$(LIB_DIR)/$(LIB): $(SPOBJ) $(OBJ)
	@mkdir -p $(LIB_DIR)
	@$(AR) cr $@ $^
	@echo "Library $(LIB) has been created"

mex:    $(HEADERS) $(SPOBJ) $(OBJ)
	$(MEX) CC#$(CC) -O KPP_mex_Fun.c     -lm $(SPOBJ) $(OBJ)
	$(MEX) CC#$(CC) -O KPP_mex_Jac_SP.c  -lm $(SPOBJ) $(OBJ)
	$(MEX) CC#$(CC) -O KPP_mex_Hessian.c -lm $(SPOBJ) $(OBJ)

build:
	@mkdir -p $(APP_DIR)
	@mkdir -p $(OBJ_DIR)
	@mkdir -p $(LIB_DIR)

clean:
	@echo "---> Making clean in KPP <---"
	-@rm -rv $(OBJ_DIR)/*
	-@rm -rv $(LIB_DIR)/*

