# Set here the desired C compiler and its optimization options
CC   = gcc
COPT = -O3 -Wall  

ifndef ROOT
	# Set root directory
	ROOT := ..
endif

# Set src folder
SRC_FOLDER := ./src

CURR_FOLDER := KPP

# Define any libraries to link into executable
LDFLAGS := -lm

ifndef BUILD
	BUILD := $(ROOT)/build
endif
ifndef APP_DIR
	APP_DIR := $(BUILD)/apps
endif
OBJ_DIR := $(BUILD)/objects/$(CURR_FOLDER)

# Define the executable file:
TARGET := KPP

# Define any directories containing header files
INCLUDE := -I./include

HEADERS = $(wildcard ./include/*.h)

# To create Matlab gateway routines
# Note: use $(CC) as the mex C compiler
MEX  = mex

SPSRC = $(SRC_FOLDER)/KPP_JacobianSP.c \
		$(SRC_FOLDER)/KPP_HessianSP.c  \
		$(SRC_FOLDER)/KPP_StoichiomSP.c

SPOBJ = ${SPSRC:$(SRC_FOLDER)/%.c=$(OBJ_DIR)/%.o}

# Define the C source files
SRC =   $(SRC_FOLDER)/KPP_Main.c      $(SRC_FOLDER)/KPP_Integrator.c   \
		$(SRC_FOLDER)/KPP_Function.c  $(SRC_FOLDER)/KPP_Initialize.c   \
		$(SRC_FOLDER)/KPP_Jacobian.c  $(SRC_FOLDER)/KPP_LinearAlgebra.c\
		$(SRC_FOLDER)/KPP_Rates.c     $(SRC_FOLDER)/KPP_Hessian.c      \
		$(SRC_FOLDER)/KPP_Stoichiom.c $(SRC_FOLDER)/KPP_Util.c         \
		$(SRC_FOLDER)/KPP_Monitor.c   $(SRC_FOLDER)/KPP_Photol.c

# Define the C object files
OBJ = ${SRC:$(SRC_FOLDER)/%.c=$(OBJ_DIR)/%.o}

all:    build $(APP_DIR)/$(TARGET)

$(OBJ_DIR)/%.o: $(SRC_FOLDER)/%.c $(HEADERS)
	@mkdir -p $(@D)
	$(CC) $(COPT) $(INCLUDE) -o $@ -c $<

$(APP_DIR)/$(TARGET): $(HEADERS) $(SPOBJ) $(OBJ)
	mkdir -p $(@D)
	$(CC) $(COPT) $(INCLUDE) -o $(APP_DIR)/$(TARGET) $(SPOBJ) $(OBJ) $(LDFLAGS)
	@echo "---> Executable file $(TARGET) has been compiled! <---"

mex:    $(HEADERS) $(SPOBJ) $(OBJ)
	$(MEX) CC#$(CC) -O KPP_mex_Fun.c     -lm $(SPOBJ) $(OBJ)
	$(MEX) CC#$(CC) -O KPP_mex_Jac_SP.c  -lm $(SPOBJ) $(OBJ)
	$(MEX) CC#$(CC) -O KPP_mex_Hessian.c -lm $(SPOBJ) $(OBJ)

.PHONY: all build clean 

build:
	@mkdir -p $(APP_DIR)
	@mkdir -p $(OBJ_DIR)

clean:
	@echo "---> Making clean in KPP <---"
	-@rm -rv $(OBJ_DIR)/*
	-@rm -rv $(APP_DIR)/*

