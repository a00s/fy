################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CU_SRCS += \
../src/aminos.cu \
../src/aula1.cu \
../src/learn.cu \
../src/pdb.cu \
../src/util.cu 

CU_DEPS += \
./src/aminos.d \
./src/aula1.d \
./src/learn.d \
./src/pdb.d \
./src/util.d 

OBJS += \
./src/aminos.o \
./src/aula1.o \
./src/learn.o \
./src/pdb.o \
./src/util.o 


# Each subdirectory must supply rules for building sources it contributes
src/%.o: ../src/%.cu
	@echo 'Building file: $<'
	@echo 'Invoking: NVCC Compiler'
	/usr/local/cuda-6.5/bin/nvcc -Imysqlcppconn -G -g -O0 -gencode arch=compute_30,code=sm_30  -odir "src" -M -o "$(@:%.o=%.d)" "$<"
	/usr/local/cuda-6.5/bin/nvcc -Imysqlcppconn -G -g -O0 --compile --relocatable-device-code=false -gencode arch=compute_30,code=compute_30 -gencode arch=compute_30,code=sm_30  -x cu -o  "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


