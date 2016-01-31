################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
C_SRCS += \
../src/kernel/bstworst.c \
../src/kernel/change.c \
../src/kernel/ckpoint.c \
../src/kernel/crossovr.c \
../src/kernel/ephem.c \
../src/kernel/eval.c \
../src/kernel/event.c \
../src/kernel/exch.c \
../src/kernel/fitness.c \
../src/kernel/genspace.c \
../src/kernel/gp.c \
../src/kernel/individ.c \
../src/kernel/main.c \
../src/kernel/memory.c \
../src/kernel/mutate.c \
../src/kernel/output.c \
../src/kernel/params.c \
../src/kernel/populate.c \
../src/kernel/pretty.c \
../src/kernel/random.c \
../src/kernel/reproduc.c \
../src/kernel/select.c \
../src/kernel/tournmnt.c \
../src/kernel/tree.c 

OBJS += \
./src/kernel/bstworst.o \
./src/kernel/change.o \
./src/kernel/ckpoint.o \
./src/kernel/crossovr.o \
./src/kernel/ephem.o \
./src/kernel/eval.o \
./src/kernel/event.o \
./src/kernel/exch.o \
./src/kernel/fitness.o \
./src/kernel/genspace.o \
./src/kernel/gp.o \
./src/kernel/individ.o \
./src/kernel/main.o \
./src/kernel/memory.o \
./src/kernel/mutate.o \
./src/kernel/output.o \
./src/kernel/params.o \
./src/kernel/populate.o \
./src/kernel/pretty.o \
./src/kernel/random.o \
./src/kernel/reproduc.o \
./src/kernel/select.o \
./src/kernel/tournmnt.o \
./src/kernel/tree.o 

C_DEPS += \
./src/kernel/bstworst.d \
./src/kernel/change.d \
./src/kernel/ckpoint.d \
./src/kernel/crossovr.d \
./src/kernel/ephem.d \
./src/kernel/eval.d \
./src/kernel/event.d \
./src/kernel/exch.d \
./src/kernel/fitness.d \
./src/kernel/genspace.d \
./src/kernel/gp.d \
./src/kernel/individ.d \
./src/kernel/main.d \
./src/kernel/memory.d \
./src/kernel/mutate.d \
./src/kernel/output.d \
./src/kernel/params.d \
./src/kernel/populate.d \
./src/kernel/pretty.d \
./src/kernel/random.d \
./src/kernel/reproduc.d \
./src/kernel/select.d \
./src/kernel/tournmnt.d \
./src/kernel/tree.d 


# Each subdirectory must supply rules for building sources it contributes
src/kernel/%.o: ../src/kernel/%.c
	@echo 'Building file: $<'
	@echo 'Invoking: Cross GCC Compiler'
	gcc -O0 -g3 -Wall -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


