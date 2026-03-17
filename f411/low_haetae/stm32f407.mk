# Always build the firmware by default (not just libopencm3 artifacts)
.DEFAULT_GOAL := all

OPENCM3_DIR ?= ../libopencm3

# STM32F411E-DISCO => STM32F411VET6 (F4 family)
LIBNAME := opencm3_stm32f4

ARCH_FLAGS := -mthumb -mcpu=cortex-m4 -mfloat-abi=hard -mfpu=fpv4-sp-d16

# Board define used by common/hal-opencm3.c
DEFINES := -DSTM32F4 -DSTM32F411xE -DSTM32F411E_DISCO
# -DSIGN_RECOMPUTE_Y=1

CFLAGS += -O3  -DBENCH_CYCLES=1 -DBENCH_ITERS=100\
  -Wall -Wextra -Wimplicit-function-declaration \
  -Wredundant-decls -Wmissing-prototypes -Wstrict-prototypes \
  -Wundef -Wshadow \
  -I$(OPENCM3_DIR)/include \
  -Icommon \
  -fno-common $(ARCH_FLAGS) -MD $(DEFINES)

# Use the linker script you actually have in common/
# (your tree shows common/stm32f411vet6.ld)
LDSCRIPT := common/stm32f411vet6.ld

LDLIBS  += -l$(LIBNAME)
LDFLAGS += -L$(OPENCM3_DIR)/lib \
  --specs=nosys.specs \
  -Wl,--wrap=_sbrk \
  -nostartfiles -ffreestanding \
  -T$(LDSCRIPT) \
  $(ARCH_FLAGS)

# Build libopencm3 before compiling anything that includes its headers
OPENCM3_A := $(OPENCM3_DIR)/lib/lib$(LIBNAME).a

$(OPENCM3_A):
	$(MAKE) -C $(OPENCM3_DIR)

obj/%.c.o: | $(OPENCM3_A)
obj/%.S.o: | $(OPENCM3_A)

# Your project needs common HAL objects linked in
LINKDEPS += obj/common/hal-opencm3.c.o obj/common/randombytes.c.o obj/common/fault.c.o

FW_NAME := stm32f411-test
all: bin/$(FW_NAME).bin

# Convenience flash target (uses ST-LINK)
FLASH_ADDR ?= 0x08000000
flash: bin/$(FW_NAME).bin
	st-flash --reset write $< $(FLASH_ADDR)