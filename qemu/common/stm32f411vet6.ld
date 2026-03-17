/* common/stm32f411vet6.ld
 * STM32F411VET6: 512K Flash, 128K SRAM
 */

MEMORY
{
  rom (rx)  : ORIGIN = 0x08000000, LENGTH = 512K
  ram (rwx) : ORIGIN = 0x20000000, LENGTH = 128K
}

/* Top of stack (used by vector table) */
PROVIDE(_stack = ORIGIN(ram) + LENGTH(ram));
PROVIDE(_estack = ORIGIN(ram) + LENGTH(ram));

SECTIONS
{
  /* Vector table + code + const data in FLASH */
  .text :
  {
    KEEP(*(.vectors))          /* libopencm3 vector table */
    *(.text*)
    *(.rodata*)
    *(.glue_7)
    *(.glue_7t)

    /* C++/libgcc metadata (harmless if empty) */
    *(.eh_frame*)

    . = ALIGN(4);

    /* Exception unwinding tables (must be in FLASH before .data LMA) */
    *(.ARM.extab* .gnu.linkonce.armextab.*)
    . = ALIGN(4);
    __exidx_start = .;
    *(.ARM.exidx* .gnu.linkonce.armexidx.*)
    __exidx_end = .;

    . = ALIGN(4);
    _etext = .;                /* End of FLASH image used for .data load */
  } > rom

  /* Init/fini arrays (optional, but safe) */
  .preinit_array :
  {
    . = ALIGN(4);
    KEEP(*(.preinit_array*))
  } > rom

  .init_array :
  {
    . = ALIGN(4);
    KEEP(*(.init_array*))
  } > rom

  .fini_array :
  {
    . = ALIGN(4);
    KEEP(*(.fini_array*))
  } > rom

  /* .data in RAM, loaded from FLASH right after .text/.exidx/... */
  .data : AT(_etext)
  {
    . = ALIGN(4);
    _sdata = .;
    _data  = .;                /* compatibility */
    *(.data*)
    . = ALIGN(4);
    _edata = .;
  } > ram

  /* Start address of initialised data in FLASH */
  _sidata = LOADADDR(.data);

  /* .bss in RAM */
  .bss (NOLOAD) :
  {
    . = ALIGN(4);
    _sbss = .;
    _bss  = .;                 /* compatibility */
    *(.bss*)
    *(COMMON)
    . = ALIGN(4);
    _ebss = .;
  } > ram

  /* Heap starts here */
  . = ALIGN(4);
  _end = .;
  PROVIDE(end = .);            /* <-- this fixes your undefined reference */
}

PROVIDE(end = .);
_estack = ORIGIN(ram) + LENGTH(ram);  /* top of stack */
_sstack = _estack - 0x1000;           /* bottom of stack, 4 KB stack */
PROVIDE(_stack = _estack);

PROVIDE(__stack_top__ = _estack);
PROVIDE(__stack_bottom__ = _sstack);
