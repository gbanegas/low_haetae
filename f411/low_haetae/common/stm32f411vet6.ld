/* STM32F411VET6 (STM32F411E-DISCO)
 * Flash: 512K @ 0x08000000
 * SRAM : 128K @ 0x20000000
 *
 * This linker script is compatible with libopencm3's cm3/vector.c startup.
 */

ENTRY(reset_handler)

MEMORY
{
  rom (rx)  : ORIGIN = 0x08000000, LENGTH = 512K
  ram (rwx) : ORIGIN = 0x20000000, LENGTH = 128K
}


SECTIONS
{
  /* Vector table (from libopencm3's vector.o) */
  .vectors :
  {
    KEEP(*(.vectors))
  } > rom

  /* Code + rodata */
  .text :
  {
    *(.text*)
    *(.rodata*)
    *(.glue_7)
    *(.glue_7t)
    *(.eh_frame*)
    KEEP(*(.init))
    KEEP(*(.fini))
    . = ALIGN(4);
  } > rom

  /* Exception unwinding tables must live in FLASH and must be placed
     BEFORE we define the .data load address, otherwise AT(...) may overlap. */
  .ARM.extab :
  {
    *(.ARM.extab* .gnu.linkonce.armextab.*)
    . = ALIGN(4);
  } > rom

  .ARM.exidx :
  {
    __exidx_start = .;
    *(.ARM.exidx* .gnu.linkonce.armexidx.*)
    __exidx_end = .;
    . = ALIGN(4);
  } > rom

  /* C runtime arrays (needed by libopencm3 reset_handler) */
  .preinit_array :
  {
    __preinit_array_start = .;
    KEEP(*(.preinit_array*))
    __preinit_array_end = .;
    . = ALIGN(4);
  } > rom

  .init_array :
  {
    __init_array_start = .;
    KEEP(*(SORT(.init_array.*)))
    KEEP(*(.init_array*))
    __init_array_end = .;
    . = ALIGN(4);
  } > rom

  .fini_array :
  {
    __fini_array_start = .;
    KEEP(*(SORT(.fini_array.*)))
    KEEP(*(.fini_array*))
    __fini_array_end = .;
    . = ALIGN(4);
  } > rom

  /* End of everything that lives in FLASH */
  . = ALIGN(4);
  _etext = .;

  /* .data lives in RAM, but its load image is in FLASH at _etext */
  .data : AT(_etext)
  {
    . = ALIGN(4);
    _data = .;
    *(.data*)
    . = ALIGN(4);
    _edata = .;
  } > ram

  /* What libopencm3's reset handler expects */
  _data_loadaddr = LOADADDR(.data);

  /* .bss lives in RAM */
  .bss (NOLOAD) :
  {
    . = ALIGN(4);
    _bss = .;
    *(.bss*)
    *(COMMON)
    . = ALIGN(4);
    _ebss = .;
  } > ram

  /* Heap start symbol (your code referenced `end` before) */
  . = ALIGN(4);
  PROVIDE(end = .);
  PROVIDE(_end = .);

  /DISCARD/ :
  {
    *(.note*)
    *(.comment*)
  }
}

PROVIDE(end = .);
_estack = ORIGIN(ram) + LENGTH(ram);  /* top of stack */
_sstack = _estack - 0x1000;           /* bottom of stack, 4 KB stack */
PROVIDE(_stack = _estack);

PROVIDE(__stack_top__ = _estack);
PROVIDE(__stack_bottom__ = _sstack);
