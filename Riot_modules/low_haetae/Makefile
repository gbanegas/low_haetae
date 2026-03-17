# name of your application
APPLICATION = blinky

CFLAGS += -DHAETAE_MODE=3


MAIN_THREAD_SIZE = THREAD_STACKSIZE_DEFAULT+70000
CFLAGS += -DTHREAD_STACKSIZE_MAIN=\($(MAIN_THREAD_SIZE)\)
CFLAGS += -DHAETAE_THREAD_STACKSIZE=70000
# If no BOARD is found in the environment, use this default:
BOARD ?= native

# This has to be the absolute path to the RIOT base directory:
RIOTBASE ?= $(CURDIR)/../../..

# Comment this out to disable code in RIOT that does safety checking
# which is not needed in a production environment but helps in the
# development process:
DEVELHELP ?= 1

# Change this to 0 show compiler invocation lines by default:
QUIET ?= 1

USEMODULE += xtimer
USEMODULE += random
USEMODULE += auto_init_random


# Use a peripheral timer for the delay, if available
FEATURES_OPTIONAL += periph_timer

include $(RIOTBASE)/Makefile.include
