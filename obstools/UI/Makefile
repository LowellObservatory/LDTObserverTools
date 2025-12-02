# -*- coding: utf-8 -*-
#
#  This file is part of LDTObserverTools.
#
#   This Source Code Form is subject to the terms of the Mozilla Public
#   License, v. 2.0. If a copy of the MPL was not distributed with this
#   file, You can obtain one at http://mozilla.org/MPL/2.0/.
#
#  Created on 18-Jul-2025
#
#  @author: tbowers
#
#
# Makefile to generate PyQt6 python files from .ui XML-formatted templates

# Define the source UI files
UI_FILES    = $(wildcard *.ui)
GEN_FILES   = $(patsubst %.ui,%.py,$(UI_FILES))
XTRA_FILES  = $(wildcard *_ui.py)

# Detect whether to disable colors
ifeq ($(NO_COLOR)$(TERM),dumb)
COLOR_ENABLED := 0
else
COLOR_ENABLED := 1
endif

# ANSI bright + bold color codes
ifeq ($(COLOR_ENABLED),1)
BRIGHT_GREEN  := \033[1;92m
BRIGHT_CYAN   := \033[1;96m
BRIGHT_YELLOW := \033[1;93m
RESET         := \033[0m
else
BRIGHT_GREEN  :=
BRIGHT_CYAN   :=
BRIGHT_YELLOW :=
RESET         :=
endif

all: $(GEN_FILES) cleanup

%.py: %.ui
<<<<<<< HEAD
	@printf "$(BRIGHT_GREEN)[BUILD]$(RESET) Converting %s → %s\n" "$<" "$@"
	@mkdir -p $(@D)
	@pyuic6 -o $@ $<

cleanup:
	@printf "$(BRIGHT_CYAN)[CLEANUP]$(RESET) Removing extraneous *_ui.py files (if any)\n"
	@if [ -n "$(XTRA_FILES)" ]; then rm $(XTRA_FILES); fi
	@printf "$(BRIGHT_CYAN)[FORMAT]$(RESET) Running black formatter\n"
	@black .
=======
	@mkdir -p $(@D) # Create the output directory if it doesn't exist
	pyuic6 -o $@ $<
# Remove extra files created by Qt Designer 
	rm $(XTRA_FILES)
# Clean up formatting in the generated Python files
	black .
>>>>>>> eca9007 (Add new RIMAS slit offset tool)

clean:
	@printf "$(BRIGHT_YELLOW)[CLEAN]$(RESET) Removing generated .py files\n"
	@rm -f *.py
	@rmdir $(PY_DIR) 2>/dev/null || true
