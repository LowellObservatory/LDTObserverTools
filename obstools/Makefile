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
UI_FILES = $(wildcard UI/*.ui)

# Define the directory for generated Python files
PY_DIR = UI

# Default target: build all UI files to Python
all: $(patsubst $(PY_DIR)/%.ui,$(PY_DIR)/%.py,$(UI_FILES))

# Pattern rule for converting a single .ui file to .py
$(PY_DIR)/%.py: $(PY_DIR)/%.ui
	@mkdir -p $(@D) # Create the output directory if it doesn't exist
	pyuic6 -o $@ $<
# Remove extra files created by Qt Designer 
	rm $(PY_DIR)/*_ui.py
# Clean up formatting in the generated Python files
	black .

# Clean target: remove all generated Python files
clean:
	rm -f $(PY_DIR)/*.py
	@rmdir $(PY_DIR) 2>/dev/null || true # Remove the directory if it's empty
