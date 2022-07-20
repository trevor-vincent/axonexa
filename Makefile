PYTHON3 := $(shell which python3 2>/dev/null)

PYTHON := python3
TESTRUNNER := -m pytest tests --tb=short

CPP_DIR := ./src/

format: format-cpp

format-cpp:
ifdef check
	./bin/format --check --cfversion $(if $(version:-=),$(version),0) ./src
else
	./bin/format --cfversion $(if $(version:-=),$(version),0) ./src
endif
