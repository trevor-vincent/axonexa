PYTHON3 := $(shell which python3 2>/dev/null)

PYTHON := python3
COVERAGE := --cov=pennylane_lightning_kokkos --cov-report term-missing --cov-report=html:coverage_html_report
TESTRUNNER := -m pytest tests --tb=short

CPP_DIR := ./src/

format: format-cpp

format-cpp:
ifdef check
	./bin/format --check --cfversion $(if $(version:-=),$(version),0) ./src
else
	./bin/format --cfversion $(if $(version:-=),$(version),0) ./src
endif
