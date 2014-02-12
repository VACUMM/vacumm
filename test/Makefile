# All tests in the right order
TESTS=test_cdat test_plot test_grid test_regrid test_dataset test_arakawa

# Verbose version
TESTSV=$(patsubst %,%/v,$(TESTS))

# Default target
all: $(TESTS)

# Verbose target
verbose: $(TESTV)

# Run a test
# - normal mode
$(TESTS):test_%:
	-python -m unittest $@
# - verbose mode
$(TESTV):test_%/v:
	-python -m unittest -v $(subst /v,,$@)

# Help
help:
	@echo "Run unittests"
	@echo ""
	@echo "    make             # run all tests"
	@echo "    make verbose     # same in verbose mode"
	@echo "    make test_NAME   # run a single test"
	@echo "    make test_NAME/v # same in verbose mode"
	@echo ""

