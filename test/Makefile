# All tests in the right order
TESTS=test_cdat test_plot test_grid test_regrid test_stats test_diag test_arakawa test_mask test_io test_dataset test_ifroco test_color test_basemap

# Verbose version
TESTSV=$(patsubst %,%/v,$(TESTS))

# Other options
PYTHONOPTIONS=-W ignore

# Default target
all: $(TESTS)

# Verbose target
verbose: $(TESTV)

# Run a test
# - normal mode
$(TESTS):test_%:
	-python $(PYTHONOPTIONS) -m unittest $@
# - verbose mode
$(TESTV):test_%/v:
	-python $(PYTHONOPTIONS) -m unittest -v $(subst /v,,$@)

# Help
help:
	@echo "Run unittests"
	@echo ""
	@echo "    make             # run all tests"
	@echo "    make verbose     # same in verbose mode"
	@echo "    make test_NAME   # run a single test"
	@echo "    make test_NAME/v # same in verbose mode"
	@echo ""

