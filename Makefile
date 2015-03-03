################################################################################
NAME=vacumm
CURDIRNAME=$(shell basename $(CURDIR))
PYTHON_PACKAGE_NAME=$(NAME)
PYTHON_VERSION=$(shell  python -c 'import sys; print "%d.%d"%sys.version_info[:2]')
VERSION=$(shell PYTHONPATH=lib/python python -c 'import '$(PYTHON_PACKAGE_NAME)'; print '$(PYTHON_PACKAGE_NAME)'.__version__')
RELEASE=$(shell PYTHONPATH=lib/python python -c 'import '$(PYTHON_PACKAGE_NAME)'; print '$(PYTHON_PACKAGE_NAME)'.__release__')
DATE=$(shell PYTHONPATH=lib/python python -c 'import '$(PYTHON_PACKAGE_NAME)'; print '$(PYTHON_PACKAGE_NAME)'.__date__')
SETUP_FILE_NAME=setup.py
SDIST_FILE_NAME=$(NAME)-$(VERSION).tar.gz
SDIST_FILE_PATH=dist/$(SDIST_FILE_NAME)
TEST_INST_DIR=$(CURDIR)/test/install
TEST_SETUP_INST_DIR=$(TEST_INST_DIR)/setup-inst
#RPM_NAME=$(NAME)-$(VERSION)-$(RELEASE)
RPM_NAME=$(NAME)-$(VERSION)-*
#RPM_FILE_NAME=$(RPM_NAME).noarch.rpm
RPM_FILE_NAME=$(RPM_NAME).*.rpm
RPM_FILE_PATH=dist/$(RPM_FILE_NAME)
TEST_RPM_INST_DIR=$(TEST_INST_DIR)/rpm-inst
#EXCLUDES=--exclude var --exclude *.pid --exclude *.log* --exclude *.out --exclude *.pyc
EXCLUDES=--exclude var

################################################################################
.PHONY: lib doc html pdf

################################################################################
all: help
help:
	@echo ""
	@echo "$(NAME) $(VERSION).$(RELEASE) ($(DATE))"
	@echo ""
	@echo "Usage: make [<target1> [<target2>...]]"
	@echo ""
	@echo "Available make targets:"
	@echo ""
	@echo "  Build distribution:"
	@echo "    sdist                     build python source distribution $(SDIST_FILE_PATH)"
	@echo "    rpm                       build rpm $(RPM_FILE_PATH)"
	@echo "    dist                      call all build targets"
	@echo ""
	@echo "  Test:"
	@echo "    test-install              test python setup in $(TEST_SETUP_INST_DIR)"
	@echo "    test-rpm-info             test rpm info of $(RPM_FILE_PATH)"
	@echo "    test-rpm-install          test rpm install in $(TEST_RPM_INST_DIR)"
	@echo "    test-rpm-update           test rpm update in $(TEST_RPM_INST_DIR)"
	@echo "    test-rpm-erase            test rpm uninstall in  $(TEST_RPM_INST_DIR)"
	@echo "    test                      call all test targets"
	@echo ""
	@echo "  Clean"
	@echo "    clean-doc                 clean documentations"
	@echo "    clean-lib                 clean source package"
	@echo "    clean-build               clean builds"
	@echo "    clean-test                clean test installations"
	@echo "    clean-all                 call all clean and uninstall targets"
	@echo ""
	@echo "  Development"
	@echo "    lib                       local compilation of modules and extensions"
	@echo "    doc                       generate documentations"
	@echo "    html                      generate html documentation"
	@echo "    pdf                       generate pdf documentation"
	@echo "    safedoc                   generate documentations and all their dependencies"
	@echo "    install                   prepare for local use (build libs, fix permissions)"
	@echo "    uninstall                 clean local installations"
	@echo "    arch                      clean all and create source archive in parent directory"
	@echo ""

################################################################################
# DVELOPMENT
################################################################################

doc:
	cd doc/sphinx && make

html:
	cd doc/sphinx && make html

pdf:
	cd doc/sphinx && make pdf

safedoc:
	cd test && make
	cd scripts/tutorials && make
	cd scripts/courses && make
	touch lib/python/vacumm/misc/color.py
	make doc

lib:
	python setup.py build_ext --inplace --force
	#cd lib/python/vacumm/misc/grid && make

# install: lib doc
install: lib

uninstall: clean-lib clean-doc

arch: clean-all
	cd .. && tar cvjf $(CURDIRNAME).tbz $(CURDIRNAME) $(EXCLUDES)

################################################################################
# BUILD DIST
################################################################################

sdist: clean-partial install
	python $(SETUP_FILE_NAME) sdist

rpm: clean-partial install
	python $(SETUP_FILE_NAME) bdist_rpm

dist: sdist rpm

################################################################################
# TEST
################################################################################

test-install: clean-build clean-test-install
	python $(SETUP_FILE_NAME) install -O1 --prefix=$(TEST_SETUP_INST_DIR)
	PYTHONPATH=$(TEST_SETUP_INST_DIR)/lib/python$(PYTHON_VERSION)/site-packages python -c "import "$(PYTHON_PACKAGE_NAME)"; "$(PYTHON_PACKAGE_NAME)".info()"

test-rpm-info:
	@echo -e "====\nInfo\n===="
	rpm -qip $(RPM_FILE_PATH)
	@echo -e "====\nList\n===="
	rpm -qlp $(RPM_FILE_PATH)
	@echo -e "========\nRequires\n========"
	rpm -qRp $(RPM_FILE_PATH)

test-rpm-install: rpm
	rpm --dbpath $(TEST_RPM_INST_DIR)/var/lib/rpm --initdb
	rpm --dbpath $(TEST_RPM_INST_DIR)/var/lib/rpm -ivh $(RPM_FILE_PATH) --nodeps --replacepkgs --prefix $(TEST_RPM_INST_DIR)/usr

test-rpm-update: rpm
	rpm --dbpath $(TEST_RPM_INST_DIR)/var/lib/rpm --initdb
	rpm --dbpath $(TEST_RPM_INST_DIR)/var/lib/rpm -Uvh $(RPM_FILE_PATH) --nodeps --replacepkgs --prefix $(TEST_RPM_INST_DIR)/usr

test-rpm-erase:
	rpm --dbpath $(TEST_RPM_INST_DIR)/var/lib/rpm -e $(RPM_NAME)

test: test-install test-rpm-install

################################################################################
# CLEAN
################################################################################

clean-doc:
	-cd doc/sphinx && make clean

clean-lib:
	-find lib/python -name '*.py[co]' -delete

clean-build:
	-rm -rf build MANIFEST MANIFEST.in setup.pyc

clean-dist:
	-rm -rf dist

clean-test-install:
	-rm -rf $(TEST_SETUP_INST_DIR)

clean-test-rpm-install:
	-rm -rf $(TEST_RPM_INST_DIR)

clean-test: clean-test-install clean-test-rpm-install
	-rmdir $(TEST_INST_DIR)

clean-partial: clean-build clean-lib

clean: clean-partial clean-test

clean-all: uninstall clean clean-dist


