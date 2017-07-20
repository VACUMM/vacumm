
uifiles=$(wildcard *.ui)
uinames=$(basename $(uifiles))
pyuifiles=$(addsuffix .py,$(uinames))
designer=$(shell which designer &>/dev/null && echo designer || echo designer-qt4)
pyuic=pyuic4

all: uic

uic: $(pyuifiles)

%.py: %.ui
	$(pyuic) $< -xo $@

clean:
	rm -f $(pyuifiles)

designer: _designer uic

_designer:
	$(designer) $(uifiles)
