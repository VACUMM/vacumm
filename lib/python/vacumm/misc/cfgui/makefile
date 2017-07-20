
all: ui

test: ui
	./cfgui

ui:
	cd resources/ui && make

designer:
	cd resources/ui && make designer

clean:
	find . -name '*.py[co]' -delete
	cd resources/ui && make clean
