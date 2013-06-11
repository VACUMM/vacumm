png = static/logo_vacumm.png static/modules.png
all: $(png)

static/%.png: static/%.pdf
	pdftoppm $< $(patsubst %.png,%,$@)
	convert $(patsubst %.png,%-1.ppm,$@) $@
	rm -rf $(patsubst %.png,%-*.ppm,$@)

static/%.pdf : %.tex
	pdflatex $< && echo pdflatex $<
	pdfcrop $(patsubst %.tex,%.pdf,$<) $@
#	rm -f $(patsubst %.tex,%.pdf,$<)

.SECONDARY: $(patsubst %.png,%.pdf,$(png))


