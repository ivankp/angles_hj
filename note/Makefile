PDF := note.pdf
TEX := pdflatex -interaction=batchmode -shell-escape

.PHONY: all clean $(CLEAN)

all: $(PDF)

# .SECONDEXPANSION:

$(PDF): %.pdf: %.tex ivanp.sty
	@$(TEX) $* > /dev/null && \
	$(TEX) $* > /dev/null || \
	printf '\033[0;31mCompilation failed\033[0m\n'
	@awk '/Warning/{a=1}/^$$/{a=0}a' $*.log | sed "s/.*Warning.*/\x1b[33m&\x1b[0m/"
	@awk '/^!/{a=1}/}$$/{print;a=0}a' $*.log | sed "s/^!.*/\x1b[31m&\x1b[0m/"

clean:
	@rm -fvr $(PDF) *.aux *.toc *.out *.log \
	  *.nav *.snm *.bbl *.blg *.bcf *.run.xml \
	  _minted-$(shell basename $(PDF) .pdf)

