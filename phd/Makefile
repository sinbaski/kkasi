ifeq ($(cygwin), 1)
	prefix := miktex-
endif


.PHONY: first
first:
	$(prefix)pdflatex MyThesis.tex

.PHONY: second
second:
	$(prefix)bibtex MyThesis.aux

