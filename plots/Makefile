NPAR = 3 4
XRANGE = 8

ifeq (,$(findstring u,$(MAKEFLAGS)))
INPUT = ~/work/bh_analysis2/H1j_angles.root
SUF = full
else
INPUT = ../data/H1j_mtop_unweighted.root
SUF = unw
endif

FIT = fits_$(SUF)$(if $1,_)$1$(if $2,_)$2.root
FITS = $(foreach n,$(NPAR),$(foreach r,$(XRANGE),$(call FIT,$n,$r)))
FITS_pdf = $(FITS:.root=.pdf)
PARS = $(FITS:fits%=pars%)
PARS_pdf = $(PARS:.root=.pdf)

LLR = $(patsubst fits%,llr%,$(foreach r,$(XRANGE),$(call FIT,,$r)))
LLR_pdf = $(LLR:.root=.pdf)

YRANGE_  = 0:2
YRANGE_5 = 0.35:0.8
YRANGE_8 = 0.25:1.25

define TOKENIZE
$(eval W = $(shell sed 's/\..\+//;s/_/ /g' <<< $(@F)))
$(eval n = $(word 3,$W))
$(eval r = $(word 4,$W))
endef

.PHONY: all

all: $(FITS_pdf) $(PARS_pdf) $(LLR_pdf)
	@echo $(LLR)

$(FITS): $(INPUT)
	$(call TOKENIZE)
	../bin/fit $< -o $@ -M 12:250:550 \
	  -p0:0:0$(shell (($n<4)) || echo ':1' ) \
	  -n$n -r0.$r --use-chi2-pars --nbins=50

$(FITS_pdf): %.pdf: %.root
	$(call TOKENIZE)
	../bin/draw $< -o $@ -y$(YRANGE_$r)

$(PARS): pars_%: fits_%
	../bin/pars $< -o $@

$(PARS_pdf): %.pdf: %.root
	$(call TOKENIZE)
	rxplot $< -o $@ -ltr \
	  -r 'gg/^.*' nl \
	  "tt/^.*/$n par fit,  cos**theta **in [-0.$r,0.$r)" \
	  'xx/^.*/Hj mass [GeV]'

$(LLR): llr_$(SUF)_%.root: $(foreach n,$(NPAR),$(call FIT,$n,%))
	../bin/llr $^ -o $@

$(LLR_pdf): %.pdf: %.root
	rxplot $< -o $@ -r 'sn/^P.*' 'xx/^.*/Hj mass [GeV]'

