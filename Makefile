# Compilation :  make            (FLANN detecte via pkg-config, sinon recherche exhaustive)
#                make FLANN=0    (force la recherche exhaustive O(N^2), sans FLANN)
#                make METIS=1    (active la numerotation METIS_NodeND, option -num metis)
# Tests :        make test
# Experiences :  make experiences
CC      ?= gcc
CFLAGS  ?= -O2 -std=c11 -Wall -Wextra -pedantic
CPPFLAGS := -Iheaders
LDLIBS  := -lm
BUILD   := build

FLANN ?= $(shell pkg-config --exists flann 2>/dev/null && echo 1 || echo 0)
ifeq ($(FLANN),1)
  ADJ     := utils/adjacences.c
  CPPFLAGS += $(shell pkg-config --cflags flann 2>/dev/null)
  LDLIBS  += $(shell pkg-config --libs flann 2>/dev/null || echo -lflann)
else
  ADJ     := utils/adjacences_exhaustif.c
endif

ifeq ($(METIS),1)
  CPPFLAGS += -DAVEC_METIS
  LDLIBS  += -lmetis
endif

SRCS := $(ADJ) utils/normales.c utils/divergence.c utils/poisson.c utils/cholesky.c \
        utils/solver.c utils/sparse_matrix.c utils/numerotation.c utils/utils.c \
        utils/reconstruction.c
OBJS := $(patsubst %.c,$(BUILD)/%.o,$(SRCS))

TESTS := test_poisson test_symbolique test_cholesky_numerique test_solver \
         test_adjacences test_numerotation test_divergence test_reconstruction test_performances

all: reconstruction
ifneq ($(FLANN),1)
	@echo "Attention : FLANN introuvable, recherche exhaustive des voisins utilisee."
endif

reconstruction: $(BUILD)/src/main.o $(OBJS)
	$(CC) $(CFLAGS) -o $@ $^ $(LDLIBS)

$(BUILD)/%.o: %.c $(wildcard headers/*.h)
	@mkdir -p $(dir $@)
	$(CC) $(CPPFLAGS) $(CFLAGS) -c -o $@ $<

$(BUILD)/test_%: $(BUILD)/test/test_%.o $(OBJS)
	$(CC) $(CFLAGS) -o $@ $^ $(LDLIBS)

test: $(addprefix $(BUILD)/,$(TESTS))
	@for t in $(TESTS); do echo "===== $$t"; ./$(BUILD)/$$t || exit 1; done

experiences: reconstruction
	@mkdir -p $(BUILD)/mplconfig
	MPLBACKEND=Agg MPLCONFIGDIR=$(BUILD)/mplconfig PYTHONDONTWRITEBYTECODE=1 python3 test/experiences/lancer_experiences.py

clean:
	rm -rf $(BUILD) reconstruction test/resultat_experience scripts/__pycache__ test/experiences/__pycache__

.PHONY: all test experiences clean
.SECONDARY:
