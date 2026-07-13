# PLoTM convenience targets. Override CONFIG / PY / RSCRIPT as needed.
CONFIG ?= config/RA.yaml
PY      ?= ./envs/cnmf/bin/python
RSCRIPT ?= ./envs/r_core/bin/Rscript
export PY RSCRIPT PYTHONNOUSERSITE=1

.PHONY: help envs all download posterior input cnmf stageB thresholds burden figures annotate analysis check
help:
	@echo "targets: envs | all | <phase> | check"
	@echo "  phases: download posterior input cnmf stageB thresholds burden figures annotate analysis"
	@echo "  set CONFIG=config/<disease>.yaml (default $(CONFIG))"

envs:
	conda env create -p ./envs/cnmf      -f environment/cnmf.yml
	conda env create -p ./envs/r_core    -f environment/r_core.yml
	conda env create -p ./envs/genebayes -f environment/genebayes.yml

all download posterior input cnmf stageB thresholds burden figures annotate analysis:
	drivers/run_disease.sh $(CONFIG) $@

# static checks (no data required)
check:
	@bash -n drivers/*.sh src/**/*.sh 2>/dev/null && echo "bash OK"
	@for f in $$(find src drivers -name '*.py'); do $(PY) -m py_compile $$f || exit 1; done && echo "py OK"
	@for f in $$(find src -name '*.R'); do $(RSCRIPT) -e "invisible(parse(file='$$f'))" >/dev/null || exit 1; done && echo "R OK"
