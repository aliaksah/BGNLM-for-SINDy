.PHONY: all install data manifests test smoke tune synth real analyze clean

R ?= Rscript
WORKERS ?= 4
THREADS ?= 2

all: install data manifests test smoke

install:
	$(R) scripts/install_dependencies.R

data:
	$(R) scripts/download_real_data.R

manifests:
	$(R) scripts/generate_manifest.R

test:
	$(R) -e 'testthat::test_dir("tests/testthat", stop_on_failure = TRUE)'

tune:
	$(R) scripts/run_tuning.R

smoke:
	$(R) scripts/generate_manifest.R --smoke_test
	bash scripts/run_remote_server.sh both $(WORKERS) $(THREADS)

synth:
	bash scripts/run_remote_server.sh synth $(WORKERS) $(THREADS)

real:
	bash scripts/run_remote_server.sh real $(WORKERS) $(THREADS)

analyze:
	$(R) scripts/aggregate_results.R
	$(R) scripts/run_mixed_models.R

clean:
	rm -rf results/checkpoints_synth results/checkpoints_real logs
