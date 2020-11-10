# Makefile for the AMOEBAE Example Workflow repository.

SHELL=/bin/bash

install:
	bash scripts/install.sh

pull:
	bash scripts/singularity_pull.sh

dry_run:
	bash scripts/dry_run.sh

run_data_setup:
	nohup bash scripts/run_data_setup.sh & echo $$! > pid_nohup.txt

run_fwd_only:
	nohup bash scripts/run_fwd_only.sh & echo $$! > pid_nohup.txt

run:
	nohup bash scripts/run_workflow.sh & echo $$! > pid_nohup.txt

archive:
	nohup bash scripts/archive_workflow.sh & echo $$! > pid_nohup.txt

killall:
	kill -9 `cat pid_nohup.txt`
	rm pid_nohup.txt
	killall -TERM snakemake
	echo Jobs may still be running on the cluster. Cancel those manually, if necessary.

clean:
	-rm nohup.out
	-rm pid_nohup.txt
	-rm slurm*

clean_results:
	rm -rf results/*

clean_virtualenvs:
	bash scripts/clean_virtualenvs.sh

unlock:
	bash scripts/unlock_workflow.sh

uninstall:
	bash scripts/uninstall.sh
	
.PHONY: list
list:
	@$(MAKE) -pRrq -f $(lastword $(MAKEFILE_LIST)) : 2>/dev/null | awk -v RS= -F: '/^# File/,/^# Finished Make data base/ {if ($$1 !~ "^[#.]") {print $$1}}' | sort | egrep -v -e '^[^[:alnum:]]' -e '^$@$$'


