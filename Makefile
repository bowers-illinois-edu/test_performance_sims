.PHONY: all

all: Paper/paper.pdf \
     Analysis/siusims.rda \


## Using renv for libraries now although it is having trouble on the cluster with seeing
## libraries installed via mamba install in the cache and keeps wanting to reinstall everything from source.

## Paper tasks
Paper/paper.pdf: Paper/paper.tex \
	figures/fwer_strong_propblocks.pdf \
	figures/power_propblocks.pdf \
	figures/sim_fwer_weak.tex \
	figures/mdrc_results.tex
	cd Paper && latexmk -pdf paper.tex

## Figures and Tables and Result Interpretations
figures/fwer_propblocks.pdf: figures/simsres_plotting.rda \
	Analysis/simsres.rda \
	figures/fwer_plot.R
	Rscript figures/fwer_plot.R

figures/sim_fwer_weak.tex: figures/simsres_plotting.rda \
	Analysis/simsres.rda \
	figures/fwer_plot.R
	Rscript figures/fwer_plot.R

figures/overalloutcome.pdf: figures/descfakeoutcome.R \
	Data/mtwrkdat.rda Data/blocks.rda
	R  --file=figures/descfakeoutcome.R

figures/outcomesb.pdf: figures/overalloutcome.pdf
	touch figures/outcomesb.pdf

figures/descfakeoutcomes.pdf: figures/overalloutcome.pdf
	touch figures/descfakeoutcomes.pdf

figures/mdrc_results.tex: Analysis/create_mdrc_results.R \
	Analysis/utility_fns.R \
	Data/origdat_dt.rda
	R --file=Analysis/create_mdrc_results.R

## Analysis tasks

Analysis/simparms.rda: Analysis/simparms.R #renv.lock
	Rscript  Analysis/simparms.R

## Make the not_done_idx file out of all of the CSVS
Analysis/not_done_idx.rda: Analysis/what_has_been_done.R
	R --file=Analysis/what_has_been_done.R

Analysis/merged_results.csv:  Analysis/concat_clean.py
	source .venv/bin/activate && python Analysis/concat_clean.py

Analysis/simsres.rda: Analysis/merged_results.csv Analysis/make_rda.R
	R --file=Analysis/make_rda.R

figures/simsres_plotting.rda: figures/results_data_for_plots.R \
	Analysis/simsres.rda
	R --file=figures/results_data_for_plots.R

## Data tasks
### We are mostly just using a simulated dataset

make_test_data.done: Data/make_test_data.R
	R --file=Data/make_test_data.R

Data/idat_equal_nb.rda: make_test_data.done
Data/idat_not_equal_nb.rda: make_test_data.done
Data/bdat_equal_nb.rda: make_test_data.done
Data/bdat_not_equal_nb.rda: make_test_data.done

## Visualize the Makefile

makefile.plots: makefile.png makefile.pdf

## see also https://stackoverflow.com/questions/14784405/how-to-set-the-output-size-in-graphviz-for-the-dot-format/20536144
## dot -Tpdf -Gsize=9,15\! -Gdpi=100 -ofoo.pdf foo.gv

## A graph of the makefile
makefile.png: Makefile make_p_to_json.py json_to_dot.py
	make -qp | python3 make_p_to_json.py | python3 json_to_dot.py | dot -Tpng >| makefile.png

makefile.pdf: Makefile make_p_to_json.py json_to_dot.py
	make -qp | python3 make_p_to_json.py | python3 json_to_dot.py | dot -Tpdf >| makefile.pdf

## see also https://stackoverflow.com/questions/14784405/how-to-set-the-output-size-in-graphviz-for-the-dot-format/20536144
## dot -Tpng -Gsize=9,15\! -Gdpi=100 -ofoo.png foo.gv
