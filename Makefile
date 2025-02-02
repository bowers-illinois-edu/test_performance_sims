.PHONY: all

all: Paper/paper.pdf \
     Analysis/siusims.rda \

## Using renv for libraries now although it is having trouble on the cluster with seeing
## libraries installed via mamba install in the cache and keeps wanting to reinstall everything from source.

## Paper tasks
Paper/paper.pdf: Paper/paper.tex \
	figures/sim_fwer_weak.tex \
	cd Paper && latexmk -pdf paper.tex

## Figures and Tables and Result Interpretations
figures/sim_fwer_weak.tex: figures/simsres_plotting.rda \
	Analysis/simsres.rda \
	figures/fwer_plot.R
	Rscript figures/fwer_plot.R

figures/overalloutcome.pdf: figures/descfakeoutcome.R \
	Data/mtwrkdat.rda Data/blocks.rda
	R  --file=figures/descfakeoutcome.R

figures/descfakeoutcomes.pdf: figures/overalloutcome.pdf
	touch figures/descfakeoutcomes.pdf

## Analysis tasks
Analysis/simparms.rda: Analysis/simparms.R
	R --file=Analysis/simparms.R

#### The long simulations are running using Analysis/siusims.R
## it creates many CSV files, so it is hard to make a Makefile target for it.

## Make the not_done_idx file out of all of the CSVS
Analysis/not_done_idx.rda: Analysis/what_has_been_done.R
	R --file=Analysis/what_has_been_done.R

## Collect the CSV files together into an R dataset
## Use python for the merging of csv files mostly because I want to learn more about python
## the python environment exists for the libraries/packages needed by concat_clean.py
Analysis/merged_results.csv:  Analysis/concat_clean.py
	source .venv/bin/activate && python Analysis/concat_clean.py

Analysis/simsres.rda: Analysis/merged_results.csv Analysis/make_rda.R
	R --file=Analysis/make_rda.R

figures/simsres_plotting.rda: figures/results_data_for_plots.R \
	Analysis/simsres.rda
	R --file=figures/results_data_for_plots.R

### this next might take too long. Basically adding it here to help us understand relationships among files.
Analysis/p_sims_res_equal_nb_2024.rda: Analysis/siusims.R \
	Analysis/simparms.rda \
	Data/equal_blocksize_data.rda
	R --file=Analysis/siusims.R

## Data tasks
### Create the simulated base data one with equal block sizes and the other not equal block sizes

Data/make_test_data.done: Data/make_test_data.R
	R --file=Data/make_test_data.R

Data/idat_equal_nb.rda: Data/make_test_data.done
Data/idat_not_equal_nb.rda: Data/make_test_data.done
Data/bdat_equal_nb.rda: Data/make_test_data.done
Data/bdat_not_equal_nb.rda: Data/make_test_data.done

Data/blocks_sampled.rda: Data/make_test_data.done \
	Data/make_block_lists.R
	R --file=Data/make_block_lists.R

Data/block_datasets.done: Data/make_block_datasets.R \
	Data/blocks_sampled.rda \
	Data/make_test_data.done
	R --file=Data/make_block_datasets.R

Data/unequal_blocksize_data.rda: Data/block_datasets.done
Data/equal_blocksize_data.rda: Data/block_datasets.done

#################### The simple analysis approach

## This next shows a dependency but really it takes so long to run that we use the CSVS below.
## We need to figure out how to represent this process: run the sims.R files on
## a cluster and generate many CSVS
## then concatenate the CSVS into the merged_results.csv, convert the csv to
## an rda file, and then use that file for graphs, tables, etc.

Simple_Analysis/simple_raw_unadj_results.rda: Simple_Analysis/simple_unadj_sims.R \
	Simple_Analysis/simple_p_draws_fns.R
	R --file=Simple_Analysis/simple_unadj_sims.R

## In theory this next should depend on the csv files themselves but there are thousands.
## So, I just create a file called ready_to_merge.done by hand at the command line

Simple_Analysis/CSVS_Unadj/ready_to_merge_unadj_csvs.done:
	touch Simple_Analysis/CSVS_Unadj/ready_to_merge_unadj_csvs.done

Simple_Analysis/simple_sims_unadj_merged_results.csv: Simple_Analysis/CSVS_Unadj/ready_to_merge_unadj_csvs.done \
	Simple_Analysis/simple_concat_clean_unadj.py
	source .venv/bin/activate && python Simple_Analysis/simple_concat_clean_unadj.py

Simple_Analysis/simple_sims_unadj_results.rda: Simple_Analysis/simple_sims_unadj_merged_results.csv
	R --file=Simple_Analysis/simple_unadj_make_rda.R

Simple_Analysis/simple_results_exploration.done : Simple_Analysis/simple_results_exploration.R \
	Simple_Analysis/simple_sims_unadj_results.rda
	R --file=Simple_Analysis/simple_results_exploration.R


####################
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
