# TACCO examples: Repository for examples about how to use TACCO

This repository contains example notebooks on the usage of [TACCO](https://github.com/simonwm/tacco) as well as a snakemake workflow to prepare and provide all necessary files for their execution.

# How to use the examples

To view the examples with their results, open the executed notebooks in the `notebooks` directory of the `main` branch or in the examples section of the [TACCO documentation](https://simonwm.github.io/tacco/examples.html).

To execute the examples locally, clone the repository, and run the workflow with `snakemake`. The `main` branch contains also the executed notebooks, i.e. it is quite big, while the `devel` branch contains only the code itself to regenerate the notebooks. To clone only the devel branch use `git clone --single-branch --branch devel git@github.com:simonwm/tacco_examples.git`.

## Setup snakemake

In case you do not have `snakemake` set up already, you can follow the [snakemake instructions](https://snakemake.readthedocs.io/en/stable/getting_started/installation.html) or just install a new environment in your existing `conda` setup using `conda env create -f workflow/envs/snakemake_env.yml`
NB: While `snakemake` recommends to install `mamba` as drop in replacement for `conda` in the conda base environment, it does also work with regular `conda` by specifying the additional command line arguments `--use-conda --conda-frontend conda` to `snakemake`. It just used to take a while longer to get the individual environments set up with older conda versions, so we also recommend to use `mamba` or a [recent version of conda](https://conda.org/blog/2023-11-06-conda-23-10-0-release).

## Run examples

If you dont want to run the full set of examples, check the main `workflow/snakemake` file for options. You can

- call `snakemake` without any targets to prepare and run all examples,
- call it with the name of a selected example as first argument (e.g. `snakemake slideseq_mouse_olfactory_bulb`) to prepare the datasets and to run the notebook for this example only, or
- call it with the name of a selected example as first argument prepended with `prepare` (e.g. `snakemake prepare_slideseq_mouse_olfactory_bulb`) to prepare the datasets for this example only but not run the notebook.

For some examples (e.g. `slideseq_mouse_olfactory_bulb`) there is also a version using just a single sample (e.g. `slideseq_mouse_olfactory_bulb_single`) which can be used in place of the names for the full example. These run faster, as they have to download and process much less data. They also come in handy if your machine has not enough memory to run the full example.

### Running on the local node

To run the workflow on the machine on which the `snakemake` command is executed, use the `--cores <N>` option. This tells snakemake to plan with `N` cpu cores for working on all the steps in the workflow. Note that `--cores 1` will not limit TACCO to use all the cores which are available on the machine as this value is not propagated to the notebook. So changing this number will change only the degree of parallelization for the preparation of the datasets, i.e. the number of parallel download and data conversion tasks.

### Running somewhere else

`snakemake` supports a wide variety of distributed execution modes. To use them for this workflow consult the [snakemake documentation](https://snakemake.readthedocs.io/en/stable/executing/cli.html?highlight=profile#profiles).
