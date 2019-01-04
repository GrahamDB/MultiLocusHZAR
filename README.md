# Simple Helper Scripts for automated processing of multiple loci and traits

## Dependencies

- hzar
- subplex
- doParallel (optional)

## Configuration

- Read and edit params.R for parameters related to the observed data 
(e.g. the range of distances sampled)

- Read and edit FIXME for parameters related to parallelism

- Read and edit FIXME for parameters related to model generation

## Usage notes

- To reduce the number of traces required, AICc is used to select the best models for each locus.
- This code will attempt to run loci in parallel.  Generation of traces for models with exponential tails can consume very large amounts of memory.  Be care when running these scripts on laptops.
- This code automatically makes several directories to contain the interim data files.  

## Execution

- Once you have everything setup, sourcing runHZAR.R should execute the entire analysis