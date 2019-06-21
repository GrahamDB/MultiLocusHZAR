# Simple Helper Scripts for automated processing of multiple loci and traits

## Dependencies

- hzar
- subplex
- doParallel (optional)

## Configuration

- Read and edit params.R for parameters related to the observed data 
(e.g. the range of distances sampled)

- Read and edit params.R for parameters related to parallelism

- Read and edit modelLoader.R for parameters related to model generation

## Usage notes

- To reduce the number of traces required, AICc is used to select the best models for each locus.
- This code will attempt to run loci in parallel.  Generation of traces for models with exponential tails can consume very large amounts of memory.  Be care when running these scripts on laptops.
- This code automatically makes several directories to contain the interim data files.  

## Execution

- Once you have everything setup, sourcing runHZAR.R should execute through model selection.
- When that script run to completion, adjust the two if(FALSE) code blocks to perform trace analysis as
needed on the relevant loci and traits

## References

The included example data is from a multi-locus, multi-trait analysis of Jacanas in the Jacana hybrid zone.

Lipshutz, S. E., Meier, J. I., Derryberry, G. E., Miller, M. J., Seehausen, O. and Derryberry, E. P. (2019), Differential introgression of a female competitive trait in a hybrid zone between sex‐role reversed species. Evolution, 73: 188-201. doi:10.1111/evo.13675