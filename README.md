# Offspring simulator

This script is for simulating synthethic offspring in Kerminen et al. 2020. More specifically, the script randomly samples parents from a given input list, models recombination between parental haplotypes and samples newly generated haplotypes for the offspring. For each generation , new offspring are simulated as long as there are unused parents left. The newly generated offspring are used as parents for the consecutive generation. 

## Inputs

```
file_name		Prefix for output file names
chr			Chromosome number to separate/parallel simulations.
recomb_file		Data frame with column "recom.rate" defining the recombination probability. See example file "Testi_recombination.txt".
input_phase_file	Parental haplotype file in a format used in ChromoPainter program (Lawson et al. 2012). See example file "Testi.phase".
input_label_file	List of parental IDs for defining their order in input_phase_file, similar to the ChromoPainter (Lawson et al. 2012) label file. See example file  "Test_label_infile.txt".
input_parent_file	List of defining parents to sample from. See example file "Testi.F0.txt".
```

## Outputs

```
.phase			ChromoPainter styled .phase file including simulated individuals.
_label_infile.txt	ID list file defining the order of simulated individuals in the .phase file. Prefix (F1) defines the generation and postfix (_1) defines the individual identifier.
_parents_FX.txt		Files listing the consecutive parental individuals.
```
