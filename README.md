# ngs-sim
Synthetic next-generation sequencing data with artefacts

## Mixtures

Multiple fasta references and mutation profiles can be provided. A proportion argument specifies the concentration of each. If no proportion is specified the default is 1.

## Mutations

Mutations are defined in an ad hoc "variant call format" file. This is a list of:

```
position	deletions	insertions
```

Where `position` is the 0 indexed position in the reference sequence, `deletions` is the number of bases to remove, and `insertions` is the string of bases to insert.

## Simulating sequencing artefacts

Optionally, the sequencing error profile is defined as a transition matrix.

## Simulation

The simulation chains the following passes:

* mutation: specified variants are spiked into the reference
* amplification: the sequencing error profile is applied
* fragmentation: mutated strands are truncated to the mean and std fragment size. Pairs of mated reads are read from the ends. The second mate is reverse complemented.
