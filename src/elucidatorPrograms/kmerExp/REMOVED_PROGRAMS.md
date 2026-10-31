# Removed k-mer programs

Removed on 2026-10-30 because they were stubs or had their bodies commented out. Every one is still in
elucidator commit `e226960`, so a starting point can be recovered with
`git show e226960:<path>`.

| program | runner | file at e226960 | what was there | intended to do |
|---|---|---|---|---|
| `kmerTestingGround` | kmerExp | `src/elucidatorPrograms/kmerExp/kmerExp_kmerTestingGround.cpp` (lines 165–729) | printed "Nothing in testing ground currently" and returned, followed by ~560 commented lines of an earlier `getKmerSharedBlocksBetweenGenomes` | scratch space for k-mer experiments; the commented code lives on as `getKmerSharedBlocksBetweenGenomes` |
| `kmerConnectionGraph` | kmerExp | `src/elucidatorPrograms/kmerExp/kmerExp_kmerConnectionGraph.cpp` (lines 20–66, with an empty `KmerSimpleConGraph` class) | parsed options (`--fastq1/--fastq2/--fastq/--fasta`, `--revCompMate`, `--kmerLength` default 4) and returned | build a k-mer connection graph from reads |
| `encodingSeqsByCommonKmers` | kmerExp | `src/elucidatorPrograms/kmerExp/encodingSeqsByCommonKmers.cpp` | parsed `--klen` (default 9) and `--outKey`, then a `@todo` | encode sequences by their common k-mers: take the most frequent k-mer, mask it, take the next most frequent, and repeat until no k-mer occurs more than once; write the reduction key |
| `compressPeptideKmersByReducedLibrary` | kmerExp | `src/elucidatorPrograms/kmerExp/compressPeptideKmersByScoring.cpp` (lines 73–231, inside `/* */`) | full body commented out, registration commented out | compress peptide k-mers (default length 16) with a reduced amino acid library, clustering with dbscan |
| `compressPeptideKmersByScoring` | kmerExp | same block as above | full body commented out, registration commented out | compress peptide k-mers (default length 16) by substitution-matrix score distance |
| `findUniqKmersFromGenomeSubRegions` | kmerExp | `src/elucidatorPrograms/kmerExp/kmerExp.cpp` line 97 | only a commented-out `addFunc`, no definition left | superseded by `kmerSetExp findUniqKmersFromGenomeSubRegionsMultiple` |
| `createMinimallyNonRedundantDownSampledSetAfterFragmentation` | kmerSetExp | `src/elucidatorPrograms/kmerSetExp/kmerSetExp.cpp` line 45, `kmerSetExp.hpp` line 41 | commented-out `addFunc` and declaration, no definition | a fragmenting version of `createMinimallyNonRedundantDownSampledSet` |
