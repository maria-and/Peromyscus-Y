# `known_y_proteins.fa`

The `braker_annotate` rule (`workflow/rules/annotate.smk`) expects a FASTA of known
mammalian Y-chromosome protein sequences at `config/known_y_proteins.fa` for homology
evidence. This isn't included in the repo (keeping it free of bulk reference data —
see `.gitignore`) — fetch it yourself, e.g. for the canonical Y genes used elsewhere
in this project (*Sry, Zfy1, Zfy2, Ddx3y, Uty, Eif2s3y, Usp9y, Kdm5d*):

```bash
# Example: pull mouse (Mus musculus) Y protein sequences from NCBI as a starting
# homology reference — works across rodents even without an exact species match.
for gene in Sry Zfy1 Zfy2 Ddx3y Uty Eif2s3y Usp9y Kdm5d; do
  esearch -db protein -query "${gene}[gene] AND Mus musculus[organism]" \
    | efetch -format fasta >> config/known_y_proteins.fa
done
```

Requires NCBI's `entrez-direct` tools (`esearch`/`efetch`), available via
`conda install -c bioconda entrez-direct`.
