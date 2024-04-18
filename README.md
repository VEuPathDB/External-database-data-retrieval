# External database retrieval
These scripts will retrieve various data from external databases, and where necessary reformat that data to suit VEuPathDB

**Current databases supported:**
- [PHIbase](http://www.phi-base.org/)
    - Linkouts
- [Pombase](https://www.pombase.org/)
    - GO annotations 
- [Flybase](http://flybase.org/)
    - GO annotations 
- [*Saccharomyces* Genome Database (SGD)](https://www.yeastgenome.org/)
    - GO annotations 
- [*Candida* Genome Database (CGD)](http://www.candidagenome.org/)
    - GO annotations
    - Phenotypes 

## GO term annotations
**Data format that data loading need:**
```

```

### Scripts

#
## Phenotype data
**Data format that data loading need:**
```

```

### Scripts

`process-CGD-phenotypes.py` - Downloads and processes CGD phenotypes

**Usage:**
```
python process-CGD-phenotypes.py
```

#
## Linkouts
**Data format that data loading need:**
A tab delimited list with two columns, `"Locus ID"` and `PHI MolConn ID`
```
Locus ID	PHI MolConn ID
LINF_080010600	PHI:2643
LMJFC_180006800	PHI:8215
LMJFC_150005300	PHI:3459
LMJFC_350052700	PHI:3839
LMJFC_040015000	PHI:12114
...
```

### Scripts

`PHIbase-make-linkouts-walk-through.py` - This will create lists of Locus IDs (Protein IDs or Gene Locus IDs) and their matching PHIbase ID (e.g. PHI:2643). 
This script currently requires user input for various file paths and has a step that must be done manually on the VEuPathDB site. 

**Usage:**
```
python PHIbase-make-linkouts-walk-through.py
```
