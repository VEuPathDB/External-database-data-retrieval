# External database retrieval
These scripts will retrieve various data from external databases, and where necessary reformat that data to suit VEuPathDB

<details open>
    <summary><b>Current databases supported:</b></summary>
    
- [PHI-base](http://www.phi-base.org/)
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
</details>

# 
## GO term annotations

<b>Data format that data loading need:</b>
    
```

```

<details>
 <summary>    
<b>Scripts:</b>
    </summary>
    
  </details>


#
## Phenotype data
<b>Data format that data loading need:</b>
    
```

```

<details>
 <summary>    
<b>Scripts:</b>
    </summary>


### `process-CGD-phenotypes.py` 
Downloads and processes CGD phenotypes

**Usage:**
```
python process-CGD-phenotypes.py
```    
  </details>


#
## Linkouts
**Data format that data loading need:**

A tab delimited list with two columns:

`Locus ID` is either the Protein ID or Gene Locus ID that can be found on both VEuPathDB _and_ the target database (e.g. PHI-base).
`External ID` is the ID in the target database that you want to link to

```
Locus ID	External ID
LINF_080010600	PHI:2643
LMJFC_180006800	PHI:8215
LMJFC_150005300	PHI:3459
LMJFC_350052700	PHI:3839
LMJFC_040015000	PHI:12114
...
```

<details>
 <summary>    
<b>Scripts:</b>
    </summary>

### `PHIbase-make-linkouts-walk-through.py` 
This will create lists of Locus IDs (Protein IDs or Gene Locus IDs) and their matching PHI-base ID (e.g. PHI:2643). 
This script currently requires user input for various file paths and has a step that must be done manually on the VEuPathDB site. 

**Usage:**
```
python PHIbase-make-linkouts-walk-through.py
```
  </details>
