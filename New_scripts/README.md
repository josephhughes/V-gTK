# Overview
New scripts with class and init method for tool merging

## Running the scripts
```shell
conda create --file vgtk_env.yml
conda activate vgtk
```

### 1. Download GenBank XML's
```shell
python GenBankFetcher.py -t 11520
```

### 2. GenBank XML to TSV
```shell
python GenBankParser.py
```

### 3. Validate metadata for Country, Host and Date
```shell
python ValidateMatrix.py
```

### 4. AddMissing data
```shell
python AddMissingData.py -f generic/fillup_table.tsv
```
#### Format for fillup table
| primary_accession | country | host | collection_date |
|----------|----------|----------|----------|
| PP706245   | Czechia   |  Dog  | 2023 |

####       OR 

```shell
python AddMissingData.py -b generic/bulk_fillup_table.tsv
```
#### Format for bulk table
| host                          | replaced_by |
|-------------------------------|-------------|
| Canis lupus familiaris brain  | Dog        |
| bovine                        | Cow     |

### 5. Filter and extract sequences
```shell
python FilterAndExtractSequences.py -g tmp/AddMissingData/gB_matrix_replaced.tsv -r generic-influenza/ref_list.txt -v Y
```

#### Example of ref list format for segmented virus
| CY005140|H9N6                          | 1 |
|-------------------------------|-------------|
| CY005371|H12N5  | 1        |
| CY067675|H7N9                        | 1     |
