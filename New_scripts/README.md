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

