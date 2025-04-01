# Overview
The directory contains the reference accession list (generic/ref_list.txt) and the master GFF file (tmp/Gff).

## Running the workflow
```shell
bash vgtk-hcv.sh
```

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

#### Example of ref list file format for segmented virus
<table>
    <tbody>
        <tr>
            <td>CY005140|H9N6</td>
            <td>1</td>
        </tr>
        <tr>
            <td>CY005371|H12N5</td>
            <td>1</td>
        </tr>
        <tr>
            <td>CY067675|H7N9</td>
            <td>1</td>
        </tr>
        <tr>
            <td>CY075051|H9N2</td>
            <td>1</td>
        </tr>
        <tr>
            <td>CY079178|H3N8</td>
            <td>1</td>
        </tr>
        <tr>
            <td>CY096645|H8N4</td>
            <td>1</td>
        </tr>
    </tbody>
</table>



####      OR
```shell
python FilterAndExtractSequences.py -g tmp/AddMissingData/gB_matrix_replaced.tsv -r generic-influenza/ref_list.txt
```

#### Example of ref list file format for non segmented virus

<table>
    <tbody>
        <tr>
            <td>KF726852</td>
        </tr>
        <tr>
            <td>KF726853</td>
        </tr>
        <tr>
            <td>KM016899</td>
        </tr>
        <tr>
            <td>KJ004416</td>
        </tr>
        <tr>
            <td>KC193267</td>
        </tr>
        <tr>
            <td>AY138549</td>
        </tr>
    </tbody>
</table>

### 6. Blast alignment
Command for segmented viruses
```shell
python BlastAlignment.py -s Y -f generic-influenza/ref_list.txt
```

###         OR
For non segmented viruses such as RABV
```shell
python BlastAlignment.py
```

### 7. NextAlign alignment
```shell
python NextalignAlignment.py
```

### 8. Pad Alignment
```shell
python PadAlignment.py
```

### 9. Generate table
```shell
python GenerateTable.py
```

### 10. Create SQlite DB
```shell
python CreateSqliteDB.py
```
