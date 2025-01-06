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
AY138549
KC193267
KJ004416
KM016899
KF726853
KF726852
KC737850
GU358653
GU345748
GU345747
GU345746

