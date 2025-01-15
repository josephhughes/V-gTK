The workflow works for both segmented and non-segmented virus.

Before running pleaase make sure that the reference list is created, refer [README](https://github.com/josephhughes/TING/tree/main/New_scripts) for better understanding of ref list file formats.

Settings for running the script, make the changes accordingly to the vgtk.sh script
```
TAX_ID=${1:-1980419} #segmented virus, hdv, hepatitis d virus

scripts_dir="$(dirname "$0")/scripts"
generic_dir="$(dirname "$0")/generic"
db_name="torno.db"
skip_fill=${2:-true}  # Use this variable to control skipping AddMissingData.py
is_segmented=${3:-Y} #segmented virus or not Y for Yes and N for Not
```
