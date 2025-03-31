#!/bin/bash

TAX_ID=${1:-3052230} # HCV

scripts_dir="$(dirname "$0")/scripts"
generic_dir="$(dirname "$0")/generic"

db_name="hcv-rgd-test.db"
skip_fill=${2:-true}  # Use this variable to control skipping AddMissingData.py
is_segmented=${3:-N} #segmented virus or not Y for Yes and N for Not
master_acc="AF009606"

#python GenBankFetcher.py
#python "${scripts_dir}/GenBankFetcher.py" -t "$TAX_ID -u tmp/GenBank-matrix/gB_matrix_raw.tsv"
if [ $? -ne 0 ]; then
  echo "Error: GenBankFetcher.py failed."
  exit 1
fi
echo "GenBankFetcher.py completed successfully."
echo ""

#python DownloadGFF.py
python "${scripts_dir}/DownloadGFF.py" -id $master_acc
if [ $? -ne 0 ]; then
  echo "Error: DownloadGFF.py failed."
  exit 1
fi
echo "DownloadGFF.py completed successfully."
echo ""

#python GenBankParser.py
python "${scripts_dir}/GenBankParser.py"
if [ $? -ne 0 ]; then
  echo "Error: GenBankParser.py failed."
  exit 1
fi
echo "GenBankParser.py completed successfully."
echo ""

# python ValidateMatrix.py
python "${scripts_dir}/ValidateMatrix.py"
if [ $? -ne 0 ]; then
  echo "Error: ValidateMatrix.py failed."
  exit 1
fi
echo "ValidateMatrix.py completed successfully."
echo ""

# python AddMissingData.py
if [ "$skip_fill" = false ]; then
  python "${scripts_dir}/AddMissingData.py" -b "${generic_dir}/bulk_fillup_table.tsv"
  if [ $? -ne 0 ]; then
    echo "Error: AddMissingData.py failed."
    exit 1
  fi
  echo "AddMissingData.py completed successfully."
else
  echo "Skipping AddMissingData.py as skip_fill is set to true."
fi
echo ""

if [ "$skip_fill" = true ]; then
  filter_matrix_path="tmp/GenBank-matrix/gB_matrix_raw.tsv"
else
  filter_matrix_path="tmp/AddMissingData/gB_matrix_replaced.tsv"
fi

#python FilterAndExtractSequences.py
if [ "$is_segmented" = "Y" ]; then
  python "${scripts_dir}/FilterAndExtractSequences.py" -g "$filter_matrix_path" -r "${generic_dir}/ref_list.txt" -v "$is_segmented"
else
  python "${scripts_dir}/FilterAndExtractSequences.py" -g "$filter_matrix_path" -r "${generic_dir}/ref_list.txt"
fi

if [ $? -ne 0 ]; then
  echo "Error: FilterAndExtractSequences.py failed."
  exit 1
fi
echo "FilterAndExtractSequences.py completed successfully."
echo ""

#python BlastAlignment.py
#if [ "$is_segmented" = "Y" ]; then
#	python "${scripts_dir}/BlastAlignment.py" -s Y -f "${generic_dir}/ref_list.txt" -m ${master_acc}
#else
#	python "${scripts_dir}/BlastAlignment.py" -f "${generic_dir}/ref_list.txt" -m ${master_acc}
#fi

if [ $? -ne 0 ]; then
  echo "Error: BlastAlignment.py failed."
  exit 1
fi
echo "BlastAlignment.py completed successfully."
echo ""

# python NextalignAlignment.py
python "${scripts_dir}/NextalignAlignment.py" -m $master_acc
if [ $? -ne 0 ]; then
  echo "Error: NextalignAlignment.py failed."
  exit 1
fi
echo "NextalignAlignment.py completed successfully."
echo ""

# python PadAlignment.py
python "${scripts_dir}/PadAlignment.py" -r "tmp/Nextalign/reference_aln/$master_acc"
if [ $? -ne 0 ]; then
  echo "Error: PadAlignment.py failed."
  exit 1
fi
echo "PadAlignment.py for query sequence is completed successfully."
echo ""

#python "${scripts_dir}/CalcAlignmentCord.py
python "${scripts_dir}/CalcAlignmentCord.py" -i "tmp/Pad-alignment/AF009606.aligned_merged_MSA.fasta" -m $master_acc -g "tmp/Gff/$master_acc.gff3"
if [ $? -ne 0 ]; then
  echo "Error: CalcAlignmentCord.py failed."
  exit 1
fi
echo "CalcAlignmentCord.py is completed successfully."
echo ""

# python SoftwareVersion.py
python "${scripts_dir}/SoftwareVersion.py"
if [ $? -ne 0 ]; then
  echo "Error: SoftwareVersion.py failed."
  exit 1
fi
echo "SoftwareVersion.py completed successfully."
echo ""


# python GenerateTables.py
python "${scripts_dir}/GenerateTables.py" -p "tmp/Pad-alignment/$master_acc.aligned_merged_MSA.fasta"
if [ $? -ne 0 ]; then
  echo "Error: GenerateTables.py failed."
  exit 1
fi
echo "GenerateTables.py completed successfully."
echo ""

# python CreateSqliteDB.py
python "${scripts_dir}/CreateSqliteDB.py" -d "${db_name}"
if [ $? -ne 0 ]; then
  echo "Error: CreateSqliteDB.py failed."
  exit 1
fi
echo "CreateSqliteDB.py completed successfully."
echo ""
