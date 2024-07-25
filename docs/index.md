# TING

TING is a toolkit for building virus genome databases, focusing on 
automating metadata curation and sequence alignment. It aims to 
facilitate the sharing and updating of multiple sequence alignments 
and metadata, making the process automated and reproducible. Validated 
metadata and sequence data are stored in a relational database, which 
can be queried via a web interface or programmatically using an API. 
The database uses a core schema from GLUE-tools, with modifications to 
how alignments are stored.


## Tools

### genBank_downloader.py

genBank_downloader.py enables the downloading of genbank xml formatted files for any 
organism by taxonomy identifier (taxid).

** Improvement of the script, will be an initial query to determine the number of sequences in the 
taxid and adjust the --max_ret automatically as users will want to download everything in the given taxid


### genBank_to_tsv.py

genBank_to_tsv.py extracts key metadata from the xml and convert it to a tab-separated value table (gB_matrix.tsv). 

### genBank_updater.py

Using the gB_matrix.tsv with the "Accession" and "Update Date" to compare to a new set of accessions numbers 
from NCBI. Anything that has a different date or any new accession is downloaded as genbank xml.


### validate_tsv.py

This script has different metadata validators. Currently the collection date, the collection country and the host are validated.
** add taxonomy hierarchy of the host based on the validation
