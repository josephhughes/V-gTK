from AnnotateGFF import AnnotateGFF

reference_list = "generic/ref_list.txt"  # if needed, or leave as ""
tmp_dir = "tmp"
ref_gff_dir = "tmp/ref_gff/raw_gff/"
output_dir = "ref_gff/annotated"
annotate_gtf = False  # if you're planning on extending for GTF
master_gff = "tmp/Gff/NC_001542.gff3"
diff_window = 15  # or whatever value you use for window overlap checks

annotator = AnnotateGFF(
    reference_list,
    tmp_dir,
    ref_gff_dir,
    output_dir,
    annotate_gtf,
    master_gff,
    diff_window
)
annotator.process()
