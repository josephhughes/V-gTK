import os
import re
import csv
import time
import sqlite3
import read_file
import subprocess
import pandas as pd
from os.path import join
from argparse import ArgumentParser
from collections import defaultdict

class CreateSqliteDB:
	def __init__(self, meta_data, features, pad_aln, gene_info, m49_countries, m49_interm_region, m49_regions, m49_sub_regions, proj_settings, tmp_dir, db_name):
		self.meta_data = meta_data
		self.features = features
		self.pad_aln = pad_aln
		self.gene_info = gene_info
		self.m49_countries = m49_countries
		self.m49_interm_region = m49_interm_region
		self.m49_regions = m49_regions
		self.m49_sub_regions = m49_sub_regions
		self.proj_settings = proj_settings
		self.tmp_dir = tmp_dir
		self.db_name = db_name

	def create_db(self):
		df_meta_data = pd.read_csv(join(self.meta_data), sep="\t")
		df_features = pd.read_csv(join(self.features), sep="\t")
		df_aln = pd.read_csv(join(self.pad_aln), sep="\t")
		df_gene = pd.read_csv(join(self.gene_info), sep="\t")
		df_m49_country = pd.read_csv(join(self.m49_countries), sep=",")
		df_m49_interm = pd.read_csv(join(self.m49_interm_region), sep=",")
		df_m49_region = pd.read_csv(join(self.m49_regions), sep=",")
		df_m49_sub_region = pd.read_csv(join(self.m49_sub_regions), sep=",")
		df_proj_setting = pd.read_csv(join(self.proj_settings), sep=",")

		conn = sqlite3.connect(self.db_name)
		cursor = conn.cursor()

		df_meta_data.to_sql("meta_data", conn, if_exists="replace", index=False)
		df_features.to_sql("features", conn, if_exists="replace", index=False)
		df_aln.to_sql("sequence", conn, if_exists="replace", index=False)
		df_gene.to_sql("genes", conn, if_exists="replace", index=False)
		df_m49_country.to_sql("m49_country", conn, if_exists="replace", index=False)
		df_m49_interm.to_sql("m49_intermediate", conn, if_exists="replace", index=False)
		df_m49_region.to_sql("m49_regions", conn, if_exists="replace", index=False)
		df_m49_sub_region.to_sql("m49_sub_regions", conn, if_exists="replace", index=False)
		df_proj_setting.to_sql("project_settings", conn, if_exists="replace", index=False)
		cursor.execute("PRAGMA foreign_keys = ON;")

		cursor.execute("""CREATE TABLE IF NOT EXISTS meta_data AS SELECT * FROM meta_data;""")
		cursor.execute("""CREATE TABLE IF NOT EXISTS features AS SELECT * FROM features;""")
		cursor.execute("""CREATE TABLE IF NOT EXISTS sequence AS SELECT * FROM sequence;""")
		cursor.execute("""CREATE TABLE IF NOT EXISTS genes AS SELECT * FROM genes;""")
		cursor.execute("""CREATE TABLE IF NOT EXISTS m49_country AS SELECT * FROM m49_country;""")
		cursor.execute("""CREATE TABLE IF NOT EXISTS m49_intermediate AS SELECT * FROM m49_intermediate;""")
		cursor.execute("""CREATE TABLE IF NOT EXISTS m49_regions AS SELECT * FROM m49_regions;""")
		cursor.execute("""CREATE TABLE IF NOT EXISTS m49_sub_regions AS SELECT * FROM m49_sub_regions;""")
		cursor.execute("""CREATE TABLE IF NOT EXISTS project_settings AS SELECT * FROM project_settings;""")

		conn.commit()
		conn.close()

def process(args):
	db_creator = CreateSqliteDB(
			args.meta_data,
			args.features,
			args.pad_aln,
			args.gene_info,
			args.m49_countries,
			args.m49_interm_region,
			args.m49_regions,
			args.m49_sub_regions,
			args.proj_settings,
			args.tmp_dir,
			args.db_name,
		)
	db_creator.create_db()


if __name__ == "__main__":
	parser = ArgumentParser(description='Creating sqlite DB')
	parser.add_argument('-m', '--meta_data', help='Meta data table', default="tmp/GenBank-matrix/gB_matrix.tsv")
	parser.add_argument('-t', '--tmp_dir', help='tmp directory to store all the db-ready tsv files', default="tmp/sqliteDB/")
	parser.add_argument('-f', '--features', help='Features table', default="tmp/Tables/features.tsv")
	parser.add_argument('-p', '--pad_aln', help='Padded alignment file', default="tmp/Tables/sequence.tsv")
	parser.add_argument('-g', '--gene_info', help='Gene table', default="generic/Tables/gene_info.csv")
	parser.add_argument('-mc', '--m49_countries', help='M49 countries', default="generic/Tables/m49_country.csv")
	parser.add_argument('-mir', '--m49_interm_region', help='M49 intermediate regions', default="generic/Tables/m49_intermediate_region.csv")
	parser.add_argument('-mr', '--m49_regions', help='M49 regions', default="generic/Tables/m49_region.csv")
	parser.add_argument('-msr', '--m49_sub_regions', help='M49 sub-regions', default="generic/Tables/m49_sub_region.csv")
	parser.add_argument('-s', '--proj_settings', help='Project settings', default="tmp/Software_info/software_info.tsv")
	parser.add_argument('-d', '--db_name', help='Name of the Sqlite database', default="rabv-gdb.db")
	args = parser.parse_args()

	process(args)

