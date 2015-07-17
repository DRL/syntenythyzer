#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
File   		: syntenythyzer.py
Author 		: Dominik R. Laetsch, dominik.laetsch at gmail dot com 
To do 		:	1. read annotations (genes) from GFF3 into contig-objs
					- store contig objects in gene_to_contig dict
				2. read orthomcl clusters (only 1:1's of Gp and Gr ...)

"""

from __future__ import division
import sys 
import os

class contigObj():
	def __init__(self, contig_id):
		self.contig_id = contig_id
		self.genes = set()
		self.gene_list = []
		self.order_fw = {}
		self.order_rv = {}

	def add_gene(self, gene):
		self.genes.add(gene)
		self.gene_list.append(gene)

	def bless(self):
		self.order_fw = {x:i for i,x in enumerate(self.gene_list)}
		self.order_rv = {x:i for i,x in enumerate(self.gene_list[::-1])}

	def get_position_of(self, gene):
		return (self.order_fw[gene], self.order_rv[gene]) 

def parse_contigs_from_gff3(gffs):
	contigs = {}
	genes = set()
	gene2contigs = {}
	for gff in gffs:

		with open(gff) as fh:
			for line in fh:
				if line.startswith("#"):
					pass
				else:
					field = line.split("\t")
					if not field[2] == 'gene':
						pass
					else:
						contig_id = field[0]
						gene = field[8].rstrip("\n").lstrip("ID=").split(";")[0]
						if not contig_id in contigs:
							contigs[contig_id] = contigObj(contig_id)
						contigs[contig_id].add_gene(gene)
						genes.add(gene)

	for contig_id, obj in contigs.items():
		obj.bless()
		for gene in obj.genes:
			gene2contigs[gene] = obj
	return gene2contigs

def parse_clusters(cluster_file):

	clusters_on_contigs = {}
	with open(cluster_file) as fh:
		for line in fh:
			field = line.rstrip("\n").split("\t")
			clustername = field[0]
			proteinA = field[1].split("|")[1]
			proteinB = field[2].split("|")[1]
			print gene2contigs[proteinA].contig_id, gene2contigs[proteinA].get_position_of(proteinA)
			print gene2contigs[proteinB].contig_id, gene2contigs[proteinB].get_position_of(proteinB)


if __name__ == "__main__":
	
	__version__ = 0.1

	gffs = [sys.argv[1], sys.argv[2]]
	cluster_file = sys.argv[3]

	gene2contigs = parse_contigs_from_gff3(gffs)

	clusters = parse_clusters(cluster_file)
