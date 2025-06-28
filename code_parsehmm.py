import os 
import re
import sys




def read_hmmsearch_results(indir):
	hit_tab={}
	for f in os.listdir(indir):
		fin = open(indir+'/'+f,'r')
		gene=re.sub(r'.res',r'',f)
		for l in fin.readlines():
			if not re.search(r'^#',l.strip()):
				a=re.split(r'\s+',l.strip())
				acc,score=a[0],float(a[5])
				if not acc in hit_tab:
					hit_tab[acc]={gene:score}
				else:
					hit_tab[acc][gene]=score
	return hit_tab

def find_besthit(hit_tab):
	besthit_tab={}
	for acc in hit_tab:
		anno = max(hit_tab[acc], key=hit_tab[acc].get)
		besthit_tab[acc]=anno
	return besthit_tab

def write_anno(besthit_tab,outdir):
	rdsr_out = open('%s/00_rdsr_all.txt'%(outdir),'w')
	shdr_out = open('%s/00_shdr_all.txt'%(outdir),'w')
	sox_out = open('%s/00_sox_all.txt'%(outdir),'w')

	for acc, anno in besthit_tab.items():
		if anno in tab['rdsr']:
			anno = re.sub(r'oxidative_',r'',anno)
			rdsr_out.write('%s\t%s\n'%(acc,anno))
		elif anno in tab['shdr']:
			shdr_out.write('%s\t%s\n'%(acc,anno))
		elif anno in tab['sox']:
			sox_out.write('%s\t%s\n'%(acc,anno))





tab = {'rdsr':['oxidative_dsrA','oxidative_dsrB','dsrE','dsrF','dsrH'],
		'shdr':['sHdrA','sHdrB1','sHdrB2','sHdrC1','sHdrC2'],
		'sox':['soxA','soxB','soxC','soxD','soxX','soxY','soxZ']}

indir = sys.argv[1]#'01_hmmsearch_res'
outdir = sys.argv[2]#'02_gene_cluster'

hit_tab = read_hmmsearch_results(indir)
besthit_tab = find_besthit(hit_tab)
write_anno(besthit_tab,outdir)

		