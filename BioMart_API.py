from pybiomart import Server
import pandas as pd
from datetime import date
import xlsxwriter
import openpyxl

# date for filename
date = date.today()

#excel = pd.read_excel("MANE.GRCh38.v1.0.summary.xlsx")
#refseq_ids = excel["RefSeq_nuc"].values.tolist()
refseq_ids = ['NM_004067','NM_001370581']

### access Ensembl server
server = Server('http://www.ensembl.org')

### access correct mart and dataset
mart = server['ENSEMBL_MART_ENSEMBL']
#ensembl = server.list_marts()

# dataset needed for query of biomart server
dataset = mart['hsapiens_gene_ensembl']

### find all attributes
ensembl_attributes = dataset.list_attributes().to_string()

### find all filters
ensembl_filters = dataset.list_filters().to_string()

chromosome_name = ['1', '2', '3', '4', '5', '6', '7', '8', '9', '10', '11', '12', '13', '14', '15', '16', '17', '18', '19', '20', '21', '22', 'X', 'Y', 'MT']

ext_attr = ['refseq_mrna', 'mim_morbid_accession', 'mim_morbid_description', 'mim_gene_accession','mim_gene_description',
            'hgnc_symbol', 'entrezgene_description', 'gene_biotype','transcript_count','transcript_length',
            'transcript_biotype', 'refseq_peptide','ccds','ucsc', 'external_transcript_name',
            'ensembl_transcript_id_version', 'ensembl_peptide_id_version', 'ens_lrg_transcript','hgnc_id',
            'entrezgene_id', 'ensembl_gene_id_version', 'uniprotswissprot', 'hpa_id','pdb','wikigene_id', 'ens_lrg_gene','chromosome_name',
            'start_position', 'end_position', 'strand', 'band', 'transcript_start', 'transcript_end', 'transcription_start_site']

ext_attr_header = ['RefSeq mRNA ID','MIM morbid accession', 'MIM morbid description', 'MIM gene accession', 'MIM gene description',
                   'HGNC symbol', 'NCBI gene (formerly Entrezgene) description', 'Gene type', 'Transcript count','Transcript length (including UTRs and CDS)', 'Transcript type',
                   'RefSeq peptide ID', 'CCDS ID', 'UCSC Stable ID', 'Transcript name', 'Transcript stable ID version', 'Protein stable ID version',
                   'LRG display in Ensembl transcript ID', 'HGNC ID', 'NCBI gene (formerly Entrezgene) ID', 'Gene stable ID version',
                   'UniProtKB/Swiss-Prot ID', 'Human Protein Atlas ID', 'PDB ID', 'WikiGene ID','LRG display in Ensembl gene ID', 'Chromosome/scaffold name',
                   'Gene start (bp)', 'Gene end (bp)', 'Strand', 'Karyotype band', 'Transcript start (bp)','Transcript end (bp)', 'Transcription start site (TSS)']

########## only needed for MANE select lists ##########
# function splitting RefSeq ID and version number
def NM_IDs(refseq_ids):
    NM_ids = []
    for i in range(0, len(refseq_ids)):
        # splitting version from transcript at '.'and only appending first element of split (transcript) to empty list
        NM_ids.append(refseq_ids[i].split('.')[0])
    return NM_ids

######################################################

# function to retrieve IDs from BioMart query
def BioMart(NM_ids, ext_attr, ext_attr_header):

    df_merge = []

    # sorting attributes based on refseq IDs from excel sheet
    for i in range(0, len(ext_attr)):
        # biomart query for each attribute filtered by chromosome name
        query = dataset.query(attributes=[ext_attr[i], 'refseq_mrna'], filters={'chromosome_name': chromosome_name})
        # filling NA with '.'
        df_attr = pd.DataFrame(query).fillna('.')
        # checking if refseq ID (excel table) is in query
        find_ext_attr = df_attr.loc[df_attr['RefSeq mRNA ID'].isin(NM_ids)]
        # grouping attributes by refseq ID and combining attributes for same ID
        df_merge.append(find_ext_attr.groupby('RefSeq mRNA ID')[ext_attr_header[i]].apply(lambda x: '|'.join(map(str, x))))
    return df_merge

# function to find CDS length
def find_CDS(ENST):
    # separate query for CDS length because it can only be accessed through ensembl transcript ID
    ensembl_transcript = pd.DataFrame(dataset.query(attributes=['cds_length', 'ensembl_transcript_id_version'],filters={'chromosome_name': chromosome_name})).fillna('.')
    # indexing cds length to keep correct order
    cds_length = ensembl_transcript.iloc[pd.Index(ensembl_transcript['Transcript stable ID version']).get_indexer(ENST)]
    return cds_length

# write output to excel
with pd.ExcelWriter(str(date)  + "_biomart_"+ "_output.xlsx") as writer:
    NM_ids = NM_IDs(refseq_ids)
    df_merge = BioMart(NM_ids, ext_attr, ext_attr_header)
    ENST = BioMart(NM_ids, ext_attr, ext_attr_header)[15].values.tolist()
    cds_length = find_CDS(ENST)

    # looping over output of different IDs
    for i in range(0, len(df_merge)):
        df = pd.DataFrame(df_merge[i])
        df = df.rename(columns={ext_attr_header[i]: ext_attr[i]})
        df.to_excel(writer,  engine="xlsxwriter", startcol=0+i, index=False)
    cds_length.to_excel(writer, columns=["CDS Length"], engine = "xlsxwriter", startcol=len(df_merge), index=False)
