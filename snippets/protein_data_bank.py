import requests
import io
import pandas as pd
from Bio import PDB
from Bio.SeqIO import PdbIO, FastaIO
# from Bio import ExPASy, SwissProt

# server = 'http://www.uniprot.org/uniprot'
#
#
# def do_request(server, ID='', **kwargs):
#     params = ''
#     req = requests.get('%s/%s%s' % (server, ID, params), params=kwargs)
#     if not req.ok:
#         req.raise_for_status()
#     return req
#
#
# req = do_request(server, query='gene:p53 AND reviewed:yes', format='tab',
#                  columns='id,entry name,length,organism,organism-id,database(PDB),database(HGNC)',
#                  limit='50')
#
#
# uniprot_list = pd.read_table(io.StringIO(req.text))
# uniprot_list.rename(columns={'Organism ID': 'ID'}, inplace=True)
# print(uniprot_list)

#
# p53_human = uniprot_list[uniprot_list.ID == 9606]['Entry'].tolist()[0]
# handle = ExPASy.get_sprot_raw(p53_human)
# sp_rec = SwissProt.read(handle)
#
# print(sp_rec.entry_name, sp_rec.sequence_length, sp_rec.gene_name)
# print(sp_rec.description)
# print(sp_rec.organism, sp_rec.seqinfo)
# print(sp_rec.sequence)
# print(sp_rec.comments)
# print(sp_rec.keywords)

repository = PDB.PDBList()
repository.retrieve_pdb_file('1TUP', pdir='.', file_format='pdb')
repository.retrieve_pdb_file('1OLG', pdir='.', file_format='pdb')
repository.retrieve_pdb_file('1YCQ', pdir='.', file_format='pdb')

parser = PDB.PDBParser()
p53_1tup = parser.get_structure('P 53 - DNA Binding', 'pdb1tup.ent')
p53_1olg = parser.get_structure('P 53 - Tetramerization', 'pdb1olg.ent')
p53_1ycq = parser.get_structure('P 53 - Transactivation', 'pdb1ycq.ent')


def print_pdb_headers(headers, indent=0):
    ind_text = ' ' * indent
    for header, content in headers.items():
        if type(content) == dict:
            print('\n%s%20s:' % (ind_text, header))
            print_pdb_headers(content, indent + 4)
            print()
        elif type(content) == list:
            print('%s%20s:' % (ind_text, header))
            for elem in content:
                print('%s%21s %s' % (ind_text, '->', elem))
        else:
            print('%s%20s: %s' % (ind_text, header, content))


print_pdb_headers(p53_1tup.header)

print(p53_1tup.header['compound'])
print(p53_1olg.header['compound'])
print(p53_1ycq.header['compound'])


def describe_model(name, pdb):
    print()
    for model in pdb:
        for chain in model:
            print('%s - Chain: %s. Number of residues: %d. Number of atoms: %d.' %
                  (name, chain.id, len(chain), len(list(chain.get_atoms()))))


describe_model('1TUP', p53_1tup)
describe_model('1OLG', p53_1olg)
describe_model('1YCQ', p53_1ycq)

#will go deep in a next recipe (bottom up)
for residue in p53_1tup.get_residues():
    if residue.id[0] in [' ', 'W']:
        continue
    print(residue.id)

res = next(p53_1tup[0]['A'].get_residues())
print(res)

for atom in res:
    print(atom, atom.serial_number, atom.element)

print(p53_1tup[0]['A'][94]['CA'])


def get_fasta(pdb_file, fasta_file, transfer_ids=None):
    fasta_writer = FastaIO.FastaWriter(fasta_file)
    fasta_writer.write_header()
    for rec in PdbIO.PdbSeqresIterator(pdb_file):
        if len(rec.seq) == 0:
            continue
        if transfer_ids is not None and rec.id not in transfer_ids:
            continue
        print(rec.id, rec.seq, len(rec.seq))
        fasta_writer.write_record(rec)


get_fasta(open('../results/pdb1tup.ent'), open('../results/1tup.fasta', 'w'),
          transfer_ids=['1TUP:B'])

get_fasta(open('../results/pdb1olg.ent'), open('../results/1olg.fasta', 'w'),
          transfer_ids=['1OLG:B'])

get_fasta(open('../results/pdb1ycq.ent'), open('../results/1ycq.fasta', 'w'),
          transfer_ids=['1YCQ:B'])
