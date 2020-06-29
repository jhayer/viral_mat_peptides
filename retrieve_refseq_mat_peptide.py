
#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse
from Bio.Seq import Seq
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import IUPAC


def main():
    parser = argparse.ArgumentParser(description="To retrieve gbk entries that\
     contain mat_peptide features from a gbk file containing several entries")
    parser.add_argument("-i", "--input_file",
        help="Genebank file to parse", required=True)
    parser.add_argument("-o", "--output_file",
        help="Genebank file containing only entries having mat_peptide \
        features", required=True)
    parser.add_argument("-f", "--fasta_mp",
        help="Ouput Fasta file containing all the mat_peptide sequences found \
        in the genbank entrie(s)", required=False)

    args = parser.parse_args()

    infile = args.input_file
    outfile = args.output_file
    fasta_o = args.fasta_mp

    filter_entries(infile, outfile, fasta_o)


def filter_entries(in_gbk, out_gbk, fasta_out):
    print(fasta_out)
    gb_handle = open(in_gbk, 'r')
    save_handle = open(out_gbk, 'w')
    fasta_handle = open(fasta_out, 'w')

    for gb_record in SeqIO.parse(gb_handle, "genbank"):
        new_rec = gb_record
        #print("Name %s, %i features" % (gb_record.name, len(gb_record.features)))

        product_mat_peptide = index_genbank_features(gb_record, "mat_peptide", "product")
        if product_mat_peptide:
        #    print(product_mat_peptide)
            SeqIO.write(new_rec, save_handle, "genbank")
            if fasta_out!=None:
                get_mat_peptide_translation(new_rec, fasta_handle)

    gb_handle.close()
    save_handle.close()
    fasta_handle.close()

def get_mat_peptide_translation(gb_record,fasta_handle):
    #fasta_handle = open(fasta_out, 'w')

    acc_num = gb_record.id
    taxid = 'none'

    for (index, feature) in enumerate(gb_record.features) :
        s_prod = "unknow_prot"
        mp_seq = ""
        # we retrieve the taxid from the source feature
        if feature.type == 'source':
            if 'db_xref' in feature.qualifiers:
                s_db_xref = feature.qualifiers['db_xref'][0]
                print(s_db_xref)
                if s_db_xref.startswith('taxon'):
                    taxid = s_db_xref.split(':')[1]

        elif feature.type == 'mat_peptide':
            # retrieve the mat_peptide seq from location and translate
            #print(feature.location)
            if feature.location.strand == +1:
                mp_seq = gb_record.seq[feature.location.start:feature.location.end].translate(table=1)
            else:
                mp_seq = gb_record.seq[feature.location.start:feature.location.end].reverse_complement().translate(table=1)

            if 'product' in feature.qualifiers:
                s_prod = feature.qualifiers['product'][0].replace(" ","-")

            key = "%s_%s_%s" % (acc_num, taxid, s_prod)

            #print(key)
            #print(mp_seq)
            mp_record = SeqRecord(Seq(str(mp_seq), IUPAC.protein),
                    id=key, name=acc_num+str(index), description="")
        #    mp_record.description(key)
        #    mp_record.seq(mp_seq)
            SeqIO.write(mp_record, fasta_handle, "fasta")

    #fasta_handle.close()


def index_genbank_features(gb_record, feature_type, qualifier) :
    answer = dict()
    for (index, feature) in enumerate(gb_record.features) :
        if feature.type==feature_type :
            if qualifier in feature.qualifiers :
                #There should only be one locus_tag per feature, but there
                #are usually several db_xref entries
                for value in feature.qualifiers[qualifier] :
                    if value in answer :
                        print("WARNING - Duplicate key %s for %s features %i and %i" \
                           % (value, feature_type, answer[value], index))
                    else :
                        answer[value] = index
    return answer


if __name__ == '__main__':
    main()
