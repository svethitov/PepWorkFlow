#!/usr/bin/python
'''Workflow for Peptide clustering and alignment'''

from pepwork.uniprotcollection import UniProtCollection

def main():
    '''Main function for the WorkFlow'''
    mycollection = UniProtCollection()
    mycollection.save_records()
    mycollection.buildguidetree()
    mycollection.cluster()
    mycollection.find_motifs()


if __name__ == '__main__':
    main()
