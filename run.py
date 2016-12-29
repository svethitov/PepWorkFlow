#!/usr/bin/python
'''Workflow for Peptide clustering and alignment'''

import os
from pepwork.uniprotcollection import UniProtCollection

def main():
    '''Main function for the WorkFlow'''
    if os.path.isfile('records.bin'):
        mycollection = UniProtCollection('bin')
    elif os.path.isfile('records.list'):
        mycollection = UniProtCollection()
    else:
        print('ERROR: No file found exiting')
        return 1
    mycollection.save_records()
    mycollection.buildguidetree()
    mycollection.cluster()
    mycollection.find_motifs()
    mycollection.plot_3dscatter()
    mycollection.plot_dendrogram()
    mycollection.save_groups()
    print()
    print('Job Done.')


if __name__ == '__main__':
    main()
