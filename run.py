'''Workflow for Peptide clustering and alignment'''

from cluster import *
from extract import *

def main():
    '''Main function for the WorkFlow'''
    print('The script expect list.list file:')
    swissrecords = getrecords('list.list')
    getseqstat(swissrecords)
    input('Press any key to continue ...')
    pdbswissrecords = getpdb(swissrecords)
    getseqstat(pdbswissrecords)
    input('Press any key to continue ...')
    savebinary('PDBRecordsSet.bin', pdbswissrecords)
    featurevector = getfeaturesvector(pdbswissrecords)
    

if __name__ == '__main__':
    main()
