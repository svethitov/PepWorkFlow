'''Group Class'''
import subprocess
import pepwork.extract

class Group:
    '''Class holding expanded records list and there motifs'''
    def __init__(self, records, idx, kword, cluster_idx):
        '''Constructor for group object'''
        self.records = records
        self.idx = idx
        self.kword = kword
        self.cluster_idx = cluster_idx

        ## Creates MEME output
        print()
        print('Finding motifs for group {} ...'.format(idx))
        pepwork.extract.writefasta(self.records, 'working.fasta')
        # Purges redundant sequences
        print('Removing redundant sequences ...')
        command = 't_coffee -other_pg seq_reformat -in working.fasta \
                    -action +trim _seq_%%80'
        with open('trimmed.fasta', 'w') as trimmedfile:
            subprocess.run(command.split(), stdout=trimmedfile)
        # Runs MEME
        min_seq = min([len(self.records[key].sequence) for key in self.records.keys()])
        print('Running MEME ...')
        command = 'meme trimmed.fasta -mod zoops -nmotifs 8 -evt 0.0005\
                    -maxw {} -maxiter 1000 -o meme_{}_{}_{}'.\
                      format(min_seq, self.kword, self.idx, self.cluster_idx)
        subprocess.run(command.split())
