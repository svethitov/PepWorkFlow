'''Group Class'''
import copy
import os
import math
import random
import subprocess
import time
from multiprocessing import pool
import numpy as np
import pandas as pd
from Bio import motifs
from Bio.motifs.jaspar import calculate_pseudocounts
from Bio import SeqIO, Seq
import pepwork.extract
from pepwork.plots import hist, plot_scattermatrix

def _distance_offset(work_chunk):
    '''Helper function for multithreding'''
    return work_chunk[0].dist_pearson(work_chunk[1])

def _run_meme(min_len, kword, idx, cluster_idx, scramble):
    '''Initinates meme motif search'''
    print('Running MEME ...')
    command = 'meme trimmed_{}.fasta -mod zoops -nmotifs 8 -evt 0.00005\
                -maxw {} -maxiter 1000 -o meme_{}_{}_{}_{}'.\
                    format(scramble, min_len, kword, idx, cluster_idx, scramble)
    subprocess.run(command.split())

def _read_motifs(kword, idx, cluster_idx, scramble):
    '''Reads motifs from MEME output'''
    os.chdir('meme_{}_{}_{}_{}'.format(kword, idx, cluster_idx, scramble))
    with open('meme.txt') as motifs_file:
        try:
            meme_motifs = motifs.parse(motifs_file, 'meme')
        except ValueError:
            meme_motifs = None
            print('No motifs found ...')
    os.chdir(os.pardir)
    return meme_motifs



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
        with open('trimmed_0.fasta', 'w') as trimmedfile:
            subprocess.run(command.split(), stdout=trimmedfile)
        # Read sequences that will be used from MEME
        self.trimmed_records = pepwork.extract.readfasta('trimmed_0.fasta')
        # Runs MEME
        min_len = min([len(self.records[key].sequence) for key in self.records.keys()])
        _run_meme(min_len, self.kword, self.idx, self.cluster_idx, 0)

        # Read motifs
        self.motifs = _read_motifs(self.kword, self.idx, self.cluster_idx, 0)

        ## Now Scramble sequences and run MEME on the scrambled ones
        self.scrambled = []
        self.motifs_scrambled = []
        for scrambling in range(0, 3):
            self.scrambled.append(copy.deepcopy(self.trimmed_records))
            for idx in range(0, len(self.scrambled[scrambling])):
                self.scrambled[scrambling][idx].seq = Seq.Seq(
                    ''.join(random.sample(
                        list(self.scrambled[scrambling][idx].seq),
                        len(self.scrambled[scrambling][idx].seq)
                    ))
                )


            SeqIO.write(
                self.scrambled[scrambling],
                'trimmed_{}.fasta'.format(scrambling + 1),
                'fasta'
            )
            _run_meme(min_len, self.kword, self.idx, self.cluster_idx, scrambling + 1)
            self.motifs_scrambled.append(
                _read_motifs(self.kword, self.idx, self.cluster_idx, scrambling + 1)
            )


class GroupCollection():
    '''Holds all groups from all UniprotCollections'''
    def __init__(self):
        '''Constructor for GroupCollection Class'''
        self.collections = pepwork.extract.loadbinary('collections.bin')
        self.background = {
            'A': 0.089352,
            'C': 0.012621,
            'D': 0.054292,
            'E': 0.061605,
            'F': 0.039367,
            'G': 0.071722,
            'H': 0.022238,
            'I': 0.057097,
            'K': 0.050686,
            'L': 0.098467,
            'M': 0.024041,
            'N': 0.039567,
            'P': 0.048783,
            'Q': 0.038465,
            'R': 0.056696,
            'S': 0.068216,
            'T': 0.055995,
            'V': 0.068216,
            'W': 0.013022,
            'Y': 0.029550
        }

        self.motifs = [] # list to hold all motifs

        # Make motifs data frame
        print('Populating dataframe ...')
        position = 0
        self.all_df = pd.DataFrame(columns=['KW', 'Group_idx', 'Motif_idx',\
            'position', 'kw_pos', 'ratio'])
        for kword in self.collections:
            kw_pos = 0
            for group in kword:
                if group.motifs is not None:
                    for idx in range(0, len(group.motifs)):
                        if group.motifs[idx] is not None:
                            ratio = len(group.motifs[idx].instances) / len(group.trimmed_records)
                            temp_df = pd.DataFrame(
                                [[group.kword, group.idx, idx, position, kw_pos, ratio]],
                                columns=\
                                ['KW', 'Group_idx', 'Motif_idx', 'position', 'kw_pos', 'ratio']
                                )
                            self.all_df = self.all_df.append(temp_df, ignore_index=True)
                            self.motifs.append(group.motifs[idx])
                            kw_pos += 2
                            position += 1

        self.all_df = self.all_df.astype({'Group_idx': int})
        self.all_df = self.all_df.astype({'Motif_idx': int})
        self.all_df = self.all_df.astype({'position': int})
        self.all_df = self.all_df.astype({'kw_pos': int})

        print('Setting background, pseudocount and alphabet on all motifs ...')
        for idx, motif in enumerate(self.motifs):
            print('Setting for motif {} {}'.\
            format(self.all_df.KW.iloc[idx], int(int(self.all_df.kw_pos.iloc[idx])/2)))
            motif.background = self.background
            motif.pseudocounts = calculate_pseudocounts(motif)

        self.distance_df = None
        self.offset_df = None
        self.stats = None


    def compare(self):
        '''Make pairwise comparison of all motifs'''

        index_list = []
        for row in range(0, len(self.all_df)):
            index_list.append('{}_{}_{}'.format(
                self.all_df.KW.iloc[row],
                self.all_df.Group_idx.iloc[row],
                self.all_df.Motif_idx.iloc[row]
            ))

        self.distance_df = pd.DataFrame(index=index_list, columns=index_list)
        self.offset_df = pd.DataFrame(index=index_list, columns=index_list)

        # Populates distance_df and offset_df
        work_list = []
        print('Starting calculation of distances and offset on {}'.format(time.ctime()))
        for idx in range(0, len(self.all_df)):
            for next_idx in range(idx, len(self.all_df)):
                row_label = '{}_{}_{}'.format(self.all_df.KW.iloc[idx], \
                    self.all_df.Group_idx.iloc[idx],\
                    self.all_df.Motif_idx.iloc[idx])
                column_label = '{}_{}_{}'.format(self.all_df.KW.iloc[next_idx], \
                    self.all_df.Group_idx.iloc[next_idx],\
                    self.all_df.Motif_idx.iloc[next_idx])
                if row_label != column_label:
                    #print('Calculating correlation for {} vs. {} ...'.\
                    #    format(row_label, column_label))
                    work_list.append((self.motifs[idx].pssm, self.motifs[next_idx].pssm,
                                      row_label, column_label))

        print('Working list populated on {}'.format(time.ctime()))

        with pool.Pool(30) as mypool:
            answers_list = mypool.map(_distance_offset, work_list)

        print('Calculation finished on {}'.format(time.ctime()))
        print('Now populating distance and offset dataframes ...')
        for idx in range(0, len(work_list)):
            self.distance_df.set_value(work_list[idx][2], work_list[idx][3], answers_list[idx][0])
            self.distance_df.set_value(work_list[idx][3], work_list[idx][2], answers_list[idx][0])
            self.offset_df.set_value(work_list[idx][2], work_list[idx][3], answers_list[idx][1])
            self.offset_df.set_value(work_list[idx][3], work_list[idx][2], (- answers_list[idx][1]))
        print('Calculation of distances and offset done on {}'.format(time.ctime()))

    def _write_links(self, handle, thickness, color, z, idx, second_idx, kw_dict):
        '''Writes a link to file'''
        handle.write('kw{} {} {} kw{} {} {} color={},z={},thickness={}\n'.format(
            kw_dict[self.all_df.KW.iloc[idx]],
            self.all_df.kw_pos.iloc[idx],
            self.all_df.kw_pos.iloc[idx] + 1,
            kw_dict[self.all_df.KW.iloc[second_idx]],
            self.all_df.kw_pos.iloc[second_idx],
            self.all_df.kw_pos.iloc[second_idx] + 1,
            color,
            z,
            thickness
            ))

    def _trim_outliers(self, data, thresh=4):
        '''Returns the same dimention vector with outliers trimmed
        This is done for better visual represetation'''
        data = np.array(data)
        outliers = self._is_outlier(data, thresh)
        inliers = []
        for idx in range(0, len(data)):
            if not outliers[idx]:
                inliers.append(data[idx])
        inliers = np.array(inliers)
        min_inliers = np.min(inliers)
        max_inliers = np.max(inliers)
        mean = np.mean(inliers)

        for idx in range(0, len(data)):
            if outliers[idx]:
                if data[idx] > mean:
                    print('Outlier found on position {} with value {}. Changed to {}'.\
                        format(idx, data[idx], max_inliers))
                    data[idx] = max_inliers
                elif data[idx] < mean:
                    print('Outlier found on position {} with value {}. Changed to {}'.\
                        format(idx, data[idx], min_inliers))
                    data[idx] = min_inliers

        return data


    def _is_outlier(self, points, thresh=4):
        '''Returns a boolean vector with outliers and inliers'''
        if len(points.shape) == 1:
            points = points[:, None]
        median = np.median(points, axis=0)
        diff = np.sum((points - median)**2, axis=-1)
        diff = np.sqrt(diff)
        med_abs_deviation = np.median(diff)
        modified_z_score = 0.6745 * diff / med_abs_deviation
        return modified_z_score > thresh

    def _write_data_line(self, handle, kw_dict, idx, value):
        '''Writes data line in file for circos'''
        handle.write('kw{} {} {} {}\n'.\
                    format(kw_dict[self.all_df.KW.iloc[idx]], self.all_df.kw_pos.iloc[idx],
                           self.all_df.kw_pos.iloc[idx] + 1, value))


    def write_circos_files(self):
        '''Writes data files ready for import in circos'''

        # init dataframe for stats of motifs
        self.stats = pd.DataFrame(
            index=range(0, len(self.motifs)),
            columns=['Length', 'Hits', 'E value', 'Ratio', 'Links', 'KWord'])

        self.stats['KWord'] = self.all_df['KW']
        # Writing 'chromosomes' file
        print('Writing "chromosomes" file ...')
        with open('chr.txt', 'w') as myfile:
            for idx, kword in enumerate(self.all_df.KW.unique()):
                myfile.write('chr - kw{} {} {} {} color_{}\n'.format(
                    idx,
                    kword,
                    0,
                    max(self.all_df.loc[lambda df: df.KW == kword, 'kw_pos']) + 1,
                    kword
                ))

        # adding 'bands' to 'chromosome' file
        print('Writing bands to file ...')
        with open('chr.txt', 'a') as myfile:
            for idx, kword in enumerate(self.all_df.KW.unique()):
                temp_df = self.all_df.loc[lambda df: df.KW == kword]
                start = 0
                color = 'lgrey_a1'
                for group in temp_df.Group_idx.unique():
                    stop = max(temp_df.loc[lambda df: df.Group_idx == group, 'kw_pos']) + 1
                    myfile.write('band kw{} group{} group{} {} {} {}\n'.format(
                        idx,
                        group,
                        group,
                        start,
                        stop,
                        color
                    ))
                    start = stop + 1

        kw_dict = {kw: idx for idx, kw in enumerate(self.all_df.KW.unique())}
        # writing eval file
        print('Writing evalue file ...')
        evalue = []
        for idx, motif in enumerate(self.motifs):
            if motif.evalue != 0:
                evalue.append(math.log10(motif.evalue))
            else:
                print('Found wrong evalue for motif {} {} {}'.\
                    format(self.all_df.KW.iloc[idx], self.all_df.Group_idx.iloc[idx],
                           motif.name))
                motif.evalue = float(input('Enter value manually: '))
                print(motif.evalue)
                if motif.evalue == 0:
                    print('Unable to set evalue')
                    motif.evalue = 1
                evalue.append(math.log10(motif.evalue))

        for idx, value in enumerate(evalue):
            if value == 0:
                evalue[idx] = min(evalue)

        # get all far outliers trimmed
        evalue = self._trim_outliers(evalue)

        evalue = [value/np.min(evalue) for value in evalue]

        with open('evalue.txt', 'w') as myfile:
            for idx, value in enumerate(evalue):
                self._write_data_line(myfile, kw_dict, idx, value)
                self.stats.set_value(idx, 'E value', value)

        # writing num_occurences file
        print('Writing num_occurences.txt file ...')
        num_occurences = [motif.num_occurrences for motif in self.motifs]
        num_occurences = self._trim_outliers(num_occurences)
        num_occurences = [value/np.max(num_occurences) for value in num_occurences]
        with open('num_occurences.txt', 'w') as myfile:
            for idx, value in enumerate(num_occurences):
                self._write_data_line(myfile, kw_dict, idx, value)
                self.stats.set_value(idx, 'Hits', value)

        # Writing lengths.txt
        print('Writing lengths.txt ...')
        lengths = [motif.length for motif in self.motifs]
        lengths = self._trim_outliers(lengths)
        lengths = [value/np.max(lengths) for value in lengths]
        with open('lengths.txt', 'w') as myfile:
            for idx, value in enumerate(lengths):
                self._write_data_line(myfile, kw_dict, idx, value)
                self.stats.set_value(idx, 'Length', value)

        # Writing text_labels.txt
        print('Writing text_labels.txt ...')
        with open('text_labels.txt', 'w') as myfile:
            for idx in range(0, len(self.motifs)):
                self._write_data_line(myfile, kw_dict, idx, idx)

        # Writing ratio.txt
        print('Writing ratio.txt ...')
        with open('ratio.txt', 'w') as myfile:
            for idx in range(0, len(self.motifs)):
                self._write_data_line(myfile, kw_dict, idx, self.all_df.iloc[idx].loc['ratio'])
                self.stats.set_value(idx, 'Ratio', self.all_df.iloc[idx].loc['ratio'])


        # Writing links.txt file
        links_dict = {key: 0 for key in range(0, len(self.motifs))}

        print('Writing links.txt file ...')
        with open('links.txt', 'w') as myfile:
            for idx in range(0, len(self.motifs)):
                for second_idx in range(idx, len(self.motifs)):
                    if self.distance_df.iloc[idx, second_idx] < 0.04:
                        self._write_links(
                            handle=myfile,
                            thickness=3,
                            color='rdylgn-9-div-1',
                            z=80,
                            idx=idx,
                            second_idx=second_idx,
                            kw_dict=kw_dict
                        )
                        links_dict[idx] += 1
                        links_dict[second_idx] += 1
                    elif self.distance_df.iloc[idx, second_idx] < 0.08:
                        self._write_links(
                            handle=myfile,
                            thickness=2,
                            color='rdylgn-9-div-2',
                            z=70,
                            idx=idx,
                            second_idx=second_idx,
                            kw_dict=kw_dict
                        )
                        links_dict[idx] += 1
                        links_dict[second_idx] += 1
                    elif self.distance_df.iloc[idx, second_idx] < 0.12:
                        self._write_links(
                            handle=myfile,
                            thickness=1,
                            color='rdylgn-9-div-3',
                            z=60,
                            idx=idx,
                            second_idx=second_idx,
                            kw_dict=kw_dict
                        )
                        links_dict[idx] += 1
                        links_dict[second_idx] += 1
                    elif self.distance_df.iloc[idx, second_idx] < 0.16:
                        self._write_links(
                            handle=myfile,
                            thickness=1,
                            color='rdylgn-9-div-4',
                            z=50,
                            idx=idx,
                            second_idx=second_idx,
                            kw_dict=kw_dict
                        )
                        links_dict[idx] += 1
                        links_dict[second_idx] += 1
                    elif self.distance_df.iloc[idx, second_idx] < 0.20:
                        self._write_links(
                            handle=myfile,
                            thickness=1,
                            color='rdylgn-9-div-5_a1',
                            z=40,
                            idx=idx,
                            second_idx=second_idx,
                            kw_dict=kw_dict
                        )
                        links_dict[idx] += 1
                        links_dict[second_idx] += 1
                    elif self.distance_df.iloc[idx, second_idx] < 0.24:
                        self._write_links(
                            handle=myfile,
                            thickness=1,
                            color='rdylgn-9-div-6_a2',
                            z=30,
                            idx=idx,
                            second_idx=second_idx,
                            kw_dict=kw_dict
                        )
                    elif self.distance_df.iloc[idx, second_idx] < 0.28:
                        self._write_links(
                            handle=myfile,
                            thickness=1,
                            color='rdylgn-9-div-7_a3',
                            z=20,
                            idx=idx,
                            second_idx=second_idx,
                            kw_dict=kw_dict
                        )
                    elif self.distance_df.iloc[idx, second_idx] < 0.32:
                        self._write_links(
                            handle=myfile,
                            thickness=1,
                            color='rdylgn-9-div-8_a4',
                            z=10,
                            idx=idx,
                            second_idx=second_idx,
                            kw_dict=kw_dict
                        )
                    elif self.distance_df.iloc[idx, second_idx] < 0.36:
                        self._write_links(
                            handle=myfile,
                            thickness=1,
                            color='rdylgn-9-div-9_a5',
                            z=0,
                            idx=idx,
                            second_idx=second_idx,
                            kw_dict=kw_dict
                        )

        # adds links_dict to self.stats
        for key in links_dict.keys():
            self.stats.set_value(key, 'Links', links_dict[key])

        ## Here start tiles file generation
        print('Writing tiles files ...')
        all_records = {} # dict holding all records
        for kword in self.collections:
            for group in kword:
                for key in group.records.keys():
                    all_records[key] = group.records[key]

        # Writes all used records to a file
        with open('used_records.txt', 'w') as myfile:
            for key in all_records.keys():
                myfile.write('{}\n'.format(key))

        all_keywords = ['Lantibiotic', 'Bacteriocin', 'Antibiotic',
                        'Bacteriolytic enzyme', 'Defensin', 'Fungicide']
        all_keywords = set(all_keywords)

        handles = {}
        for kword in all_keywords:
            handles[kword] = open('{}_tile.txt'.format(kword), 'w')

        for idx, motif in enumerate(self.motifs):
            motif_sequences_names = [instance.sequence_name for instance in motif.instances]
            motif_kwords = set()
            for sequence in motif_sequences_names:
                for kword in all_records[sequence].keywords:
                    motif_kwords.add(kword)
            motif_kwords = motif_kwords.intersection(all_keywords)
            motif_kwords = motif_kwords.difference(set([self.all_df.KW.iloc[idx]]))

            for motif_kword in motif_kwords:
                self._write_data_line(handles[motif_kword], kw_dict, idx, '')

        for handle in handles.values():
            handle.close()


    def plot_hist(self, mode='distance', filename='distance.html'):
        '''Plots histogram of distances between motifs'''
        values = []
        for idx in range(0, len(self.motifs)):
            for second_idx in range(idx, len(self.motifs)):
                if mode == 'distance':
                    values.append(self.distance_df.iloc[idx, second_idx])
                elif mode == 'offset':
                    values.append(self.offset_df.iloc[idx, second_idx])
        hist(values, filename)

    def plot_scattermatrix(self):
        '''Plots Scatterplot Matrix'''
        plot_scattermatrix(self.stats, 'KWord', 'ScatterMatrix.html')
