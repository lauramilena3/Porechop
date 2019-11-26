"""
Copyright 2017 Ryan Wick (rrwick@gmail.com)
https://github.com/rrwick/Porechop

This module contains the class and sequences for known adapters used in Oxford Nanopore library
preparation kits.

This file is part of Porechop. Porechop is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by the Free Software Foundation,
either version 3 of the License, or (at your option) any later version. Porechop is distributed in
the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
details. You should have received a copy of the GNU General Public License along with Porechop. If
not, see <http://www.gnu.org/licenses/>.
"""


class Adapter(object):

    def __init__(self, name, start_sequence=None, end_sequence=None, both_ends_sequence=None):
        self.name = name
        self.start_sequence = start_sequence if start_sequence else []
        self.end_sequence = end_sequence if end_sequence else []
        if both_ends_sequence:
            self.start_sequence = both_ends_sequence
            self.end_sequence = both_ends_sequence
        self.best_start_score, self.best_end_score = 0.0, 0.0

    def best_start_or_end_score(self):
        return max(self.best_start_score, self.best_end_score)

    def is_barcode(self):
        return self.name.startswith('Barcode ')

    def barcode_direction(self):
        if '_rev' in self.start_sequence[0]:
            return 'reverse'
        else:
            return 'forward'

    def get_barcode_name(self):
        """
        Gets the barcode name for the output files. We want a concise name, so it looks at all
        options and chooses the shortest.
        """
        possible_names = [self.name, self.start_sequence[0]]
        if self.end_sequence:
            possible_names.append(self.end_sequence[0])
        barcode_name = sorted(possible_names, key=lambda x: len(x))[0]
        return barcode_name.replace(' ', '_')


# INSTRUCTIONS FOR ADDING CUSTOM ADAPTERS
# ---------------------------------------
# If you need Porechop to remove adapters that aren't included, you can add your own my modifying
# the ADAPTERS list below.
#
# Here is the format for a normal adapter:
#     Adapter('Adapter_set_name',
#             start_sequence=('Start_adapter_name', 'AAAACCCCGGGGTTTTAAAACCCCGGGGTTTT'),
#             end_sequence=('End_adapter_name', 'AACCGGTTAACCGGTTAACCGGTTAACCGGTT'))
#
# You can exclude start_sequence and end_sequence as appropriate.
#
# If you have custom Barcodes, make sure that the adapter set name starts with 'Barcode '. Also,
# remove the existing barcode sequences from this file to avoid conflicts:
#     Adapter('Barcode 1',
#             start_sequence=('Barcode_1_start', 'AAAAAAAACCCCCCCCGGGGGGGGTTTTTTTT'),
#             end_sequence=('Barcode_1_end', 'AAAAAAAACCCCCCCCGGGGGGGGTTTTTTTT')),
#     Adapter('Barcode 2',
#             start_sequence=('Barcode_2_start', 'TTTTTTTTGGGGGGGGCCCCCCCCAAAAAAAA'),
#             end_sequence=('Barcode_2_end', 'TTTTTTTTGGGGGGGGCCCCCCCCAAAAAAAA'))


ADAPTERS = [Adapter('SQK-NSK007',
                    start_sequence=('SQK-NSK007_Y_Top', 'AATGTACTTCGTTCAGTTACGTATTGCT'),
                    end_sequence=('SQK-NSK007_Y_Bottom', 'GCAATACGTAACTGAACGAAGT')),

            Adapter('Rapid',
                    start_sequence=('Rapid_adapter', 'GTTTTCGCATTTATCGTGAAACGCTTTCGCGTTTTTCGTGCGCCGCTTCA')),

            Adapter('SQK-MAP006',
                            start_sequence=('SQK-MAP006_Y_Top_SK63',    'GGTTGTTTCTGTTGGTGCTGATATTGCT'),
                            end_sequence=  ('SQK-MAP006_Y_Bottom_SK64', 'GCAATATCAGCACCAACAGAAA')),

            Adapter('SQK-MAP006 Short',
                    start_sequence=('SQK-MAP006_Short_Y_Top_LI32',    'CGGCGTCTGCTTGGGTGTTTAACCT'),
                    end_sequence=  ('SQK-MAP006_Short_Y_Bottom_LI33', 'GGTTAAACACCCAAGCAGACGCCG')),

            Adapter('PCR adapters 1',
                    start_sequence=('PCR_1_start', 'ACTTGCCTGTCGCTCTATCTTC'),
                    end_sequence=  ('PCR_1_end', 'GAAGATAGAGCGACAGGCAAGT')),

            Adapter('PCR tail 1',
                    start_sequence=('PCR_tail_1_start', 'TTAACCTTTCTGTTGGTGCTGATATTGC'),
                    end_sequence=  ('PCR_tail_1_end',   'GCAATATCAGCACCAACAGAAAGGTTAA')),

            Adapter('PCR tail 2',
                    start_sequence=('PCR_tail_2_start', 'TTAACCTACTTGCCTGTCGCTCTATCTTC'),
                    end_sequence=  ('PCR_tail_2_end',   'GAAGATAGAGCGACAGGCAAGTAGGTTAA')),


            # 1D^2 kit adapters are interesting. ONT provided the following sequences on their site:
            #   start: GGCGTCTGCTTGGGTGTTTAACCTTTTTGTCAGAGAGGTTCCAAGTCAGAGAGGTTCCT
            #   end:   GGAACCTCTCTCTGACTTGGAACCTCTCTGACAAAAAGGTTAAACACCCAAGCAGACGCCAGCAAT
            # But when looking at actual reads, I found two parts. The first corresponds to one end
            # of the provided sequences (through slightly different):
            Adapter('1D^2 part 1',
                    start_sequence=('1D2_part_1_start', 'GAGAGGTTCCAAGTCAGAGAGGTTCCT'),
                    end_sequence=  ('1D2_part_1_end',   'AGGAACCTCTCTGACTTGGAACCTCTC')),
            # and the second part corresponds to the other end, combined with a bit of standard 1D
            # adapter:
            Adapter('1D^2 part 2',
                    start_sequence=('1D2_part_2_start', 'CTTCGTTCAGTTACGTATTGCTGGCGTCTGCTT'),
                    end_sequence=  ('1D2_part_2_end',   'CACCCAAGCAGACGCCAGCAATACGTAACT')),
            # The middle part of the provided sequences is less common, so I've left it out of the
            # adapter sequences here.


            # Some barcoding kits (like the native barcodes) use the rev comp barcode at the start
            # of the read and the forward barcode at the end of the read.

            # Other barcoding kits (like the PCR and rapid barcodes) use the reverse barcode at the

            # start of the read and the rev comp barcode at the end of the read.

            Adapter('Barcode 1 (reverse)',

            start_sequence=('repBC01_rev', 'CACAAAGACACCGACAACTTTCTT'),

            end_sequence=('repBC01', 'AAGAAAGTTGTCGGTGTCTTTGTG')),

            Adapter('Barcode 2 (reverse)',

            start_sequence=('repBC02_rev', 'ACAGACGACTACAAACGGAATCGA'),

            end_sequence=('repBC02', 'TCGATTCCGTTTGTAGTCGTCTGT')),

            Adapter('Barcode 3 (reverse)',

            start_sequence=('repBC03_rev', 'CCTGGTAACTGGGACACAAGACTC'),

            end_sequence=('repBC03', 'GAGTCTTGTGTCCCAGTTACCAGG')),

            Adapter('Barcode 4 (reverse)',

            start_sequence=('repBC04_rev', 'TAGGGAAACACGATAGAATCCGAA'),

            end_sequence=('repBC04', 'TTCGGATTCTATCGTGTTTCCCTA')),

            Adapter('Barcode 5 (reverse)',

            start_sequence=('repBC05_rev', 'AAGGTTACACAAACCCTGGACAAG'),

            end_sequence=('repBC05', 'CTTGTCCAGGGTTTGTGTAACCTT')),

            Adapter('Barcode 6 (reverse)',

            start_sequence=('repBC06_rev', 'GACTACTTTCTGCCTTTGCGAGAA'),

            end_sequence=('repBC06', 'TTCTCGCAAAGGCAGAAAGTAGTC')),

            Adapter('Barcode 7 (reverse)',

            start_sequence=('repBC07_rev', 'AAGGATTCATTCCCACGGTAACAC'),

            end_sequence=('repBC07', 'GTGTTACCGTGGGAATGAATCCTT')),

            Adapter('Barcode 8 (reverse)',

            start_sequence=('repBC08_rev', 'ACGTAACTTGGTTTGTTCCCTGAA'),

            end_sequence=('repBC08', 'TTCAGGGAACAAACCAAGTTACGT')),

            Adapter('Barcode 9 (reverse)',

            start_sequence=('repBC09_rev', 'AACCAAGACTCGCTGTGCCTAGTT'),

            end_sequence=('repBC09', 'AACTAGGCACAGCGAGTCTTGGTT')),

            Adapter('Barcode 10 (reverse)',

            start_sequence=('repBC10_rev', 'GAGAGGACAAAGGTTTCAACGCTT'),

            end_sequence=('repBC10', 'AAGCGTTGAAACCTTTGTCCTCTC')),

            Adapter('Barcode 11 (reverse)',

            start_sequence=('repBC11_rev', 'TCCATTCCCTCCGATAGATGAAAC'),

            end_sequence=('repBC11', 'GTTTCATCTATCGGAGGGAATGGA')),

            Adapter('Barcode 12 (reverse)',

            start_sequence=('repBC12_rev', 'TCCGATTCTGCTTCTTTCTACCTG'),

            end_sequence=('repBC12', 'CAGGTAGAAAGAAGCAGAATCGGA')),

            Adapter('Barcode 13 (reverse)',

            start_sequence=('repBC13_rev', ' TCACACGAGTATGGAAGTCGTTCT'),

            end_sequence=('repBC13', 'AGAACGACTTCCATACTCGTGTGA')),

            Adapter('Barcode 14 (reverse)',

            start_sequence=('repBC14_rev', 'TCTATGGGTCCCAAGAGACTCGTT'),

            end_sequence=('repBC14', 'AACGAGTCTCTTGGGACCCATAGA')),

            Adapter('Barcode 15 (reverse)',

            start_sequence=('repBC15_rev', 'CAGTGGTGTTAGCGAGGTAGACCT'),

            end_sequence=('repBC15', 'AGGTCTACCTCGCTAACACCACTG')),

            Adapter('Barcode 16 (reverse)',

            start_sequence=('repBC16_rev', 'AGTACGAACCACTGTCAGTTGACG'),

            end_sequence=('repBC16', 'CGTCAACTGACAGTGGTTCGTACT')),

            Adapter('Barcode 17 (reverse)',

            start_sequence=('repBC17_rev', 'ATCAGAGGTACTTTCCTGGAGGGT'),

            end_sequence=('repBC17', 'ACCCTCCAGGAAAGTACCTCTGAT')),

            Adapter('Barcode 18 (reverse)',

            start_sequence=('repBC18_rev', 'GCCTATCTAGGTTGTTGGGTTTGG'),

            end_sequence=('repBC18', 'CCAAACCCAACAACCTAGATAGGC')),

            Adapter('Barcode 19 (reverse)',

            start_sequence=('repBC19_rev', 'ATCTCTTGACACTGCACGAGGAAC'),

            end_sequence=('repBC19', 'GTTCCTCGTGCAGTGTCAAGAGAT')),

            Adapter('Barcode 20 (reverse)',

            start_sequence=('repBC20_rev', 'ATGAGTTCTCGTAACAGGACGCAA'),

            end_sequence=('repBC20', 'TTGCGTCCTGTTACGAGAACTCAT')),

            Adapter('Barcode 21 (reverse)',

            start_sequence=('repBC21_rev', 'TAGAGAACGGACAATGAGAGGCTC'),

            end_sequence=('repBC21', 'GAGCCTCTCATTGTCCGTTCTCTA')),

            Adapter('Barcode 22 (reverse)',

            start_sequence=('repBC22_rev', 'CGTACTTTGATACATGGCAGTGGT'),

            end_sequence=('repBC22', 'ACCACTGCCATGTATCAAAGTACG')),

            Adapter('Barcode 23 (reverse)',

            start_sequence=('repBC23_rev', 'CGAGGAGGTTCACTGGGTAGTAAG'),

            end_sequence=('repBC23', 'CTTACTACCCAGTGAACCTCCTCG')),

            Adapter('Barcode 24 (reverse)',

            start_sequence=('repBC24_rev', 'CTAACCCATCATGCAGAACTATGC'),

            end_sequence=('repBC24', 'GCATAGTTCTGCATGATGGGTTAG')),

            Adapter('Barcode 25 (reverse)',

            start_sequence=('repBC25_rev', 'CATTGCGTTGCATACCCAACTTAC'),

            end_sequence=('repBC25', 'GTAAGTTGGGTATGCAACGCAATG')),

            Adapter('Barcode 26 (reverse)',

            start_sequence=('repBC26_rev', 'ATGAGAATGCGTAGTCGCTGTATG'),

            end_sequence=('repBC26', 'CATACAGCGACTACGCATTCTCAT')),

            Adapter('Barcode 27 (reverse)',

            start_sequence=('repBC27_rev', 'TGTAAGAGGTGAATCTAACCGTCG'),

            end_sequence=('repBC27', 'CGACGGTTAGATTCACCTCTTACA')),

            Adapter('Barcode 28 (reverse)',

            start_sequence=('repBC28_rev', 'GATACGGTGCCTTCTTAGGTTTCA'),

            end_sequence=('repBC28', 'TGAAACCTAAGAAGGCACCGTATC')),

            Adapter('Barcode 29 (reverse)',

            start_sequence=('repBC29_rev', 'GGTCTGTCAACCCAAGGTGTCTAG'),

            end_sequence=('repBC29', 'CTAGACACCTTGGGTTGACAGACC')),

            Adapter('Barcode 30 (reverse)',

            start_sequence=('repBC30_rev', 'TGGGTCGAAGTAGATCCTCACTGA'),

            end_sequence=('repBC30', 'TCAGTGAGGATCTACTTCGACCCA')),

            Adapter('Barcode 31 (reverse)',

            start_sequence=('repBC31_rev', 'CAATGTAACTGATTGCTGTACGCA'),

            end_sequence=('repBC31', 'TGCGTACAGCAATCAGTTACATTG')),

            Adapter('Barcode 32 (reverse)',

            start_sequence=('repBC32_rev', 'ATGACGTTGTCGGACTTCTACTGG'),

            end_sequence=('repBC32', 'CCAGTAGAAGTCCGACAACGTCAT')),

            Adapter('Barcode 33 (reverse)',

            start_sequence=('repBC33_rev', 'AGTTACCCAACCGTACCAAGTCTG'),

            end_sequence=('repBC33', 'CAGACTTGGTACGGTTGGGTAACT')),

            Adapter('Barcode 34 (reverse)',

            start_sequence=('repBC34_rev', 'GCCTTTGACTTGAGTTCTTCGTCC'),

            end_sequence=('repBC34', 'GGACGAAGAACTCAAGTCAAAGGC')),

            Adapter('Barcode 35 (reverse)',

            start_sequence=('repBC35_rev', 'GCAGTCCCTCAGCTTCGTAAGTAG'),

            end_sequence=('repBC35', 'CTACTTACGAAGCTGAGGGACTGC')),

            Adapter('Barcode 36 (reverse)',

            start_sequence=('repBC36_rev', 'TGTTTCCTCCTCTAACTGGGACAT'),

            end_sequence=('repBC36', 'ATGTCCCAGTTAGAGGAGGAAACA')),

            Adapter('Barcode 37 (reverse)',

            start_sequence=('repBC37_rev', 'TGATACTAAGCATCAATCGCAAGC'),

            end_sequence=('repBC37', 'GCTTGCGATTGATGCTTAGTATCA')),

            Adapter('Barcode 38 (reverse)',

            start_sequence=('repBC38_rev', 'TTCTCTGTATCGTCCTCCTGTGGT'),

            end_sequence=('repBC38', 'ACCACAGGAGGACGATACAGAGAA')),

            Adapter('Barcode 39 (reverse)',

            start_sequence=('repBC39_rev', 'GAGAGGCTCTAGTTGACACTGTGG'),

            end_sequence=('repBC39', 'CCACAGTGTCAACTAGAGCCTCTC')),

            Adapter('Barcode 40 (reverse)',

            start_sequence=('repBC40_rev', 'GGCTATCCTTGGTCATCCAAACTA'),

            end_sequence=('repBC40', 'TAGTTTGGATGACCAAGGATAGCC')),

            Adapter('Barcode 41 (reverse)',

            start_sequence=('repBC41_rev', 'CGTGTACTTCTCTGGACGAACTCC'),

            end_sequence=('repBC41', 'GGAGTTCGTCCAGAGAAGTACACG')),

            Adapter('Barcode 42 (reverse)',

            start_sequence=('repBC42_rev', 'CTGGCAGGTATGCCTTACACGTAG'),

            end_sequence=('repBC42', 'CTACGTGTAAGGCATACCTGCCAG')),

            Adapter('Barcode 43 (reverse)',

            start_sequence=('repBC43_rev', 'CTACCGTCGAGTCAACAACGAAAG'),

            end_sequence=('repBC43', 'CTTTCGTTGTTGACTCGACGGTAG')),

            Adapter('Barcode 44 (reverse)',

            start_sequence=('repBC44_rev', 'GAGTGGGAAGGAACCCTTTCTACT'),

            end_sequence=('repBC44', 'AGTAGAAAGGGTTCCTTCCCACTC')),

            Adapter('Barcode 45 (reverse)',

            start_sequence=('repBC45_rev', 'CACTGAAGGCATCTCTGTTGGATC'),

            end_sequence=('repBC45', ' GATCCAACAGAGATGCCTTCAGTG ')),

            Adapter('Barcode 46 (reverse)',

            start_sequence=('repBC46_rev', 'CAGGAGAATGAAGTGGAACACAGC'),

            end_sequence=('repBC46', 'GCTGTGTTCCACTTCATTCTCCTG')),

            Adapter('Barcode 47 (reverse)',

            start_sequence=('repBC47_rev', 'GAACTACCTGTGGGAAAGTTGCAC'),

            end_sequence=('repBC47', 'GTGCAACTTTCCCACAGGTAGTTC')),

            Adapter('Barcode 48 (reverse)',

            start_sequence=('repBC48_rev', 'TACAGGTGTACCACGTTCCAGATG'),

            end_sequence=('repBC48', 'CATCTGGAACGTGGTACACCTGTA')),

            Adapter('Barcode 49 (reverse)',

            start_sequence=('repBC49_rev', 'CTAGATGTTCAAAGCTGCACCAGT'),

            end_sequence=('repBC49', 'ACTGGTGCAGCTTTGAACATCTAG')),

            Adapter('Barcode 50 (reverse)',

            start_sequence=('repBC50_rev', 'ACGCAGGAAGTTACCAAAGTCCAT'),

            end_sequence=('repBC50', 'ATGGACTTTGGTAACTTCCTGCGT')),

            Adapter('Barcode 51 (reverse)',

            start_sequence=('repBC51_rev', 'GAGGACCCAGTAGGCTCATTCAAC'),

            end_sequence=('repBC51', 'GTTGAATGAGCCTACTGGGTCCTC')),

            Adapter('Barcode 52 (reverse)',

            start_sequence=('repBC52_rev', 'GTCCACGAACAATCTTGTCTCTCA'),

            end_sequence=('repBC52', 'TGAGAGACAAGATTGTTCGTGGAC')),

            Adapter('Barcode 53 (reverse)',

            start_sequence=('repBC53_rev', 'CTTTGCATGAGACGGTCTGAATCT'),

            end_sequence=('repBC53', 'AGATTCAGACCGTCTCATGCAAAG')),

            Adapter('Barcode 54 (reverse)',

            start_sequence=('repBC54_rev', 'CATGCTCCTTAGTCAAAGCTCTTG'),

            end_sequence=('repBC54', 'CAAGAGCTTTGACTAAGGAGCATG')),

            Adapter('Barcode 55 (reverse)',

            start_sequence=('repBC55_rev', 'CGTAGATCAGGGTCTCATCTTCCA'),

            end_sequence=('repBC55', 'TGGAAGATGAGACCCTGATCTACG')),

            Adapter('Barcode 56 (reverse)',

            start_sequence=('repBC56_rev', 'TTCATGCCACCTGTTGAGTAGTGA'),

            end_sequence=('repBC56', 'TCACTACTCAACAGGTGGCATGAA')),

            Adapter('Barcode 57 (reverse)',

            start_sequence=('repBC57_rev', 'ACTTCCGAAGGAGATTGACCTAGC'),

            end_sequence=('repBC57', 'GCTAGGTCAATCTCCTTCGGAAGT')),

            Adapter('Barcode 58 (reverse)',

            start_sequence=('repBC58_rev', 'TCAGACTCACGGAGGAGTAACCTG'),

            end_sequence=('repBC58', 'CAGGTTACTCCTCCGTGAGTCTGA')),

            Adapter('Barcode 59 (reverse)',

            start_sequence=('repBC59_rev', 'ACCTTGCTTTCCCTTCTTGATTGA'),

            end_sequence=('repBC59', 'TCAATCAAGAAGGGAAAGCAAGGT')),

            Adapter('Barcode 60 (reverse)',

            start_sequence=('repBC60_rev', 'CCATAGAAGCCTTGGTTGAACATG'),

            end_sequence=('repBC60', 'CATGTTCAACCAAGGCTTCTATGG')),

            Adapter('Barcode 61 (reverse)',

            start_sequence=('repBC61_rev', 'GTGCTGAGGCACATAGTACCCTCT'),

            end_sequence=('repBC61', 'AGAGGGTACTATGTGCCTCAGCAC')),

            Adapter('Barcode 62 (reverse)',

            start_sequence=('repBC62_rev', 'TACGTCCTGAAGTAAGTGTGGGTG'),

            end_sequence=('repBC62', 'CACCCACACTTACTTCAGGACGTA')),

            Adapter('Barcode 63 (reverse)',

            start_sequence=('repBC63_rev', 'GTTCAAGACCCAGGAACTTCAGAA'),

            end_sequence=('repBC63', 'TTCTGAAGTTCCTGGGTCTTGAAC')),

            Adapter('Barcode 64 (reverse)',

            start_sequence=('repBC64_rev', 'GAAAGTCGATGAACGGTGTCTGTC'),

            end_sequence=('repBC64', 'GACAGACACCGTTCATCGACTTTC')),

            Adapter('Barcode 65 (reverse)',

            start_sequence=('repBC65_rev', 'CCTTGTCTGGAGGAAGACTGAGAA'),

            end_sequence=('repBC65', 'TTCTCAGTCTTCCTCCAGACAAGG')),

            Adapter('Barcode 66 (reverse)',

            start_sequence=('repBC66_rev', 'GAAGTTAGAAGCCACAAGGATCGG'),

            end_sequence=('repBC66', 'CCGATCCTTGTGGCTTCTAACTTC')),

            Adapter('Barcode 67 (reverse)',

            start_sequence=('repBC67_rev', 'GGTGAGCACACGAGTATGACAAAC'),

            end_sequence=('repBC67', 'GTTTGTCATACTCGTGTGCTCACC')),

            Adapter('Barcode 68 (reverse)',

            start_sequence=('repBC68_rev', 'CCACCTTCGTGTTTGCTTAGATTC'),

            end_sequence=('repBC68', 'GAATCTAAGCAAACACGAAGGTGG')),

            Adapter('Barcode 69 (reverse)',

            start_sequence=('repBC69_rev', 'AGATCACATGAGGCTCGGACTGTA'),

            end_sequence=('repBC69', 'TACAGTCCGAGCCTCATGTGATCT')),

            Adapter('Barcode 70 (reverse)',

            start_sequence=('repBC70_rev', 'ACACTCCATTCGTAGGATCTCGGT'),

            end_sequence=('repBC70', 'ACCGAGATCCTACGAATGGAGTGT')),

            Adapter('Barcode 71 (reverse)',

            start_sequence=('repBC71_rev', 'CTGTTACTACCTGATGCTCCCAGG'),

            end_sequence=('repBC71', 'CCTGGGAGCATCAGGTAGTAACAG')),

            Adapter('Barcode 72 (reverse)',

            start_sequence=('repBC72_rev', ' GTCGGTATGGAAGACAGTCAGCTA'),

            end_sequence=('repBC72', ' TAGCTGACTGTCTTCCATACCGAC ')),

            Adapter('Barcode 73 (reverse)',

            start_sequence=('repBC73_rev', 'GAGGGTTCTGTCATCCTGTTTCTT'),

            end_sequence=('repBC73', 'AAGAAACAGGATGACAGAACCCTC')),

            Adapter('Barcode 74 (reverse)',

            start_sequence=('repBC74_rev', 'AGTGGAAGTGTTGGGATGCTTGTA'),

            end_sequence=('repBC74', 'TACAAGCATCCCAACACTTCCACT')),

            Adapter('Barcode 75 (reverse)',

            start_sequence=('repBC75_rev', 'ACAACAGGGTTCATCACAATGGTC'),

            end_sequence=('repBC75', 'GACCATTGTGATGAACCCTGTTGT')),

            Adapter('Barcode 76 (reverse)',

            start_sequence=('repBC76_rev', 'GTCCAGGGTTGATGTAACAAGCAT'),

            end_sequence=('repBC76', 'ATGCTTGTTACATCAACCCTGGAC')),

            Adapter('Barcode 77 (reverse)',

            start_sequence=('repBC77_rev', 'GTTGTATCCCTGAGAAACAGGTCG'),

            end_sequence=('repBC77', 'CGACCTGTTTCTCAGGGATACAAC')),

            Adapter('Barcode 78 (reverse)',

            start_sequence=('repBC78_rev', 'TTCTGATTCAAAGGTTCGGTTGTT'),

            end_sequence=('repBC78', 'AACAACCGAACCTTTGAATCAGAA')),

            Adapter('Barcode 79 (reverse)',

            start_sequence=('repBC79_rev', 'CAGCAGTGAGAACTATCTCCGAGA'),

            end_sequence=('repBC79', 'TCTCGGAGATAGTTCTCACTGCTG')),

            Adapter('Barcode 80 (reverse)',

            start_sequence=('repBC80_rev', 'GAATCGCTATCCTATGTTCATCCG'),

            end_sequence=('repBC80', 'CGGATGAACATAGGATAGCGATTC')),

            Adapter('Barcode 81 (reverse)',

            start_sequence=('repBC81_rev', 'CCGAAACAACTTCACAAGATGAGG'),

            end_sequence=('repBC81', 'CCTCATCTTGTGAAGTTGTTTCGG')),

            Adapter('Barcode 82 (reverse)',

            start_sequence=('repBC82_rev', 'TAGTCCTGGAACTCGACATACCGT'),

            end_sequence=('repBC82', 'ACGGTATGTCGAGTTCCAGGACTA')),

            Adapter('Barcode 83 (reverse)',

            start_sequence=('repBC83_rev', 'TTCGACCTTACCTAGATCAAGCCA'),

            end_sequence=('repBC83', 'TGGCTTGATCTAGGTAAGGTCGAA')),

            Adapter('Barcode 84 (reverse)',

            start_sequence=('repBC84_rev', 'TGGCACAGGTTCTAGGTCCACTAC'),

            end_sequence=('repBC84', 'GTAGTGGACCTAGAACCTGTGCCA')),

            Adapter('Barcode 85 (reverse)',

            start_sequence=('repBC85_rev', 'GATCATCCAACTAACTCCTCCGTT'),

            end_sequence=('repBC85', 'AACGGAGGAGTTAGTTGGATGATC')),

            Adapter('Barcode 86 (reverse)',

            start_sequence=('repBC86_rev', 'TACTTACGCTTGTTGGGATCACCT'),

            end_sequence=('repBC86', 'AGGTGATCCCAACAAGCGTAAGTA')),

            Adapter('Barcode 87 (reverse)',

            start_sequence=('repBC87_rev', 'CCTCCCTAACAACAGGAGCATGTA'),

            end_sequence=('repBC87', 'TACATGCTCCTGTTGTTAGGGAGG')),

            Adapter('Barcode 88 (reverse)',

            start_sequence=('repBC88_rev', 'CTGCTTCGGATCGGTAGTAGAAGA'),

            end_sequence=('repBC88', 'TCTTCTACTACCGATCCGAAGCAG')),

            Adapter('Barcode 89 (reverse)',

            start_sequence=('repBC89_rev', 'CAACTAGCCAAACATTGATGCTGT'),

            end_sequence=('repBC89', 'ACAGCATCAATGTTTGGCTAGTTG')),

            Adapter('Barcode 90 (reverse)',

            start_sequence=('repBC90_rev', 'GCCTCAAACCGTACCCTCTACATC'),

            end_sequence=('repBC90', 'GATGTAGAGGGTACGGTTTGAGGC')),

            Adapter('Barcode 91 (reverse)',

            start_sequence=('repBC91_rev', 'AGTAGCGTGAGTTCCTATGGAGCC'),

            end_sequence=('repBC91', 'GGCTCCATAGGAACTCACGCTACT')),

            Adapter('Barcode 92 (reverse)',

            start_sequence=('repBC92_rev', 'GGTCCTGTATCTTTCCACTCACAA'),

            end_sequence=('repBC92', 'TTGTGAGTGGAAAGATACAGGACC')),

            Adapter('Barcode 93 (reverse)',

            start_sequence=('repBC93_rev', 'CCCAAGTCTGAAGTGATGGAAACT'),

            end_sequence=('repBC93', 'AGTTTCCATCACTTCAGACTTGGG')),

            Adapter('Barcode 94 (reverse)',

            start_sequence=('repBC94_rev', 'GTAGGTGGCAGTTTGAGGACAATC'),

            end_sequence=('repBC94', 'GATTGTCCTCAAACTGCCACCTAC')),

            Adapter('Barcode 95 (reverse)',

            start_sequence=('repBC95_rev', 'AAGTCCATTCTTCTTCCAGACAGG'),

            end_sequence=('repBC95', 'CCTGTCTGGAAGAAGAATGGACTT')),

            Adapter('Barcode 96 (reverse)',

            start_sequence=('repBC96_rev', 'ATGGTGGACTCTATGACCGTTCAG'),

            end_sequence=('repBC96', 'CTGAACGGTCATAGAGTCCACCAT'))
            
            Adapter('Barcode 97(reverse)',

            start_sequence=('repBC97_rev', 'TTCTTTCAACAGCCACAGAAACAC'),

            end_sequence=('repBC97', 'GTGTTTCTGTGGCTGTTGAAAGAA')),

            Adapter('Barcode 98(reverse)',

            start_sequence=('repBC98_rev', 'AGCTAAGGCAAACATCAGCAGACA'),

            end_sequence=('repBC98', 'TGTCTGCTGATGTTTGCCTTAGCT')),

            Adapter('Barcode 99(reverse)',

            start_sequence=('repBC99_rev', 'CTCAGAACACAGGGTCAATGGTCC'),

            end_sequence=('repBC99', 'GGACCATTGACCCTGTGTTCTGAG')),

            Adapter('Barcode 100(reverse)',

            start_sequence=('repBC100_rev', 'AAGCCTAAGATAGCACAAAGGGAT'),

            end_sequence=('repBC100', 'ATCCCTTTGTGCTATCTTAGGCTT')),

            Adapter('Barcode 101(reverse)',

            start_sequence=('repBC101_rev', 'GAACAGGTCCCAAACACATTGGAA'),

            end_sequence=('repBC101', 'TTCCAATGTGTTTGGGACCTGTTC')),

            Adapter('Barcode 102(reverse)',

            start_sequence=('repBC102_rev', 'AAGAGCGTTTCCGTCTTTCATCAG'),

            end_sequence=('repBC102', 'CTGATGAAAGACGGAAACGCTCTT')),

            Adapter('Barcode 103(reverse)',

            start_sequence=('repBC103_rev', 'CACAATGGCACCCTTACTTAGGAA'),

            end_sequence=('repBC103', 'TTCCTAAGTAAGGGTGCCATTGTG')),

            Adapter('Barcode 104(reverse)',

            start_sequence=('repBC104_rev', 'AAGTCCCTTGTTTGGTTCAATGCA'),

            end_sequence=('repBC104', 'TGCATTGAACCAAACAAGGGACTT')),

            Adapter('Barcode 105(reverse)',

            start_sequence=('repBC105_rev', 'TTGATCCGTGTCGCTCAGAACCAA'),

            end_sequence=('repBC105', 'TTGGTTCTGAGCGACACGGATCAA')),

            Adapter('Barcode 106(reverse)',

            start_sequence=('repBC106_rev', 'TTCGCAACTTTGGAAACAGGAGAG'),

            end_sequence=('repBC106', 'CTCTCCTGTTTCCAAAGTTGCGAA')),

            Adapter('Barcode 107(reverse)',

            start_sequence=('repBC107_rev', 'CAAAGTAGATAGCCTCCCTTACCT'),

            end_sequence=('repBC107', 'AGGTAAGGGAGGCTATCTACTTTG')),

            Adapter('Barcode 108(reverse)',

            start_sequence=('repBC108_rev', 'GTCCATCTTTCTTCGTCTTAGCCT'),

            end_sequence=('repBC108', 'AGGCTAAGACGAAGAAAGATGGAC')),

            Adapter('Barcode 109(reverse)',

            start_sequence=('repBC109_rev', 'TCTTGCTGAAGGTATGAGCACACT'),

            end_sequence=('repBC109', 'AGTGTGCTCATACCTTCAGCAAGA')),

            Adapter('Barcode 110(reverse)',

            start_sequence=('repBC110_rev', 'TTGCTCAGAGAACCCTGGGTATCT'),

            end_sequence=('repBC110', 'AGATACCCAGGGTTCTCTGAGCAA')),

            Adapter('Barcode 111(reverse)',

            start_sequence=('repBC111_rev', 'TCCAGATGGAGCGATTGTGGTGAC'),

            end_sequence=('repBC111', 'GTCACCACAATCGCTCCATCTGGA')),

            Adapter('Barcode 112(reverse)',

            start_sequence=('repBC112_rev', 'GCAGTTGACTGTCACCAAGCATGA'),

            end_sequence=('repBC112', 'TCATGCTTGGTGACAGTCAACTGC')),

            Adapter('Barcode 113(reverse)',

            start_sequence=('repBC113_rev', 'TGGGAGGTCCTTTCATGGAGACTA'),

            end_sequence=('repBC113', 'TAGTCTCCATGAAAGGACCTCCCA')),

            Adapter('Barcode 114(reverse)',

            start_sequence=('repBC114_rev', 'GGTTTGGGTTGTTGGATCTATCCG'),

            end_sequence=('repBC114', 'CGGATAGATCCAACAACCCAAACC')),

            Adapter('Barcode 115(reverse)',

            start_sequence=('repBC115_rev', 'CAAGGAGCACGTCACAGTTCTCTA'),

            end_sequence=('repBC115', 'TAGAGAACTGTGACGTGCTCCTTG')),

            Adapter('Barcode 116(reverse)',

            start_sequence=('repBC116_rev', 'AACGCAGGACAATGCTCTTGAGTA'),

            end_sequence=('repBC116', 'TACTCAAGAGCATTGTCCTGCGTT')),

            Adapter('Barcode 117(reverse)',

            start_sequence=('repBC117_rev', 'CTCGGAGAGTAACAGGCAAGAGAT'),

            end_sequence=('repBC117', 'ATCTCTTGCCTGTTACTCTCCGAG')),

            Adapter('Barcode 118(reverse)',

            start_sequence=('repBC118_rev', 'TGGTGACGGTACATAGTTTCATGC'),

            end_sequence=('repBC118', 'GCATGAAACTATGTACCGTCACCA')),

            Adapter('Barcode 119(reverse)',

            start_sequence=('repBC119_rev', 'GAATGATGGGTCACTTGGAGGAGC'),

            end_sequence=('repBC119', 'GCTCCTCCAAGTGACCCATCATTC')),

            Adapter('Barcode 120(reverse)',

            start_sequence=('repBC120_rev', 'CGTATCAAGACGTACTACCCAATC'),

            end_sequence=('repBC120', 'GATTGGGTAGTACGTCTTGATACG')),

            Adapter('Barcode 121(reverse)',

            start_sequence=('repBC121_rev', 'CATTCAACCCATACGTTGCGTTAC'),

            end_sequence=('repBC121', 'GTAACGCAACGTATGGGTTGAATG')),

            Adapter('Barcode 122(reverse)',

            start_sequence=('repBC122_rev', 'GTATGTCGCTGATGCGTAAGAGTA'),

            end_sequence=('repBC122', 'TACTCTTACGCATCAGCGACATAC')),

            Adapter('Barcode 123(reverse)',

            start_sequence=('repBC123_rev', 'GCTGCCAATCTAAGTGGAGAATGT'),

            end_sequence=('repBC123', 'ACATTCTCCACTTAGATTGGCAGC')),

            Adapter('Barcode 124(reverse)',

            start_sequence=('repBC124_rev', 'ACTTTGGATTCTTCCGTGGCATAG'),

            end_sequence=('repBC124', 'CTATGCCACGGAAGAATCCAAAGT')),

            Adapter('Barcode 125(reverse)',

            start_sequence=('repBC125_rev', 'GATCTGTGGAACCCAACTGTCTGG'),

            end_sequence=('repBC125', 'CCAGACAGTTGGGTTCCACAGATC')),

            Adapter('Barcode 126(reverse)',

            start_sequence=('repBC126_rev', 'AGTCACTCCTAGATGAAGCTGGGT'),

            end_sequence=('repBC126', 'ACCCAGCTTCATCTAGGAGTGACT')),

            Adapter('Barcode 127(reverse)',

            start_sequence=('repBC127_rev', 'ACGCATGTCGTTAGTCAATGTAAC'),

            end_sequence=('repBC127', 'GTTACATTGACTAACGACATGCGT')),

            Adapter('Barcode 128(reverse)',

            start_sequence=('repBC128_rev', 'GGTCATCTTCAGGCTGTTGCAGTA'),

            end_sequence=('repBC128', 'TACTGCAACAGCCTGAAGATGACC')),

            Adapter('Barcode 129(reverse)',

            start_sequence=('repBC129_rev', 'GTCTGAACCATGCCAACCCATTGA'),

            end_sequence=('repBC129', 'TCAATGGGTTGGCATGGTTCAGAC')),

            Adapter('Barcode 130(reverse)',

            start_sequence=('repBC130_rev', 'CCTGCTTCTTGAGTTCAGTTTCCG'),

            end_sequence=('repBC130', 'CGGAAACTGAACTCAAGAAGCAGG')),

            Adapter('Barcode 131(reverse)',

            start_sequence=('repBC131_rev', 'GATGAATGCTTCGACTCCCTGACG'),

            end_sequence=('repBC131', 'CGTCAGGGAGTCGAAGCATTCATC')),

            Adapter('Barcode 132(reverse)',

            start_sequence=('repBC132_rev', 'TACAGGGTCAATCTCCTCCTTTGT'),

            end_sequence=('repBC132', 'ACAAAGGAGGAGATTGACCCTGTA')),

            Adapter('Barcode 133(reverse)',

            start_sequence=('repBC133_rev', 'CGAACGCTAACTACGAATCATAGT'),

            end_sequence=('repBC133', 'ACTATGATTCGTAGTTAGCGTTCG')),

            Adapter('Barcode 134(reverse)',

            start_sequence=('repBC134_rev', 'TGGTGTCCTCCTGCTATGTCTCTT'),

            end_sequence=('repBC134', 'AAGAGACATAGCAGGAGGACACCA')),

            Adapter('Barcode 135(reverse)',

            start_sequence=('repBC135_rev', 'GGTGTCACAGTTGATCTCGGAGAG'),

            end_sequence=('repBC135', 'CTCTCCGAGATCAACTGTGACACC')),

            Adapter('Barcode 136(reverse)',

            start_sequence=('repBC136_rev', 'ATCAAACCTACTGGTTCCTATCGG'),

            end_sequence=('repBC136', 'CCGATAGGAACCAGTAGGTTTGAT')),

            Adapter('Barcode 137(reverse)',

            start_sequence=('repBC137_rev', 'CCTCAAGCAGGTCTCTTCATGTGC'),

            end_sequence=('repBC137', 'GCACATGAAGAGACCTGCTTGAGG')),

            Adapter('Barcode 138(reverse)',

            start_sequence=('repBC138_rev', 'GATGCACATTCCGTATGGACGGTC'),

            end_sequence=('repBC138', 'GACCGTCCATACGGAATGTGCATC')),

            Adapter('Barcode 139(reverse)',

            start_sequence=('repBC139_rev', 'GAAAGCAACAACTGAGCTGCCATC'),

            end_sequence=('repBC139', 'GATGGCAGCTCAGTTGTTGCTTTC')),

            Adapter('Barcode 140(reverse)',

            start_sequence=('repBC140_rev', 'TCATCTTTCCCAAGGAAGGGTGAG'),

            end_sequence=('repBC140', 'CTCACCCTTCCTTGGGAAAGATGA')),

            Adapter('Barcode 141(reverse)',

            start_sequence=('repBC141_rev', 'CTAGGTTGTCTCTACGGAAGTCAC'),

            end_sequence=('repBC141', 'GTGACTTCCGTAGAGACAACCTAG')),

            Adapter('Barcode 142(reverse)',

            start_sequence=('repBC142_rev', 'CGACACAAGGTGAAGTAAGAGGAC'),

            end_sequence=('repBC142', 'GTCCTCTTACTTCACCTTGTGTCG')),

            Adapter('Barcode 143(reverse)',

            start_sequence=('repBC143_rev', 'CACGTTGAAAGGGTGTCCATCAAG'),

            end_sequence=('repBC143', 'CTTGATGGACACCCTTTCAACGTG')),

            Adapter('Barcode 144(reverse)',

            start_sequence=('repBC144_rev', 'GTAGACCTTGCACCATGTGGACAT'),

            end_sequence=('repBC144', 'ATGTCCACATGGTGCAAGGTCTAC')),

            Adapter('Barcode 145(reverse)',

            start_sequence=('repBC145_rev', 'TGACCACGTCGAAACTTGTAGATC'),

            end_sequence=('repBC145', 'GATCTACAAGTTTCGACGTGGTCA')),

            Adapter('Barcode 146(reverse)',

            start_sequence=('repBC146_rev', 'TACCTGAAACCATTGAAGGACGCA'),

            end_sequence=('repBC146', 'TGCGTCCTTCAATGGTTTCAGGTA')),

            Adapter('Barcode 147(reverse)',

            start_sequence=('repBC147_rev', 'CAACTTACTCGGATGACCCAGGAG'),

            end_sequence=('repBC147', 'CTCCTGGGTCATCCGAGTAAGTTG')),

            Adapter('Barcode 148(reverse)',

            start_sequence=('repBC148_rev', 'ACTCTCTGTTCTAACAAGCACCTG'),

            end_sequence=('repBC148', 'CAGGTGCTTGTTAGAACAGAGAGT')),

            Adapter('Barcode 149(reverse)',

            start_sequence=('repBC149_rev', 'TCTAAGTCTGGCAGAGTACGTTTC'),

            end_sequence=('repBC149', 'GAAACGTACTCTGCCAGACTTAGA')),

            Adapter('Barcode 150(reverse)',

            start_sequence=('repBC150_rev', 'GTTCTCGAAACTGATTCCTCGTAC'),

            end_sequence=('repBC150', 'GTACGAGGAATCAGTTTCGAGAAC')),

            Adapter('Barcode 151(reverse)',

            start_sequence=('repBC151_rev', 'ACCTTCTACTCTGGGACTAGATGC'),

            end_sequence=('repBC151', 'GCATCTAGTCCCAGAGTAGAAGGT')),

            Adapter('Barcode 152(reverse)',

            start_sequence=('repBC152_rev', 'AGTGATGAGTTGTCCACCGTACTT'),

            end_sequence=('repBC152', 'AAGTACGGTGGACAACTCATCACT')),

            Adapter('Barcode 153(reverse)',

            start_sequence=('repBC153_rev', 'CGATCCAGTTAGAGGAAGCCTTCA'),

            end_sequence=('repBC153', 'TGAAGGCTTCCTCTAACTGGATCG')),

            Adapter('Barcode 154(reverse)',

            start_sequence=('repBC154_rev', 'GTCCAATGAGGAGGCACTCAGACT'),

            end_sequence=('repBC154', 'AGTCTGAGTGCCTCCTCATTGGAC')),

            Adapter('Barcode 155(reverse)',

            start_sequence=('repBC155_rev', 'AGTTAGTTCTTCCCTTTCGTTCCA'),

            end_sequence=('repBC155', 'TGGAACGAAAGGGAAGAACTAACT')),

            Adapter('Barcode 156(reverse)',

            start_sequence=('repBC156_rev', 'GTACAAGTTGGTTCCGAAGATACC'),

            end_sequence=('repBC156', 'GGTATCTTCGGAACCAACTTGTAC')),

            Adapter('Barcode 157(reverse)',

            start_sequence=('repBC157_rev', 'TCTCCCATGATACACGGAGTCGTG'),

            end_sequence=('repBC157', 'CACGACTCCGTGTATCATGGGAGA')),

            Adapter('Barcode 158(reverse)',

            start_sequence=('repBC158_rev', 'GTGGGTGTGAATGAAGTCCTGCAT'),

            end_sequence=('repBC158', 'ATGCAGGACTTCATTCACACCCAC')),

            Adapter('Barcode 159(reverse)',

            start_sequence=('repBC159_rev', 'AAGACTTCAAGGACCCAGAACTTG'),

            end_sequence=('repBC159', 'CAAGTTCTGGGTCCTTGAAGTCTT')),

            Adapter('Barcode 160(reverse)',

            start_sequence=('repBC160_rev', 'CTGTCTGTGGCAAGTAGCTGAAAG'),

            end_sequence=('repBC160', 'CTTTCAGCTACTTGCCACAGACAG')),

            Adapter('Barcode 161(reverse)',

            start_sequence=('repBC161_rev', 'AAGAGTCAGAAGGAGGTCTGTTCC'),

            end_sequence=('repBC161', 'GGAACAGACCTCCTTCTGACTCTT')),

            Adapter('Barcode 162(reverse)',

            start_sequence=('repBC162_rev', 'GGCTAGGAACACCGAAGATTGAAG'),

            end_sequence=('repBC162', 'CTTCAATCTTCGGTGTTCCTAGCC')),

            Adapter('Barcode 163(reverse)',

            start_sequence=('repBC163_rev', 'CAAACAGTATGAGCACACGAGTGG'),

            end_sequence=('repBC163', 'CCACTCGTGTGCTCATACTGTTTG')),

            Adapter('Barcode 164(reverse)',

            start_sequence=('repBC164_rev', 'CTTAGATTCGTTTGTGCTTCCACC'),

            end_sequence=('repBC164', 'GGTGGAAGCACAAACGAATCTAAG')),

            Adapter('Barcode 165(reverse)',

            start_sequence=('repBC165_rev', 'ATGTCAGGCTCGGAGTACACTAGA'),

            end_sequence=('repBC165', 'TCTAGTGTACTCCGAGCCTGACAT')),

            Adapter('Barcode 166(reverse)',

            start_sequence=('repBC166_rev', 'TGGCTCTAGGATGCTTACCTCACA'),

            end_sequence=('repBC166', 'TGTGAGGTAAGCATCCTAGAGCCA')),

            Adapter('Barcode 167(reverse)',

            start_sequence=('repBC167_rev', 'GGACCCTCGTAGTCCATCATTGTC'),

            end_sequence=('repBC167', 'GACAATGATGGACTACGAGGGTCC')),

            Adapter('Barcode 168(reverse)',

            start_sequence=('repBC168_rev', 'ATCGACTGACAGAAGGTATGGCTG'),

            end_sequence=('repBC168', 'CAGCCATACCTTCTGTCAGTCGAT')),

            Adapter('Barcode 169(reverse)',

            start_sequence=('repBC169_rev', 'TTCTTTGTCCTACTGTCTTGGGAG'),

            end_sequence=('repBC169', 'CTCCCAAGACAGTAGGACAAAGAA')),

            Adapter('Barcode 170(reverse)',

            start_sequence=('repBC170_rev', 'ATGTTCGTAGGGTTGTGAAGGTGA'),

            end_sequence=('repBC170', 'TCACCTTCACAACCCTACGAACAT')),

            Adapter('Barcode 171(reverse)',

            start_sequence=('repBC171_rev', 'CTGGTAACACTACTTGGGACAACA'),

            end_sequence=('repBC171', 'TGTTGTCCCAAGTAGTGTTACCAG')),

            Adapter('Barcode 172(reverse)',

            start_sequence=('repBC172_rev', 'TACGAACAATGTAGTTGGGACCTG'),

            end_sequence=('repBC172', 'CAGGTCCCAACTACATTGTTCGTA')),

            Adapter('Barcode 173(reverse)',

            start_sequence=('repBC173_rev', 'GCTGGACAAAGAGTCCCTATGTTG'),

            end_sequence=('repBC173', 'CAACATAGGGACTCTTTGTCCAGC')),

            Adapter('Barcode 174(reverse)',

            start_sequence=('repBC174_rev', 'TTGTTGGCTTGGAAACTTAGTCTT'),

            end_sequence=('repBC174', 'AAGACTAAGTTTCCAAGCCAACAA')),

            Adapter('Barcode 175(reverse)',

            start_sequence=('repBC175_rev', 'AGAGCCTCTATCAAGAGTGACGAC'),

            end_sequence=('repBC175', 'GTCGTCACTCTTGATAGAGGCTCT')),

            Adapter('Barcode 176(reverse)',

            start_sequence=('repBC176_rev', 'GCCTACTTGTATCCTATCGCTAAG'),

            end_sequence=('repBC176', 'CTTAGCGATAGGATACAAGTAGGC')),

            Adapter('Barcode 177(reverse)',

            start_sequence=('repBC177_rev', 'GGAGTAGAACACTTCAACAAAGCC'),

            end_sequence=('repBC177', 'GGCTTTGTTGAAGTGTTCTACTCC')),

            Adapter('Barcode 178(reverse)',

            start_sequence=('repBC178_rev', 'TGCCATACAGCTCAAGGTCCTGAT'),

            end_sequence=('repBC178', 'ATCAGGACCTTGAGCTGTATGGCA')),

            Adapter('Barcode 179(reverse)',

            start_sequence=('repBC179_rev', 'ACCGAACTAGATCCATTCCAGCTT'),

            end_sequence=('repBC179', 'AAGCTGGAATGGATCTAGTTCGGT')),

            Adapter('Barcode 180(reverse)',

            start_sequence=('repBC180_rev', 'CATCACCTGGATCTTGGACACGGT'),

            end_sequence=('repBC180', 'ACCGTGTCCAAGATCCAGGTGATG')),

            Adapter('Barcode 181(reverse)',

            start_sequence=('repBC181_rev', 'TTGCCTCCTCAATCAACCTACTAG'),

            end_sequence=('repBC181', 'CTAGTAGGTTGATTGAGGAGGCAA')),

            Adapter('Barcode 182(reverse)',

            start_sequence=('repBC182_rev', 'TCCACTAGGGTTGTTCGCATTCAT'),

            end_sequence=('repBC182', 'ATGAATGCGAACAACCCTAGTGGA')),

            Adapter('Barcode 183(reverse)',

            start_sequence=('repBC183_rev', 'ATGTACGAGGACAACAATCCCTCC'),

            end_sequence=('repBC183', 'GGAGGGATTGTTGTCCTCGTACAT')),

            Adapter('Barcode 184(reverse)',

            start_sequence=('repBC184_rev', 'AGAAGATGATGGCTAGGCTTCGTC'),

            end_sequence=('repBC184', 'GACGAAGCCTAGCCATCATCTTCT')),

            Adapter('Barcode 185(reverse)',

            start_sequence=('repBC185_rev', 'TGTCGTAGTTACAAACCGATCAAC'),

            end_sequence=('repBC185', 'GTTGATCGGTTTGTAACTACGACA')),

            Adapter('Barcode 186(reverse)',

            start_sequence=('repBC186_rev', 'CTACATCTCCCATGCCAAACTCCG'),

            end_sequence=('repBC186', 'CGGAGTTTGGCATGGGAGATGTAG')),

            Adapter('Barcode 187(reverse)',

            start_sequence=('repBC187_rev', 'CCGAGGTATCCTTGAGTGCGATGA'),

            end_sequence=('repBC187', 'TCATCGCACTCAAGGATACCTCGG')),

            Adapter('Barcode 188(reverse)',

            start_sequence=('repBC188_rev', 'AACACTCACCTTTCTATGTCCTGG'),

            end_sequence=('repBC188', 'CCAGGACATAGAAAGGTGAGTGTT')),

            Adapter('Barcode 189(reverse)',

            start_sequence=('repBC189_rev', 'TCAAAGGTAGTGAAGTCTGAACCC'),

            end_sequence=('repBC189', 'GGGTTCAGACTTCACTACCTTTGA')),

            Adapter('Barcode 190(reverse)',

            start_sequence=('repBC190_rev', 'CTAACAGGAGTTTGACGGTGGATG'),

            end_sequence=('repBC190', 'CATCCACCGTCAAACTCCTGTTAG')),

            Adapter('Barcode 191(reverse)',

            start_sequence=('repBC191_rev', 'GGACAGACCTTCTTCTTACCTGAA'),

            end_sequence=('repBC191', 'TTCAGGTAAGAAGAAGGTCTGTCC')),

            Adapter('Barcode 192(reverse)',

            start_sequence=('repBC192_rev', 'GACTTGCCAGTATCTCAGGTGGTA'),

            end_sequence=('repBC192', 'TACCACCTGAGATACTGGCAAGTC'))]

def make_full_native_barcode_adapter(barcode_num):
    barcode = [x for x in ADAPTERS if x.name == 'Barcode ' + str(barcode_num) + ' (reverse)'][0]
    start_barcode_seq = barcode.start_sequence[1]
    end_barcode_seq = barcode.end_sequence[1]

    start_full_seq = 'AATGTACTTCGTTCAGTTACGTATTGCTAAGGTTAA' + start_barcode_seq + 'CAGCACCT'
    end_full_seq = 'AGGTGCTG' + end_barcode_seq + 'TTAACCTTAGCAATACGTAACTGAACGAAGT'

    return Adapter('Native barcoding ' + str(barcode_num) + ' (full sequence)',
                   start_sequence=('NB' + '%02d' % barcode_num + '_start', start_full_seq),
                   end_sequence=('NB' + '%02d' % barcode_num + '_end', end_full_seq))


def make_full_rapid_barcode_adapter(barcode_num):
    barcode = [x for x in ADAPTERS if x.name == 'Barcode ' + str(barcode_num) + ' (forward)'][0]
    start_barcode_seq = barcode.start_sequence[1]

    start_full_seq = 'AATGTACTTCGTTCAGTTACGTATTGCT' + start_barcode_seq + 'GTTTTCGCATTTATCGTGAAACGCTTTCGCGTTTTTCGTGCGCCGCTTCA'

    return Adapter('Rapid barcoding ' + str(barcode_num) + ' (full sequence)',
                   start_sequence=('RB' + '%02d' % barcode_num + '_full', start_full_seq))
