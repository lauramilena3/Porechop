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

            start_sequence=('repBC36_rev', 'ATGTCCCAGTTAGAGGAGGAAACA'),

            end_sequence=('repBC36', 'TGTTTCCTCCTCTAACTGGGACAT')),

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

            start_sequence=('repBC74', 'AGTGGAAGTGTTGGGATGCTTGTA'),

            end_sequence=('repBC74_rev', 'TACAAGCATCCCAACACTTCCACT')),

            Adapter('Barcode 75 (reverse)',

            start_sequence=('repBC75', 'ACAACAGGGTTCATCACAATGGTC'),

            end_sequence=('repBC75_rev', 'GACCATTGTGATGAACCCTGTTGT')),

            Adapter('Barcode 76 (reverse)',

            start_sequence=('repBC76', 'GTCCAGGGTTGATGTAACAAGCAT'),

            end_sequence=('repBC76_rev', 'ATGCTTGTTACATCAACCCTGGAC')),

            Adapter('Barcode 77 (reverse)',

            start_sequence=('repBC77', 'GTTGTATCCCTGAGAAACAGGTCG'),

            end_sequence=('repBC77_rev', 'CGACCTGTTTCTCAGGGATACAAC')),

            Adapter('Barcode 78 (reverse)',

            start_sequence=('repBC78', 'TTCTGATTCAAAGGTTCGGTTGTT'),

            end_sequence=('repBC78_rev', 'AACAACCGAACCTTTGAATCAGAA')),

            Adapter('Barcode 79 (reverse)',

            start_sequence=('repBC79', 'CAGCAGTGAGAACTATCTCCGAGA'),

            end_sequence=('repBC79_rev', 'TCTCGGAGATAGTTCTCACTGCTG')),

            Adapter('Barcode 80 (reverse)',

            start_sequence=('repBC80', 'GAATCGCTATCCTATGTTCATCCG'),

            end_sequence=('repBC80_rev', 'CGGATGAACATAGGATAGCGATTC')),

            Adapter('Barcode 81 (reverse)',

            start_sequence=('repBC81', 'CCGAAACAACTTCACAAGATGAGG'),

            end_sequence=('repBC81_rev', 'CCTCATCTTGTGAAGTTGTTTCGG')),

            Adapter('Barcode 82 (reverse)',

            start_sequence=('repBC82', 'TAGTCCTGGAACTCGACATACCGT'),

            end_sequence=('repBC82_rev', 'ACGGTATGTCGAGTTCCAGGACTA')),

            Adapter('Barcode 83 (reverse)',

            start_sequence=('repBC83', 'TTCGACCTTACCTAGATCAAGCCA'),

            end_sequence=('repBC83_rev', 'TGGCTTGATCTAGGTAAGGTCGAA')),

            Adapter('Barcode 84 (reverse)',

            start_sequence=('repBC84', 'TGGCACAGGTTCTAGGTCCACTAC'),

            end_sequence=('repBC84_rev', 'GTAGTGGACCTAGAACCTGTGCCA')),

            Adapter('Barcode 85 (reverse)',

            start_sequence=('repBC85', 'GATCATCCAACTAACTCCTCCGTT'),

            end_sequence=('repBC85_rev', 'AACGGAGGAGTTAGTTGGATGATC')),

            Adapter('Barcode 86 (reverse)',

            start_sequence=('repBC86', 'TACTTACGCTTGTTGGGATCACCT'),

            end_sequence=('repBC86_rev', 'AGGTGATCCCAACAAGCGTAAGTA')),

            Adapter('Barcode 87 (reverse)',

            start_sequence=('repBC87', 'CCTCCCTAACAACAGGAGCATGTA'),

            end_sequence=('repBC87_rev', 'TACATGCTCCTGTTGTTAGGGAGG')),

            Adapter('Barcode 88 (reverse)',

            start_sequence=('repBC88', 'CTGCTTCGGATCGGTAGTAGAAGA'),

            end_sequence=('repBC88_rev', 'TCTTCTACTACCGATCCGAAGCAG')),

            Adapter('Barcode 89 (reverse)',

            start_sequence=('repBC89', 'CAACTAGCCAAACATTGATGCTGT'),

            end_sequence=('repBC89_rev', 'ACAGCATCAATGTTTGGCTAGTTG')),

            Adapter('Barcode 90 (reverse)',

            start_sequence=('repBC90', 'GCCTCAAACCGTACCCTCTACATC'),

            end_sequence=('repBC90_rev', 'GATGTAGAGGGTACGGTTTGAGGC')),

            Adapter('Barcode 91 (reverse)',

            start_sequence=('repBC91', 'AGTAGCGTGAGTTCCTATGGAGCC'),

            end_sequence=('repBC91_rev', 'GGCTCCATAGGAACTCACGCTACT')),

            Adapter('Barcode 92 (reverse)',

            start_sequence=('repBC92', 'GGTCCTGTATCTTTCCACTCACAA'),

            end_sequence=('repBC92_rev', 'TTGTGAGTGGAAAGATACAGGACC')),

            Adapter('Barcode 93 (reverse)',

            start_sequence=('repBC93', 'CCCAAGTCTGAAGTGATGGAAACT'),

            end_sequence=('repBC93_rev', 'AGTTTCCATCACTTCAGACTTGGG')),

            Adapter('Barcode 94 (reverse)',

            start_sequence=('repBC94', 'GTAGGTGGCAGTTTGAGGACAATC'),

            end_sequence=('repBC94_rev', 'GATTGTCCTCAAACTGCCACCTAC')),

            Adapter('Barcode 95 (reverse)',

            start_sequence=('repBC95', 'AAGTCCATTCTTCTTCCAGACAGG'),

            end_sequence=('repBC95_rev', 'CCTGTCTGGAAGAAGAATGGACTT')),

            Adapter('Barcode 96 (reverse)',

            start_sequence=('repBC96', 'ATGGTGGACTCTATGACCGTTCAG'),

            end_sequence=('repBC96_rev', 'CTGAACGGTCATAGAGTCCACCAT'))]

def make_full_native_barcode_adapter(barcode_num):
    barcode = [x for x in ADAPTERS if x.name == 'Barcode ' + str(barcode_num) + ' (reverse)'][0]
    start_barcode_seq = barcode.start_sequence[1]
    end_barcode_seq = barcode.end_sequence[1]

    start_full_seq = 'AATGTACTTCGTTCAGTTACGTATTGCTAAAGGTTAA' + start_barcode_seq + 'CAGCACC'
    end_full_seq = 'GGTGCTG' + end_barcode_seq + 'TTAACCTTTGCAATACGTAACTGAACGAAGT'

    return Adapter('Native barcoding ' + str(barcode_num) + ' (full sequence)',
                   start_sequence=('NB' + '%02d' % barcode_num + '_start', start_full_seq),
                   end_sequence=('NB' + '%02d' % barcode_num + '_end', end_full_seq))


def make_full_rapid_barcode_adapter(barcode_num):
    barcode = [x for x in ADAPTERS if x.name == 'Barcode ' + str(barcode_num) + ' (forward)'][0]
    start_barcode_seq = barcode.start_sequence[1]

    start_full_seq = 'AATGTACTTCGTTCAGTTACGTATTGCT' + start_barcode_seq + 'GTTTTCGCATTTATCGTGAAACGCTTTCGCGTTTTTCGTGCGCCGCTTCA'

    return Adapter('Rapid barcoding ' + str(barcode_num) + ' (full sequence)',
                   start_sequence=('RB' + '%02d' % barcode_num + '_full', start_full_seq))
