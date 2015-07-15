#!/usr/bin/env python

#=======================================================================
# Authors: Ben Woodcroft, Joel Boyd
#
# Unit tests.
#
# Copyright
#
# This is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# Distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.	See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License.
# If not, see <http://www.gnu.org/licenses/>.
#=======================================================================

import unittest
import os.path
import sys
import tempfile
path_to_script = os.path.join(os.path.dirname(os.path.realpath(__file__)),'..','bin','graftM')
path_to_data = os.path.join(os.path.dirname(os.path.realpath(__file__)),'data')

sys.path = [os.path.join(os.path.dirname(os.path.realpath(__file__)),'..')]+sys.path

from graftm.cluster import readsClusterer

class Tests(unittest.TestCase):

    def test_reads_cluster(self):
        truth_clusters={'cluster1_1': ['cluster1_2'],
                         'shouldnt_cluster1': [],
                         'shouldnt_cluster10': [],
                         'shouldnt_cluster11': [],
                         'shouldnt_cluster12': [],
                         'shouldnt_cluster13': [],
                         'shouldnt_cluster14': [],
                         'shouldnt_cluster15': [],
                         'shouldnt_cluster16': [],
                         'shouldnt_cluster2': [],
                         'shouldnt_cluster3': [],
                         'shouldnt_cluster4': [],
                         'shouldnt_cluster5': [],
                         'shouldnt_cluster6': [],
                         'shouldnt_cluster7': [],
                         'shouldnt_cluster8': [],
                         'shouldnt_cluster9': [],
                         'shouldnt_cluster_end': [],
                         'shouldnt_cluster_start': []}

        reads='''>cluster1_1
CGGCTTTCATGATGTA
>cluster1_2
CGGCTTTCATGATGTA
>shouldnt_cluster_start
CGGCTTTC
>shouldnt_cluster_end
ATGATGTA
>shouldnt_cluster1
NGGCTTTCATGATGTA
>shouldnt_cluster2
CNGCTTTCATGATGTA
>shouldnt_cluster3
CGNCTTTCATGATGTA
>shouldnt_cluster4
CGGNTTTCATGATGTA
>shouldnt_cluster5
CGGCNTTCATGATGTA
>shouldnt_cluster6
CGGCTNTCATGATGTA
>shouldnt_cluster7
CGGCTTNCATGATGTA
>shouldnt_cluster8
CGGCTTTNATGATGTA
>shouldnt_cluster9
CGGCTTTCNTGATGTA
>shouldnt_cluster10
CGGCTTTCANGATGTA
>shouldnt_cluster11
CGGCTTTCATNATGTA
>shouldnt_cluster12
CGGCTTTCATGNTGTA
>shouldnt_cluster13
CGGCTTTCATGANGTA
>shouldnt_cluster14
CGGCTTTCATGATNTA
>shouldnt_cluster15
CGGCTTTCATGATGNA
>shouldnt_cluster16
CGGCTTTCATGATGTN'''
        with tempfile.NamedTemporaryFile(suffix='.fa') as fasta:
            fasta.write(reads)
            fasta.flush()
            data = fasta.name
            clust=readsClusterer()

            clusters, uc = clust.cluster(data)
            
            result=[x.strip() for x in open(clusters).readlines()]
            
            for i in range(1,17):
                i='cluster'+str(i)
                self.assertTrue(any([x for x in result if x.endswith(i)]))
                
            for i in ['end', 'start']:
                self.assertTrue(any([x for x in result if x.endswith(i)]))
                
            self.assertFalse(any([x for x in result if x == '>cluster1_2']))
            self.assertTrue(uc==truth_clusters)              

if __name__ == "__main__":
    unittest.main()



