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
import os
import sys
from graftm.pplacer import Pplacer


sys.path = [os.path.join(os.path.dirname(os.path.realpath(__file__)),'..')]+sys.path

class HmmsearcherTests(unittest.TestCase):
    path_to_data = os.path.join(os.path.dirname(os.path.realpath(__file__)),'data')
    
    def test_basic_split(self):
        input_alias_hash = {'0':{'place':[]}, '1':{'place':[]}}
        placement=['p__Proteobacteria', 'k__Bacteria', 'p__Proteobacteria']

        test_json = {"fields":
                      ["classification", "distal_length", "edge_num", "like_weight_ratio",
                        "likelihood", "pendant_length"
                      ], "tree":
                      "((696036:0.2205{0},229854:0.20827{1})1.000:0.14379{2},3190878:0.23845{3},2107103:0.32104{4}){5};",
                      "placements":
                      [
                        {"p":
                          [
                            ["p__Proteobacteria", 0.107586583111, 1, 0.970420466541,
                              -614.032176075, 0.22226616471
                            ],
                            ["k__Bacteria", 0.220493270874, 0, 0.0147918928965, -618.21582627,
                              0.248671444337
                            ],
                            ["p__Proteobacteria", 8.77624511719e-06, 2, 0.0147876405626,
                              -618.216113788, 0.248672441848
                            ]
                          ], "nm":
                          [["test_read1_0", 1],
                            ["test_read2_1", 1]
                          ]
                        }
                      ], "version": 3, "metadata":
                      {"invocation":
                        "pplacer -c test_16S.gpkg\/test_16S.gpkg.refpkg\/ GraftM_output\/combined_alignment.aln.fa"
                      }
                    }
        
        pplacer = Pplacer("refpkg_decoy")
        
        output_alias_hash = pplacer.jplace_split(test_json, input_alias_hash)
        
        for alias_idx, alias_placement_entry in output_alias_hash.iteritems():
            
            placement_dict = alias_placement_entry['place']
            self.assertEqual(len(placement_dict), 1)
            self.assertEqual([x[0] for x in placement_dict[0]['p']],
                             placement)
        
    def test_basic_split_with_different_placements(self):
        input_alias_hash = {'0':{'place':[]}, '1':{'place':[]}}
        placements1=[['p__Crenarchaeota', 'k__Archaea', 'p__Crenarchaeota'],
                     ['p__Proteobacteria', 'k__Bacteria', 'c__Alphaproteobacteria']]
        placements2=[['p__Crenarchaeota', 'k__Archaea', 'p__Crenarchaeota'],
                     ['k__Bacteria', 'p__Cyanobacteria', 'c__Chloroplast']]
        placement_results = {'0': placements1, '1':placements2}
        test_json = {"fields":
                      ["classification", "distal_length", "edge_num", "like_weight_ratio",
                        "likelihood", "pendant_length"
                      ], "tree":
                      "((696036:0.2205{0},229854:0.20827{1})1.000:0.14379{2},3190878:0.23845{3},2107103:0.32104{4}){5};",
                      "placements":
                      [
                        {"p":
                          [
                            ["p__Crenarchaeota", 0.107586583111, 1, 0.970420466541,
                              -614.032176075, 0.22226616471
                            ],
                            ["k__Archaea", 0.220493270874, 0, 0.0147918928965, -618.21582627,
                              0.248671444337
                            ],
                            ["p__Crenarchaeota", 8.77624511719e-06, 2, 0.0147876405626,
                              -618.216113788, 0.248672441848
                            ]
                          ], "nm":
                          [["test_read1_0", 1],
                            ["test_read2_1", 1]
                          ]
                        }, 
                        {"p":
                          [
                            ["p__Proteobacteria", 0.107586583111, 1, 0.970420466541,
                              -614.032176075, 0.22226616471
                            ],
                            ["k__Bacteria", 0.220493270874, 0, 0.0147918928965, -618.21582627,
                              0.248671444337
                            ],
                            ["c__Alphaproteobacteria", 8.77624511719e-06, 2, 0.0147876405626,
                              -618.216113788, 0.248672441848
                            ]
                          ], "nm":
                          [["test_read3_0", 1]]
                        },
                        {"p":
                          [
                            ["k__Bacteria", 0.107586583111, 1, 0.970420466541,
                              -614.032176075, 0.22226616471
                            ],
                            ["p__Cyanobacteria", 0.220493270874, 0, 0.0147918928965, -618.21582627,
                              0.248671444337
                            ],
                            ["c__Chloroplast", 8.77624511719e-06, 2, 0.0147876405626,
                              -618.216113788, 0.248672441848
                            ]
                          ], "nm":
                          [["test_read4_1", 1]]
                        }
                      ], "version": 3, "metadata":
                      {"invocation":
                        "pplacer -c test_16S.gpkg\/test_16S.gpkg.refpkg\/ GraftM_output\/combined_alignment.aln.fa"
                      }
                    }
        
        pplacer = Pplacer("refpkg_decoy")
        output_alias_hash = pplacer.jplace_split(test_json, input_alias_hash)
        for alias_idx, alias_placement_entry in output_alias_hash.iteritems():
            placement_dict = alias_placement_entry['place']
            self.assertEqual(len(placement_dict), 2)
            self.assertEqual([x[0] for x in placement_dict[0]['p']],
                             placement_results[alias_idx][0])
            self.assertEqual([x[0] for x in placement_dict[1]['p']],
                             placement_results[alias_idx][1])

if __name__ == "__main__":
    unittest.main()
