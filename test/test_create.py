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
import tempdir
import sys
import extern
import logging
import tempfile
from Bio import SeqIO
from Bio.Seq import Seq
from dendropy import Tree

sys.path = [os.path.join(os.path.dirname(os.path.realpath(__file__)),'..')]+sys.path
from graftm.create import Create
from graftm.graftm_package import GraftMPackageVersion2, GraftMPackage
from graftm.sequence_io import Sequence
from graftm.external_program_suite import ExternalProgramSuite

prerequisites = ExternalProgramSuite(['taxit', 'FastTreeMP', 'hmmalign', 'mafft'])
path_to_script = os.path.join(os.path.dirname(os.path.realpath(__file__)),'..','bin','graftM')
path_to_data = os.path.join(os.path.dirname(os.path.realpath(__file__)),'data')


class Tests(unittest.TestCase):

    def test_hello_world(self):
        with tempdir.TempDir() as tmp:
            with tempdir.TempDir() as tmp2:
                cmd1 = "%s create --verbosity 2 --sequences %s --alignment %s --taxonomy %s --rerooted_tree %s --output %s" \
                    %(path_to_script,
                      os.path.join(path_to_data,'create','homologs.trimmed.unaligned.faa'),
                      os.path.join(path_to_data,'create','homologs.trimmed.aligned.faa'),
                      os.path.join(path_to_data,'create','homologs.tax2tree.rerooted.decorated.tree-consensus-strings'),
                      os.path.join(path_to_data,'create','homologstre.tree'),
                      tmp+".gpkg")
                extern.run(cmd1)
                cmd2 = "%s graft --verbosity 2 --graftm_package %s --forward %s --output_directory %s" \
                    % (path_to_script,
                       "%s.gpkg" % tmp,
                       os.path.join(path_to_data,'create','test.faa'),
                       tmp2+"_")
                extern.run(cmd2)

    def test_rerooted_tree_with_node_names(self):
        with tempdir.TempDir() as tmp:
            with tempdir.TempDir() as tmp2:
                
                cmd1 = "%s create --verbosity 2 --sequences %s --alignment %s --taxonomy %s --rerooted_tree %s --output %s" \
                    %(path_to_script,
                      os.path.join(path_to_data,'create','homologs.trimmed.unaligned.faa'),
                      os.path.join(path_to_data,'create','homologs.trimmed.aligned.faa'),
                      os.path.join(path_to_data,'create','homologs.tax2tree.rerooted.decorated.tree-consensus-strings'),
                      os.path.join(path_to_data,'create','decorated.tree'),
                      tmp+".gpkg")
                extern.run(cmd1)
                cmd2 = "%s graft --verbosity 2 --graftm_package %s --forward %s --output_directory %s" \
                    % (path_to_script,
                       "%s.gpkg" % tmp,
                       os.path.join(path_to_data,'create','test.faa'),
                       tmp2+"_")
                
                extern.run(cmd2)
        
    def test_min_aligned_percent(self):
        # test it doesn't raise with a lower check limit
        with tempdir.TempDir() as tmp:
            Create(prerequisites).main(alignment=os.path.join(path_to_data,'create','homologs.trimmed.aligned.faa'),
                          taxonomy=os.path.join(path_to_data,'create','homologs.tax2tree.rerooted.decorated.tree-consensus-strings'),
                          sequences=os.path.join(path_to_data,'create','homologs.trimmed.unaligned.faa'),
                          rerooted_tree=os.path.join(path_to_data,'create','decorated.tree'),
                          min_aligned_percent=0.5,
                          prefix=tmp+".gpkg",
                          threads=5)
            original_alignment_length = len(open(os.path.join(tmp+'.gpkg',os.path.basename(tmp)+'.gpkg.refpkg','homologs_deduplicated_aligned.fasta')).readlines())

        with tempdir.TempDir() as tmp:
            Create(prerequisites).main(alignment=os.path.join(path_to_data,'create','homologs.trimmed.aligned.faa'),
                      taxonomy=os.path.join(path_to_data,'create','homologs.tax2tree.rerooted.decorated.tree-consensus-strings'),
                      sequences=os.path.join(path_to_data,'create','homologs.trimmed.unaligned.faa'),
                      #rerooted_tree=os.path.join(path_to_data,'create','decorated.tree'), 
                      min_aligned_percent=0.9,
                      prefix=tmp+".gpkg",
                      threads=5)
            self.assertEqual(original_alignment_length-4, # 2 sequences get removed
                             len(open(os.path.join(tmp+'.gpkg',os.path.basename(tmp)+'.gpkg.refpkg','homologs_deduplicated_aligned.fasta')).readlines()))
            
    def test_hmm_input(self):
        with tempdir.TempDir() as tmp:
            gpkg = tmp+".gpkg"
            Create(prerequisites).main(sequences=os.path.join(path_to_data,'create','homologs.trimmed.unaligned.faa'),
                          taxonomy=os.path.join(path_to_data,'create','homologs.tax2tree.rerooted.decorated.tree-consensus-strings'),
                          hmm=os.path.join(path_to_data, 'create', 'first5.hmm'), # an HMM created from just the first 5 sequences
                          prefix=gpkg,
                          threads=5)
            self.assertEqual('NAME  first10\n', open(GraftMPackageVersion2.acquire(gpkg).alignment_hmm_path()).readlines()[1])
            
    def test_search_hmms_input(self):
        with tempdir.TempDir() as tmp:
            gpkg = tmp+".gpkg"
            Create(prerequisites).main(sequences=os.path.join(path_to_data,'create','homologs.trimmed.unaligned.faa'),
                          taxonomy=os.path.join(path_to_data,'create','homologs.tax2tree.rerooted.decorated.tree-consensus-strings'),
                          hmm=os.path.join(path_to_data, 'create', 'first5.hmm'), # an HMM created from just the first 5 sequences
                          search_hmm_files=[os.path.join(path_to_data, 'create', 'homologs.hmm')],
                          prefix=gpkg,
                          threads=5)
            self.assertEqual('NAME  first10\n', open(GraftMPackageVersion2.acquire(gpkg).alignment_hmm_path()).readlines()[1])
            self.assertEqual(1, len(GraftMPackageVersion2.acquire(gpkg).search_hmm_paths()))
            self.assertEqual('NAME  homologs.trimmed.aligned\n', open(GraftMPackageVersion2.acquire(gpkg).search_hmm_paths()[0]).readlines()[1])
            
    def test_dna_package(self):
        with tempdir.TempDir() as tmp:
            gpkg = tmp+".gpkg"
            Create(prerequisites).main(sequences=os.path.join(path_to_data,'create','61_otus.fasta'),
                          taxtastic_taxonomy=os.path.join(path_to_data,'61_otus.gpkg','61_otus.refpkg','61_otus_taxonomy.csv'),
                          taxtastic_seqinfo=os.path.join(path_to_data,'61_otus.gpkg','61_otus.refpkg','61_otus_seqinfo.csv'),
                          alignment=os.path.join(path_to_data,'61_otus.gpkg','61_otus.refpkg','61_otus.aln.fa'),
                          prefix=gpkg,
                          threads=5)
            pkg = GraftMPackageVersion2.acquire(gpkg)
            self.assertEqual('NAME  61_otus.aln\n', open(pkg.alignment_hmm_path()).readlines()[1])
            self.assertEqual(pkg.diamond_database_path(), None)
            
                
                
    def test_create_no_alignment(self):
        with tempdir.TempDir() as tmp:
            gpkg = tmp+".gpkg"
            Create(prerequisites).main(sequences=os.path.join(path_to_data,'create','homologs.trimmed.unaligned.faa'),
                          taxonomy=os.path.join(path_to_data,'create','homologs.tax2tree.rerooted.decorated.tree-consensus-strings'),
                          prefix=gpkg,
                          threads=5)
            self.assertTrue(os.path.exists(GraftMPackageVersion2.acquire(gpkg).alignment_hmm_path()))

    def test_create_stars(self):
        with tempdir.TempDir() as tmp:
            gpkg = tmp+".gpkg"
            Create(prerequisites).main(sequences=os.path.join(path_to_data,'create','homologs.trimmed.unaligned.with_stars.faa'),
                          taxonomy=os.path.join(path_to_data,'create','homologs.tax2tree.rerooted.decorated.tree-consensus-strings'),
                          prefix=gpkg)
            self.assertTrue(os.path.exists(GraftMPackageVersion2.acquire(gpkg).alignment_hmm_path()))

    def test_create_dereplication(self):
        reads = '''>r1
MSHFLDRLTYFSAPKESFSGGHGLTNAEDRTWEDAYRNRWAHDKIVRSTHGVNCTGSCSWKIYVKGGIVTWETQQTDYPRTRWDMPNHEPRGCARGASYSWYLYSANRVKYPMVRGRLLERWRAALKVAKSPVDAWALIQQDPEARRDYQQVRGMGGFVRSTWDEVNQLIAAANVYTIKQHGPDRVIGFSPIPAMSMVSYAAGSRYLSLIGGVCMSFYDWYCDLPPSSPQVWGEQTDVPESADWYNSTFIIAWGSNVPQTRTPDAHFFTEVRYKGAKTVAVTPDYSEVAKLSDLWLHPKQGTDAALAMAMGHVALKEFYFPGNGKPRSAYFDDYARRYTDLPLLVMLKELTLPDGKVTLVPDRYVRASDFNGKMGQKSNPEWKTVAFDTAGKAVLPNGSIGFRWGEPGREDEGKWNLESKEARSNEDVKLKLSVLEDGEQTHQVVDVAFPYFGGQHNPHFPTNNQKGDVTVAKVPAVQLRLGKAGEERYAMVATVFDLQVAQYGINRGLGSGAKDFDDNAPYTPAWQESITGVPREQVITVARQFAENADKTHGKSMVIIGAAMNHWYHADMNYRGVINLLMMCGCIGQSGGGWSHYVGQEKLRPQTGWVPLAFATDWIRPSRQQNSTSFFYAHTDQWRYEKLGMGEIVSPTLSAQDRQAYSGSMIDFNVKAERMGWLPSAPQLKTNPLQVVKAAEAAGMDAKDYVVKALKDGTLEMSCEDPDASQNWPRNMFVWRSNLLGSSGKGHEYFCKHLLGTENGVQGKDLGKDDARPTEVKWHKDAPQGKLDLLVTLDFRMSTTCLYSDIVLPTASWYEKNDLNTSDMHPFIHPLSAAVDPVWQARSDWEIYKGFAKAFSDMCVGHLGVEKEVVLTPIMHDTAGEMAQPIDVREWKKGQCDLIPGKTAPQITVVERDYPNTFKRFTALGPLLDKLGNGGKGIGWQTGIEVEQLGQLNGLTQEEGISKGRPRIVSDIDACETILQLAPETNGHVAVKAWEALGKQTGRDHTHLALYREDEKIRYRDIQAQPRKIISSPTWSGIESETVSYNAGYTNVHELIPWRTLTGRQQFYQDHKWMIAFGEGFSTYRPPVNLKTTDALQNKRPNGHKEIVLNFITPHQKWGIHSTYSDNLLMLTLNRGGTVIWLSENDAKVAGIEDNDWIELFNVNGAVAGRAVVSQRVNDGMCLMYHSQEKIINTPGSEITGVRGGIHNSVTRIVTKPTHMIGGYAQLSYGFNYYGTIGTNRDEFVIIRKMVKVDWLDTPADDHLVRAYQSQGENP
>r2
MSHFLDRLRYFTAPRPQFSDGHGAVTDEDRQWENAYRQRWQHDKIVRSTHGVNCTGSCSWKIYVKGGIVTWETQQTDYPRTRPDMPNHEPRGCSRGASYSWYLYSANRLKYPLVRSALVKLWRERRVSLEPVQAWQSIVTDEAARRGYQSRRGLGGFIRSTWDEVNEIVAAANVYTVKQHGPDRVVGFSPIPAMSMVSYAAGSRYLSLIGGVCLSFYDWYCDLPPASPQTWGEQTDVPESADWYNSTFIVMWGSNVPQTRTPDAHFMTEVRYKGAKIVSIFPDYAEGAKFGDIWLHPKQGTDAALALAMGHVILKEFHLTGKSDYFVDYCRRFTDMPCLVRLVPHGDAYVPERLVRASDFDDALGQANNAEWKTVMIDAASNEFVVPLGSVGFRWGQKEGEDQGKWNLKSEAATGEPLMPRLSLANAHDDVVAVLFPYFGNIAHPHFHHTTHGSQLSRKIGVRRIATRNGEMLVATVHDLFVANYGLDQGLGGEHIARDYDDDLPYTPAWQEAITGVKRADVISVARQFAENAHKTQGKSMVIIGAGINHWFHMDMSYRAIINMLIMCGCVGKSGGGWSHYVGQEKLRPQTGWTALAFALDWHRPPRHMNATSFFYAHTDQWRYDPMNPTALLSPLADASRFHGAPIDYNVRAERMGWLPSAPQLAVNPLDYGRALNDPAQAAATVVRDLKSGSLRMACEDPDHPDNFPRNLFVWRSNLLGSSGKGHEYFLKHLLGTINGVQGEEVVKSGHPRPEEIAWHDEAPRGKLDLLVTLDFRMSTTCMYSDVVLPTATWYEKDDMNTSDMHPFIHPLSAAVDPAWQSKSDWEIFKGIAKRFSELSEGHLGVERDIVLAPIAHDSPAELAQPFDVQDWKRGECEPVPGKTMPSVLVVERDYPATYAKFTSLGPLMDKLGNGGKGVSWDTKDEVALLGDLNYRVADAGCAKGRPRIDTALDAAEVILSLAPETNGAVAVKAWSAVTAMTGIEHAHLADARADEKIRFRDIPGPAAKIISSPTWSGIESEHVSYNAGYTNVHELIPWRTLSGRQQLYQDHLWMRAFGESLCVYKPPIETGSYEHMAGARSNGNPEIVLNFITPHQKWGIHSTYTDNLLMLTLSRGGPIVWLSEADAQAIGLVDNDWIECYNANGALCARAVVSQRVPRGMVMMYHAQEKIVNTPGSEITGTRGGIHNSVTRVALKPTHMIGGYAQLAYGFNYYGTVGSNRDEFLIVRKMKHIDWLDDTPPLTMSPPVQPLRQGEAS
>r3
MSHFLDRLKFMSKVKSTFSNGHGAVVKEDRQWEDAYRQRWQHDKVVRSTHGVNCTGSCSWKVYVKNGLITWETQQTDYPRTRADLPNHEPRGCPRGASYSWYVYSAQRVKYPMMRGKLAQMWREARKTMSAVEAWEYISQDPERAREYKSRRGQGGFVRATWEEANEMVAAANAYTIKNYGPDRVIGFSPIPAMSMVSYAAGSRYLSLIGGVPLSFYDWYCDLPPASPQIWGEQTDVAESADWYNSTYLMVWGSNVPMTRTPDAHFYTEVRYKGTKTVAVSSDFGEMAKFGDIWLAPKQGTDAALAMAMGHVIFKEFHLDNPSEYFTDYIRRLTDMPMLVRLKEEGGHYLPEYFLRASQLADNLGEDNNPDWKTLMIDELTGDVVAPNGSIGFRWGQAEGKTGKWNLETREGTTDREVKGQLSLLGRHDEVVGVAFPYFGAEHDDQLTRNIPAKRVQLADGSSALVTTVFDLMAGNYGIDRGLGGGNVATSYMDDVPYTPAWQQKHTGVKPEMVIQVAREFAQNADQTQGKSMVIVGAGLNHWYHMDMSYRGIINMLMLCGCIGQSGGGWCHYVGQEKLRPQSGWAPLAFAQDWNRPARQMNGTSFFYAHTSQWRHEKLGVNEILCPTAAGNMEHMSLIDYNAKAERMGWLPSAPQLETNPLDLTRMAAEAGMDPAAYTVQGLKDGSIDMSCNDPDNPKNFPRNLFVWRSNILGSSGKGHEYFLKYLLGTQNAVMSDETQCIQPSEITVRPAAEGKLDLLVVLDFRMSTTCLYGDIVLPTATWYEKDDLNTSDMHPFIHPLSEAVQPLWQSKTDWEIYKGFAKAISEVGKDYLGVQKDLVLTPLMHDSPQELGQPFDPRDWKKGECDPIPGKTMPAMTVVERDYGAIYQKFTSIGPLMEKVGNGGKGMAWKTEREVEQLRKMNKVVADGVAKGQPRLDTAIDAAEMIMTLAPETNGHVAVKAWESLSKITGRDHTHLAIPREHDHITFRDIQAQPRKIISSPIWSGLESEEVSYNAGYTNVHELIPWRTLTGRQQFYQDHQWMRAFGEGMCVYRPAVDLKTTAAVHNHKPNGNKEILLNFLTPHQKWGIHSTYSENLRMLTLSRGGPHIWISEKDAREAGIVDNDWVEVFNVNGTLTARVVVSQRIPQGMTLMYHAQEKIVNVPGAEMSGKRGGIHNSVTRAVTKPTHMIGGYAQLAYGFNYYGTVGCNRDEFIVMRKMKNVDWMDQPLSK
>r6
MSRNDASHTDETAGESDDKRRAVDQSDDEQPKKRGANEERTDGGPAVGDPPGGDNPSRRGFLKGVGLASVLGIGSASASDDALFSMDGLQPVGDPIGEYPYRDWEDLYREQWDWDSVSRSTHSVNCTGSCSWNVYVKNGQVWREEQSGDYPRFDESLPDPNPRGCQKGACYTDYVNAEQRIKHPLKRVGERGEGKWKRITWDEALTEIAEHVVDEVEAGRYDAISGFTPIPAMSPVSFASGSRLINLLGGVSHSFYDWYSDLPPGQPITWGTQTDNAESADWYNADYIIAWGSNINVTRIPDAKYFLESGYNGTKRVGIFTDYSQTAIHTDEWLSPDPGSDTALALGMAQTIVSEGLYDEAHLKEQTDMPLLVRQDTGKFLRASDVPSVNSSADRPEWMLLMLDSNGQLREAPGSLGERDGQKDYSKSIDLDFDPRLDAETTVQTDDGSVQVRSVWAELRDELAQYDPETVTKMTGVGTETYQRVAREFADVERAKIIHGKGVNDWYHNDLGNRAIQLLVTLTGNLGRQGTGVDHYVGQEKIWTFHGWKTLSFPTGKVRGVPTTLWTYYHAGILDNTDADTAAKIRESIDKGWMPVYPEEREDGSRPDPTTMFVWRGNYFNQAKGNVAVEEELWPKLDLVVDINFRMDSTALYSDIVLPTASHYEKHDLSMTDMHTYVHPFTPAVEPLGESKTDWQIFKDLAEKIQEVATERGVEPISDRTFDREIDLQSVYDDYVRDWETETDGALAEDRAAAEYILEHSAESNPAGTDEQITFADTVEQPQRLLEAGDHWTSDIEDGEAYAPWKDFVQDKNPWPTVTGRQQYYIDHDWFLELGEQLPTHKEGPTNTGGDYPMEYNTPHGRWAIHSTWRDSEKMLRLQRGEPIVYINPDDAAERGIEDGDTVEVFNDLGAVEVQAKIYPSSERGTLRHFFSWEKFQYASRNNFNTLVPMYMKPTQLVQYPEDTGEHLYFFPNYWGPTGVNSDVRVDVRKKGGDGG
>r5
MMNNQELNPKTTRRSFIKASGMTALTVSAGINLPLQVKAETIPIKNMQDDKEVLSLCSVNCGSRCVLRLHVKNDEVRWVETDNSGDDEYGNHQIRACLRGRSMRYRMNHPDRLKYPMKRIGKRGEGKFQRISWDEALDTIADNLKRIVKDYGNEAVYNNYASGIVGGNMTRSSPFASLFTRLMNCYGGLLSYYGSYSTAQIARAMPFTYGSNVGNSTSDIVNSKLVVFFGNNPMETRMSGAGITYHLEQARERSNAKLIVIDPRYTDTASGREDEWIPIRPGTDAALVAGMAHVMITENLVDQVFLDKYCVGYDEKTLPASAPANGHYKAYILGQGDDGVEKTPEWAAKITGIPAAKIIRLAREIGTTKPCYIAQGWGLQRQSNGELACRAVAMLAILTGNVGISGGNSGAREGSFYPPIQRLPVLENPVKAKISVFSWTDAIDHGDQMTALKDGVQGKDKLDVPIKFMWNYAGNCITNQHSDINRTHEILSDDKKCEMIVVIENFMTDSAKYADILLPDLMTTEQEDIVPNDYAGNMGYLIFSQPATSAKFERRGIYDICCEIAKRLGQDVYDKFTEGGRTQEQWLQFLYSRMQEKDATLPSYDELKQMGVYKRQDPNGHLVAYQKFREDPEANPLKTPSGKIEIYSERLAEIAATWELKEDEVIHPLPIYHSGFNGWDDPKREKYPFQLVSFHYKSRTHSTYGNIDVLESACPQEMWINPIDAQRLTLTDGETVRVFNEIGETRIPVKITPRIMPGVLGMGQGAWHKANMFGDKVDHNGCVNVLTTLRPSPLAKGNPQHSNLVRIEKL
>r7
MSKIDLYLDNPAFHIVTGACPHDCPDTCTWQVAVDRITGAAVDIWGHPAHPITQGKLCGKVDRYLERTYHKDRLTVPLRRVGPKGSGRFERVSWEEAIADIARRLQAIIEEYGPEAVLPYSYSGTLGFLQGEGMASRFFNRLGASQLARTICAEAGFQGYLYTIGAAEGMLPEDFAHSRLILIWGSNTLTSNLHLWPFIQQARKQGARVLVIDPANTRTAQAGDEWIPIRPGTDAALALAMMQVIIAEERYDADYVARYTLGFKQLAERVRDWTPARSAEITGIPAERIAALAREYATIRPAGIRINYGLQRHYGGGMAVRTIACLPALVGAWRQHGGGIQLSSSGLFRHLDKRGLHRPDLLEKAKAELQKGNAGAVDQTRPFPRVINMNRLGDALSLDPARLARAHYRPRPIDPLPKPEEAGPPVKALIVYNCNPAAVAPDQTAVLAGLRREELFTVVLEHFQTDTADYADYLLPATTQLEHWDIHRAYGHAYLALNRPAIAPVGESLPNSEIFRRLALAMGYDDPCFYESDEEILRAFVAAQTHERFAGLTWERLLNEAFFRLNLPDPYLPFAEGNFPTPSGKCEFFSARMAEDGYDPLPTYTPPRWQDRSGASALSSFDGSLVCISPPAHSFLNSTFANLPRFLQREDRPMLWIHPQDARPRHIEDGTRVRVCNQHGEVMLTARVTDAVVPGTVLAPGVWWPKLSPDGRNINRTTPPDETDMGAGALFYDATVWVEPLPVMANG
>r4
MSESVVVRAACPHDCPDTCAMEIGVRDGRAVSVRGAADMPFTNGALCTKVAHYLERVYSDQRLMHPLRRVGAKGDGRFERIRWDDALDIIAERFRAIAADDPRAILPYSYAGTMGLVQSESLDRRFFNRLGAARLDRTICAAAGAMGWKATVGASIGADPEAVVDARLILIWGGNPVVSNVHGWRYLQEAKRRGAKLICIDPRRSETAAKCHQHVAPLPGTDGALALAMMQVLIAEDLLDHDYIARCTLGFAELAERVRGCTPEWAAGVTGLDAATIRELAQAYGRSAPTLIRLNYGLNRSGGGGMAVRNIACLPALTGAWRHAGGGALLSTSGNFPLDRHALERPDLYRQDPSRPPRTINMTTIGEALLETRDPPIRALYVYNANPVAVAPYSAKVRAGFRREDLFCVVHELFQTDTADYADILLPATSQLEHHDLHTSYGHFYLVANRPAIDAQGEARSNTEVFRQLAARLGFDDTALFETDEAIAAAALRREDGRSLGVGARLAEHGWARLNLPRPFAPFARGGFPTPSGKCEFLSSFLAAQGLDPLPAWHPPYESVVSNPALAARYPLALISPPARNFLNSSFANLPRFTAGEGGPRLDIHRDDAAARGLTDGMRVRIHNDRGAFHAYARVTEAVRRGVVAAPSIWWQKRSGDGENANAVTSDALTDMGGGAVFYDCLVEVAGAPA\n'''
        taxonomy = '''r1\td__Bacteria;p__Proteobacteria;c__Betaproteobacteria;o__Burkholderiales;f__Comamonadaceae;g__Hydrogenophaga;s__
r3\td__Bacteria;p__Proteobacteria;c__Betaproteobacteria;o__Burkholderiales;f__Comamonadaceae;g__Brachymonas;s__
r2\td__Bacteria;p__Proteobacteria;c__Betaproteobacteria;o__Burkholderiales;f__Burkholderiaceae;g__Burkholderia;s__
r4\td__Bacteria;p__Proteobacteria;c__Betaproteobacteria;o__Rhodocyclales;f__Rhodocyclaceae;g__;s__
r5\td__Bacteria;p__Proteobacteria;c__Gammaproteobacteria;o__Pasteurellales;f__Pasteurellaceae;g__Basfia;s__
r7\td__Bacteria;p__Chloroflexi;c__Caldilineae;o__Caldilineales;f__Caldilineaceae;g__Caldilinea;s__aerophila
r6\td__Archaea;p__Euryarchaeota;c__Methanomicrobia;o__Halobacteriales;f__Halobacteriaceae;g__Halogeometricum;s__\n'''
    
        with tempfile.NamedTemporaryFile(suffix='.fa') as fasta:
            fasta.write(reads)
            fasta.flush()
            with tempfile.NamedTemporaryFile(suffix='.tsv') as tax:
                tax.write(taxonomy)
                tax.flush()
                sequences_file = fasta.name
                taxonomy_file  = tax.name
                expected_list = [7, 2 ,3, 4, 5, 6, 7]
                for expected, i in zip( expected_list, range(0,6) ):

                    with tempdir.TempDir() as package:
                        Create(prerequisites).main(sequences = sequences_file,
                                      taxonomy = taxonomy_file,
                                      prefix = package,
                                      dereplication_level = i,
                                      force = True,
                                      threads = 5)
                        base = os.path.basename(package)
                        gpkg = GraftMPackageVersion2.acquire(package)
                        seqinfo = open(gpkg.taxtastic_seqinfo_path())\
                                                                .readlines()

                        hmm = gpkg.search_hmm_paths()[0]
                        nseq = [int(x.strip().split()[1]) 
                                for x in open(hmm).readlines()
                                if x.startswith("NSEQ")][0]

                        self.assertEqual(nseq, expected)
                        
    def test_strange_character_replace(self):
        create = Create(prerequisites)
        seqs = [Sequence('namer','SEQWENCE')]
        create._mask_strange_sequence_letters(seqs, Create._PROTEIN_PACKAGE_TYPE)
        self.assertEqual(1, len(seqs))
        self.assertEqual('SEQWENCE', str(seqs[0].seq))
        
        seqs = [Sequence('namer','SEQUENCE')]
        create._mask_strange_sequence_letters(seqs, Create._PROTEIN_PACKAGE_TYPE)
        self.assertEqual(1, len(seqs))
        self.assertEqual('SEQXENCE', str(seqs[0].seq))
        
    def test_dna_strange_character_replace(self):
        create = Create(prerequisites)
        seqs = [Sequence('namer','ATGC')]
        create._mask_strange_sequence_letters(seqs, Create._NUCLEOTIDE_PACKAGE_TYPE)
        self.assertEqual(1, len(seqs))
        self.assertEqual('ATGC', str(seqs[0].seq))
        
        seqs = [Sequence('namer','ATGCRTWU')]
        create._mask_strange_sequence_letters(seqs, Create._NUCLEOTIDE_PACKAGE_TYPE)
        self.assertEqual(1, len(seqs))
        self.assertEqual('ATGCNTNT', str(seqs[0].seq))
        
    def test_remove_strange_characters_integration_test(self):
        with tempdir.TempDir() as tmp:
            gpkg = tmp+".gpkg"
            first_seq = None
            with tempfile.NamedTemporaryFile(suffix='61_otus.Rs.fasta') as f:
                for record in SeqIO.parse(open(os.path.join(path_to_data,'create','61_otus.fasta')), 'fasta'):
                    if not first_seq:
                        first_seq = str(record.seq)
                    record.seq = Seq(str(record.seq).replace('A','R',5)) #don't replace too many otherwise hmmbuild fails
                    SeqIO.write(record, f, 'fasta')
                f.flush()

                Create(prerequisites).main(sequences=f.name,
                              taxtastic_taxonomy=os.path.join(path_to_data,'61_otus.gpkg','61_otus.refpkg','61_otus_taxonomy.csv'),
                              taxtastic_seqinfo=os.path.join(path_to_data,'61_otus.gpkg','61_otus.refpkg','61_otus_seqinfo.csv'),
                              alignment=os.path.join(path_to_data,'61_otus.gpkg','61_otus.refpkg','61_otus.aln.fa'),
                              prefix=gpkg,
                              threads=5)
                pkg = GraftMPackage.acquire(gpkg)
                self.assertEqual('NAME  61_otus.aln\n', open(pkg.alignment_hmm_path()).readlines()[1])
                self.assertEqual(pkg.diamond_database_path(), None)
                for record in SeqIO.parse(open(pkg.alignment_fasta_path()), 'fasta'):
                    self.assertEqual(str(record.seq).replace('R','N'), str(record.seq))
                    break

    def test_for_duplicate_sequence_names(self):
        with tempfile.NamedTemporaryFile(suffix='.fa') as bad_alignment:
            with tempdir.TempDir() as tmp:
                with self.assertRaises(Exception) as context:
                    Create(prerequisites).main(
                        alignment=os.path.join(path_to_data,'create','homologs.trimmed.aligned.duplicate.faa'),
                        taxonomy=os.path.join(path_to_data,'create','homologs.tax2tree.rerooted.decorated.tree-consensus-strings'),
                        sequences=os.path.join(path_to_data,'create','homologs.trimmed.unaligned.faa'),
                        rerooted_tree=os.path.join(path_to_data,'create','decorated.tree'),
                        min_aligned_percent=0.5,
                        prefix=tmp+".gpkg",
                        threads=5)

    def test_for_duplicate_alignment_names(self):
        with tempfile.NamedTemporaryFile(suffix='.fa') as bad_alignment:
            with tempdir.TempDir() as tmp:
                with self.assertRaises(Exception) as context:
                    Create(prerequisites).main(
                        alignment=os.path.join(path_to_data,'create','homologs.trimmed.aligned.faa'),
                        taxonomy=os.path.join(path_to_data,'create','homologs.tax2tree.rerooted.decorated.tree-consensus-strings'),
                        sequences=os.path.join(path_to_data,'create','homologs.trimmed.unaligned.duplicate.faa'),
                        rerooted_tree=os.path.join(path_to_data,'create','decorated.tree'),
                        min_aligned_percent=0.5,
                        prefix=tmp+".gpkg",
                        threads=5)

    def test_input_unrooted_tree(self):
        otu61 = os.path.join(path_to_data, '61_otus.gpkg','61_otus.refpkg')
        with tempfile.NamedTemporaryFile(suffix='.fa') as bad_alignment:
            with tempdir.TempDir() as tmp:
                Create(prerequisites).main(
                    taxtastic_taxonomy=os.path.join(otu61,'61_otus_taxonomy.csv'),
                    taxtastic_seqinfo=os.path.join(otu61,'61_otus_seqinfo.csv'),
                    # created with newick_utils:
                    # nw_prune test/data/61_otus.gpkg/61_otus.refpkg/61_otus.tre 4459468 >test/data/61_otus.without_4459468.tre
                    unrooted_tree=os.path.join(path_to_data,'create','61_otus.without_4459468.tre'),
                    sequences=os.path.join(path_to_data,'create','61_otus.without_4459468.fasta'),
                    alignment=os.path.join(path_to_data,'create','61_otus.without_4459468.aln.fasta'),
                    prefix=tmp, force=True)
                gpkg = GraftMPackage.acquire(tmp)
                tree=Tree.get(schema='newick', data=open(gpkg.reference_package_tree_path()).readline())
                self.assertEqual(21, len(tree.leaf_nodes()))


if __name__ == "__main__":
    logging.basicConfig(level=logging.ERROR)
    unittest.main()
