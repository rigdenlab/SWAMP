import unittest
from swamp.parsers.phenixparser import PhenixParser


class MyTestCase(unittest.TestCase):

    def test_1(self):
        stdout_contents = """#                       get_cc_mtz_pdb
#
# Get correlation between atoms in a PDB file and map 
# offsetting the PDB file by allowed origin shifts

# Type phenix.doc for help
Values of all params:
get_cc_mtz_pdb {
  pdb_in = "/data1/filo/results/SWAMP_benchmarking/SWAMP_0/swamp_mr/search_461/run_1/refmac/refmac_out.pdb"
  atom_selection = None
  mtz_in = "/data1/filo/results/SWAMP_benchmarking/3zux_phases.mtz"
  map_in = None
  labin = ""
  offset_pdb = "offset.pdb"
  resolution = 0
  use_only_refl_present_in_mtz = False
  scale = False
  split_conformers = False
  fix_xyz = False
  fix_rad_max = False
  rad_max = None
  any_offset = False
  chain_type = *PROTEIN DNA RNA
  temp_dir = "temp_dir"
  output_dir = ""
  gui_output_dir = None
  verbose = True
  quick = False
  raise_sorry = False
  debug = False
  dry_run = False
  resolve_command_list = None
  job_title = None
}
Get_cc_mtz_pdb: correlation of map and model allowing origin offsets
Copied as PDB temp_dir/TEMP_refmac_out.pdb 
Map from:  TEMP_3zux_phases.mtz  using labin  FP=FWT PHIB=PHWT 
Model from:  TEMP_refmac_out.pdb
LABIN LINE:  FP=FWT PHIB=PHWT 
Offsetting  TEMP_refmac_out.pdb  to match  TEMP_3zux_phases.mtz
Getting FC from self.pdb...
FC is in  temp_dir/resolve.mtz
Offset pdb file is in  offset.pdb
Getting CC of  TEMP_refmac_out.pdb  (offset to  offset.pdb ) to match  TEMP_3zux_phases.mtz
Detailed analysis of correlation of map and model are in:  cc.log
Offset PDB file is in offset.pdb

Correlation in region of model:  0.393 ...overall:  0.139
overall CC:  0.139
local CC:  0.393

Citations for get_cc_mtz_pdb:

Adams PD, Afonine PV, Bunkoczi G, Chen VB, Davis IW, Echols N, Headd JJ, Hung
LW, Kapral GJ, Grosse-Kunstleve RW, McCoy AJ, Moriarty NW, Oeffner R, Read RJ,
Richardson DC, Richardson JS, Terwilliger TC, Zwart PH. (2010) PHENIX: a
comprehensive Python-based system for macromolecular structure solution. Acta
Cryst. D66:213-221."""

        parser = PhenixParser(stdout=stdout_contents)
        parser.parse()

        self.assertEqual('0.139', parser.overall_CC)
        self.assertEqual('0.393', parser.local_CC)
        self.assertTupleEqual(('0.393', '0.139'), parser.summary)


if __name__ == '__main__':
    unittest.main()
