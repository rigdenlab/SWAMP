import os
import unittest
import collections
from swamp.utils import create_tempfile
from swamp.parsers.pdbtmxmlparser import PdbtmXmlParser


class PdbtmParserTestCase(unittest.TestCase):

    def test_1(self):
        file_contents = """<?xml version="1.0" encoding="iso-8859-1"?>
<pdbtm xmlns="http://pdbtm.enzim.hu" xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance" xsi:schemaLocation="http://pdbtm.enzim.hu/data/pdbtm.xsd pdbtm.xsd" ID="6e7p" TMP="yes">
  <COPYRIGHT>
     All  information, data  and  files are copyright.  PDBTM database is
     produced in the  Institute of  Enzymology,  Budapest, Hungary. There
     are no restrictions on  its  use by  non-profit institutions as long
     as its content  is  in  no  way  modified  and this statement is not
     removed from entries. Usage by and for  commercial entities requires
     a license agreement (send an email to pdbtm at enzim dot hu).
  </COPYRIGHT>
  <CREATE_DATE>2019-04-12</CREATE_DATE>
  <MODIFICATION>
    <DATE>2019-04-12</DATE>
    <DESCR>Interfacial helix detected.</DESCR>
  </MODIFICATION>
  <RAWRES>
    <TMRES>65.90</TMRES>
    <TMTYPE>Tm_Alpha</TMTYPE>
    <SPRES>Unknown</SPRES>
    <PDBKWRES>yes</PDBKWRES>
  </RAWRES>
  <MEMBRANE>
    <NORMAL X="0.01583958" Y="0.02306218" Z="17.99997902"/>
    <TMATRIX>
      <ROWX X="1.00000000" Y="0.00000000" Z="0.00000000" T="-160.51856995"/>
      <ROWY X="0.00000000" Y="1.00000000" Z="0.00000000" T="-160.51933289"/>
      <ROWZ X="0.00000000" Y="0.00000000" Z="1.00000000" T="-162.26676941"/>
    </TMATRIX>
  </MEMBRANE>
  <CHAIN CHAINID="A" NUM_TM="6" TYPE="alpha">
    <SEQ>
      MTAPAGPRGS ETERLLTPNP GYGTQAGPSP APPTPPEEED LRRRLKYFFM 
      SPCDKFRAKG RKPCKLMLQV VKILVVTVQL ILFGLSNQLA VTFREENTIA 
      FRHLFLLGYS DGADDTFAAY TREQLYQAIF HAVDQYLALP DVSLGRYAYV 
      RGGGDPWTNG SGLALCQRYY HRGHVDPAND TFDIDPMVVT DCIQVDPPER 
      PPPPPSDDLT LLESSSSYKN LTLKFHKLVN VTIHFRLKTI NLQSLINNEI 
      PDCYTFSVLI TFDNKAHSGR IPISLETQAH IQECKHPSVF QHGDNSFRLL 
      FDVVVILTCS LSFLLCARSL LRGFLLQNEF VGFMWRQRGR VISLWERLEF 
      VNGWYILLVT SDVLTISGTI MKIGIEAKNL ASYDVCSILL GTSTLLVWVG 
      VIRYLTFFHN YNILIATLRV ALPSVMRFCC CVAVIYLGYC FCGWIVLGPY 
      HVKFRSLSMV SECLFSLING DDMFVTFAAM QAQQGRSSLV WLFSQLYLYS 
      FISLFIYMVL SLFIALITGA YDTIKHPGGA GAEESELQAY IAQCQDSPTS 
      GKFRRGSGSA CSLLCCCGRD PSEEHSLLVN 
    </SEQ>
    <REGION seq_beg="1" pdb_beg="1"  seq_end="37" pdb_end="37"  type="U"/>
    <REGION seq_beg="38" pdb_beg="38"  seq_end="62" pdb_end="62"  type="2"/>
    <REGION seq_beg="63" pdb_beg="63"  seq_end="84" pdb_end="84"  type="H"/>
    <REGION seq_beg="85" pdb_beg="85"  seq_end="201" pdb_end="201"  type="1"/>
    <REGION seq_beg="202" pdb_beg="202"  seq_end="215" pdb_end="215"  type="U"/>
    <REGION seq_beg="216" pdb_beg="216"  seq_end="297" pdb_end="297"  type="1"/>
    <REGION seq_beg="298" pdb_beg="298"  seq_end="324" pdb_end="324"  type="H"/>
    <REGION seq_beg="325" pdb_beg="325"  seq_end="351" pdb_end="351"  type="2"/>
    <REGION seq_beg="352" pdb_beg="352"  seq_end="373" pdb_end="373"  type="H"/>
    <REGION seq_beg="374" pdb_beg="374"  seq_end="384" pdb_end="384"  type="1"/>
    <REGION seq_beg="385" pdb_beg="385"  seq_end="408" pdb_end="408"  type="H"/>
    <REGION seq_beg="409" pdb_beg="409"  seq_end="418" pdb_end="418"  type="2"/>
    <REGION seq_beg="419" pdb_beg="419"  seq_end="447" pdb_end="447"  type="H"/>
    <REGION seq_beg="448" pdb_beg="448"  seq_end="456" pdb_end="456"  type="1"/>
    <REGION seq_beg="457" pdb_beg="457"  seq_end="477" pdb_end="477"  type="L"/>
    <REGION seq_beg="478" pdb_beg="478"  seq_end="489" pdb_end="489"  type="1"/>
    <REGION seq_beg="490" pdb_beg="490"  seq_end="518" pdb_end="518"  type="H"/>
    <REGION seq_beg="519" pdb_beg="519"  seq_end="526" pdb_end="526"  type="2"/>
    <REGION seq_beg="527" pdb_beg="527"  seq_end="580" pdb_end="580"  type="U"/>
  </CHAIN>
  <CHAIN CHAINID="B" NUM_TM="6" TYPE="alpha">
    <SEQ>
      MTAPAGPRGS ETERLLTPNP GYGTQAGPSP APPTPPEEED LRRRLKYFFM 
      SPCDKFRAKG RKPCKLMLQV VKILVVTVQL ILFGLSNQLA VTFREENTIA 
      FRHLFLLGYS DGADDTFAAY TREQLYQAIF HAVDQYLALP DVSLGRYAYV 
      RGGGDPWTNG SGLALCQRYY HRGHVDPAND TFDIDPMVVT DCIQVDPPER 
      PPPPPSDDLT LLESSSSYKN LTLKFHKLVN VTIHFRLKTI NLQSLINNEI 
      PDCYTFSVLI TFDNKAHSGR IPISLETQAH IQECKHPSVF QHGDNSFRLL 
      FDVVVILTCS LSFLLCARSL LRGFLLQNEF VGFMWRQRGR VISLWERLEF 
      VNGWYILLVT SDVLTISGTI MKIGIEAKNL ASYDVCSILL GTSTLLVWVG 
      VIRYLTFFHN YNILIATLRV ALPSVMRFCC CVAVIYLGYC FCGWIVLGPY 
      HVKFRSLSMV SECLFSLING DDMFVTFAAM QAQQGRSSLV WLFSQLYLYS 
      FISLFIYMVL SLFIALITGA YDTIKHPGGA GAEESELQAY IAQCQDSPTS 
      GKFRRGSGSA CSLLCCCGRD PSEEHSLLVN 
    </SEQ>
    <REGION seq_beg="1" pdb_beg="1"  seq_end="37" pdb_end="37"  type="U"/>
    <REGION seq_beg="38" pdb_beg="38"  seq_end="62" pdb_end="62"  type="2"/>
    <REGION seq_beg="63" pdb_beg="63"  seq_end="84" pdb_end="84"  type="H"/>
    <REGION seq_beg="85" pdb_beg="85"  seq_end="201" pdb_end="201"  type="1"/>
    <REGION seq_beg="202" pdb_beg="202"  seq_end="215" pdb_end="215"  type="U"/>
    <REGION seq_beg="216" pdb_beg="216"  seq_end="297" pdb_end="297"  type="1"/>
    <REGION seq_beg="298" pdb_beg="298"  seq_end="324" pdb_end="324"  type="H"/>
    <REGION seq_beg="325" pdb_beg="325"  seq_end="351" pdb_end="351"  type="2"/>
    <REGION seq_beg="352" pdb_beg="352"  seq_end="373" pdb_end="373"  type="H"/>
    <REGION seq_beg="374" pdb_beg="374"  seq_end="384" pdb_end="384"  type="1"/>
    <REGION seq_beg="385" pdb_beg="385"  seq_end="408" pdb_end="408"  type="H"/>
    <REGION seq_beg="409" pdb_beg="409"  seq_end="418" pdb_end="418"  type="2"/>
    <REGION seq_beg="419" pdb_beg="419"  seq_end="448" pdb_end="448"  type="H"/>
    <REGION seq_beg="449" pdb_beg="449"  seq_end="456" pdb_end="456"  type="1"/>
    <REGION seq_beg="457" pdb_beg="457"  seq_end="477" pdb_end="477"  type="L"/>
    <REGION seq_beg="478" pdb_beg="478"  seq_end="489" pdb_end="489"  type="1"/>
    <REGION seq_beg="490" pdb_beg="490"  seq_end="518" pdb_end="518"  type="H"/>
    <REGION seq_beg="519" pdb_beg="519"  seq_end="526" pdb_end="526"  type="2"/>
    <REGION seq_beg="527" pdb_beg="527"  seq_end="580" pdb_end="580"  type="U"/>
  </CHAIN>
  <CHAIN CHAINID="C" NUM_TM="6" TYPE="alpha">
    <SEQ>
      MTAPAGPRGS ETERLLTPNP GYGTQAGPSP APPTPPEEED LRRRLKYFFM 
      SPCDKFRAKG RKPCKLMLQV VKILVVTVQL ILFGLSNQLA VTFREENTIA 
      FRHLFLLGYS DGADDTFAAY TREQLYQAIF HAVDQYLALP DVSLGRYAYV 
      RGGGDPWTNG SGLALCQRYY HRGHVDPAND TFDIDPMVVT DCIQVDPPER 
      PPPPPSDDLT LLESSSSYKN LTLKFHKLVN VTIHFRLKTI NLQSLINNEI 
      PDCYTFSVLI TFDNKAHSGR IPISLETQAH IQECKHPSVF QHGDNSFRLL 
      FDVVVILTCS LSFLLCARSL LRGFLLQNEF VGFMWRQRGR VISLWERLEF 
      VNGWYILLVT SDVLTISGTI MKIGIEAKNL ASYDVCSILL GTSTLLVWVG 
      VIRYLTFFHN YNILIATLRV ALPSVMRFCC CVAVIYLGYC FCGWIVLGPY 
      HVKFRSLSMV SECLFSLING DDMFVTFAAM QAQQGRSSLV WLFSQLYLYS 
      FISLFIYMVL SLFIALITGA YDTIKHPGGA GAEESELQAY IAQCQDSPTS 
      GKFRRGSGSA CSLLCCCGRD PSEEHSLLVN 
    </SEQ>
    <REGION seq_beg="1" pdb_beg="1"  seq_end="37" pdb_end="37"  type="U"/>
    <REGION seq_beg="38" pdb_beg="38"  seq_end="62" pdb_end="62"  type="2"/>
    <REGION seq_beg="63" pdb_beg="63"  seq_end="84" pdb_end="84"  type="H"/>
    <REGION seq_beg="85" pdb_beg="85"  seq_end="201" pdb_end="201"  type="1"/>
    <REGION seq_beg="202" pdb_beg="202"  seq_end="215" pdb_end="215"  type="U"/>
    <REGION seq_beg="216" pdb_beg="216"  seq_end="297" pdb_end="297"  type="1"/>
    <REGION seq_beg="298" pdb_beg="298"  seq_end="324" pdb_end="324"  type="H"/>
    <REGION seq_beg="325" pdb_beg="325"  seq_end="351" pdb_end="351"  type="2"/>
    <REGION seq_beg="352" pdb_beg="352"  seq_end="373" pdb_end="373"  type="H"/>
    <REGION seq_beg="374" pdb_beg="374"  seq_end="384" pdb_end="384"  type="1"/>
    <REGION seq_beg="385" pdb_beg="385"  seq_end="408" pdb_end="408"  type="H"/>
    <REGION seq_beg="409" pdb_beg="409"  seq_end="418" pdb_end="418"  type="2"/>
    <REGION seq_beg="419" pdb_beg="419"  seq_end="448" pdb_end="448"  type="H"/>
    <REGION seq_beg="449" pdb_beg="449"  seq_end="456" pdb_end="456"  type="1"/>
    <REGION seq_beg="457" pdb_beg="457"  seq_end="477" pdb_end="477"  type="L"/>
    <REGION seq_beg="478" pdb_beg="478"  seq_end="489" pdb_end="489"  type="1"/>
    <REGION seq_beg="490" pdb_beg="490"  seq_end="518" pdb_end="518"  type="H"/>
    <REGION seq_beg="519" pdb_beg="519"  seq_end="526" pdb_end="526"  type="2"/>
    <REGION seq_beg="527" pdb_beg="527"  seq_end="580" pdb_end="580"  type="U"/>
  </CHAIN>
  <CHAIN CHAINID="D" NUM_TM="6" TYPE="alpha">
    <SEQ>
      MTAPAGPRGS ETERLLTPNP GYGTQAGPSP APPTPPEEED LRRRLKYFFM 
      SPCDKFRAKG RKPCKLMLQV VKILVVTVQL ILFGLSNQLA VTFREENTIA 
      FRHLFLLGYS DGADDTFAAY TREQLYQAIF HAVDQYLALP DVSLGRYAYV 
      RGGGDPWTNG SGLALCQRYY HRGHVDPAND TFDIDPMVVT DCIQVDPPER 
      PPPPPSDDLT LLESSSSYKN LTLKFHKLVN VTIHFRLKTI NLQSLINNEI 
      PDCYTFSVLI TFDNKAHSGR IPISLETQAH IQECKHPSVF QHGDNSFRLL 
      FDVVVILTCS LSFLLCARSL LRGFLLQNEF VGFMWRQRGR VISLWERLEF 
      VNGWYILLVT SDVLTISGTI MKIGIEAKNL ASYDVCSILL GTSTLLVWVG 
      VIRYLTFFHN YNILIATLRV ALPSVMRFCC CVAVIYLGYC FCGWIVLGPY 
      HVKFRSLSMV SECLFSLING DDMFVTFAAM QAQQGRSSLV WLFSQLYLYS 
      FISLFIYMVL SLFIALITGA YDTIKHPGGA GAEESELQAY IAQCQDSPTS 
      GKFRRGSGSA CSLLCCCGRD PSEEHSLLVN 
    </SEQ>
    <REGION seq_beg="1" pdb_beg="1"  seq_end="37" pdb_end="37"  type="U"/>
    <REGION seq_beg="38" pdb_beg="38"  seq_end="62" pdb_end="62"  type="2"/>
    <REGION seq_beg="63" pdb_beg="63"  seq_end="84" pdb_end="84"  type="H"/>
    <REGION seq_beg="85" pdb_beg="85"  seq_end="201" pdb_end="201"  type="1"/>
    <REGION seq_beg="202" pdb_beg="202"  seq_end="215" pdb_end="215"  type="U"/>
    <REGION seq_beg="216" pdb_beg="216"  seq_end="297" pdb_end="297"  type="1"/>
    <REGION seq_beg="298" pdb_beg="298"  seq_end="324" pdb_end="324"  type="H"/>
    <REGION seq_beg="325" pdb_beg="325"  seq_end="351" pdb_end="351"  type="2"/>
    <REGION seq_beg="352" pdb_beg="352"  seq_end="373" pdb_end="373"  type="H"/>
    <REGION seq_beg="374" pdb_beg="374"  seq_end="384" pdb_end="384"  type="1"/>
    <REGION seq_beg="385" pdb_beg="385"  seq_end="408" pdb_end="408"  type="H"/>
    <REGION seq_beg="409" pdb_beg="409"  seq_end="418" pdb_end="418"  type="2"/>
    <REGION seq_beg="419" pdb_beg="419"  seq_end="448" pdb_end="448"  type="H"/>
    <REGION seq_beg="449" pdb_beg="449"  seq_end="456" pdb_end="456"  type="1"/>
    <REGION seq_beg="457" pdb_beg="457"  seq_end="477" pdb_end="477"  type="L"/>
    <REGION seq_beg="478" pdb_beg="478"  seq_end="489" pdb_end="489"  type="1"/>
    <REGION seq_beg="490" pdb_beg="490"  seq_end="518" pdb_end="518"  type="H"/>
    <REGION seq_beg="519" pdb_beg="519"  seq_end="526" pdb_end="526"  type="2"/>
    <REGION seq_beg="527" pdb_beg="527"  seq_end="580" pdb_end="580"  type="U"/>
  </CHAIN>
</pdbtm>
"""
        AnnotationInfo = collections.namedtuple("AnnotationInfo",
                                                ["pdb_start", "pdb_stop", "seq_start", "seq_end", "type", "chain",
                                                 "index",
                                                 "length", "pdb_region"])
        fname = create_tempfile(content=file_contents)
        self.addCleanup(os.remove, fname)
        parser = PdbtmXmlParser(fname=fname)
        parser.parse()

        ss2_annot = [
            AnnotationInfo(pdb_start=1, pdb_stop=37, seq_start=1, seq_end=37, type='U', chain='A', index=0, length=36,
                           pdb_region=[1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22,
                                       23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37]),
            AnnotationInfo(pdb_start=38, pdb_stop=62, seq_start=38, seq_end=62, type='2', chain='A', index=1, length=24,
                           pdb_region=[38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54, 55, 56, 57,
                                       58, 59, 60, 61, 62]),
            AnnotationInfo(pdb_start=63, pdb_stop=84, seq_start=63, seq_end=84, type='H', chain='A', index=2, length=21,
                           pdb_region=[63, 64, 65, 66, 67, 68, 69, 70, 71, 72, 73, 74, 75, 76, 77, 78, 79, 80, 81, 82,
                                       83, 84]),
            AnnotationInfo(pdb_start=85, pdb_stop=201, seq_start=85, seq_end=201, type='1', chain='A', index=3,
                           length=116,
                           pdb_region=[85, 86, 87, 88, 89, 90, 91, 92, 93, 94, 95, 96, 97, 98, 99, 100, 101, 102, 103,
                                       104, 105, 106, 107, 108, 109, 110, 111, 112, 113, 114, 115, 116, 117, 118, 119,
                                       120, 121, 122, 123, 124, 125, 126, 127, 128, 129, 130, 131, 132, 133, 134, 135,
                                       136, 137, 138, 139, 140, 141, 142, 143, 144, 145, 146, 147, 148, 149, 150, 151,
                                       152, 153, 154, 155, 156, 157, 158, 159, 160, 161, 162, 163, 164, 165, 166, 167,
                                       168, 169, 170, 171, 172, 173, 174, 175, 176, 177, 178, 179, 180, 181, 182, 183,
                                       184, 185, 186, 187, 188, 189, 190, 191, 192, 193, 194, 195, 196, 197, 198, 199,
                                       200, 201]),
            AnnotationInfo(pdb_start=202, pdb_stop=215, seq_start=202, seq_end=215, type='U', chain='A', index=4,
                           length=13,
                           pdb_region=[202, 203, 204, 205, 206, 207, 208, 209, 210, 211, 212, 213, 214, 215]),
            AnnotationInfo(pdb_start=216, pdb_stop=297, seq_start=216, seq_end=297, type='1', chain='A', index=5,
                           length=81,
                           pdb_region=[216, 217, 218, 219, 220, 221, 222, 223, 224, 225, 226, 227, 228, 229, 230, 231,
                                       232, 233, 234, 235, 236, 237, 238, 239, 240, 241, 242, 243, 244, 245, 246, 247,
                                       248, 249, 250, 251, 252, 253, 254, 255, 256, 257, 258, 259, 260, 261, 262, 263,
                                       264, 265, 266, 267, 268, 269, 270, 271, 272, 273, 274, 275, 276, 277, 278, 279,
                                       280, 281, 282, 283, 284, 285, 286, 287, 288, 289, 290, 291, 292, 293, 294, 295,
                                       296, 297]),
            AnnotationInfo(pdb_start=298, pdb_stop=324, seq_start=298, seq_end=324, type='H', chain='A', index=6,
                           length=26,
                           pdb_region=[298, 299, 300, 301, 302, 303, 304, 305, 306, 307, 308, 309, 310, 311, 312, 313,
                                       314, 315, 316, 317, 318, 319, 320, 321, 322, 323, 324]),
            AnnotationInfo(pdb_start=325, pdb_stop=351, seq_start=325, seq_end=351, type='2', chain='A', index=7,
                           length=26,
                           pdb_region=[325, 326, 327, 328, 329, 330, 331, 332, 333, 334, 335, 336, 337, 338, 339, 340,
                                       341, 342, 343, 344, 345, 346, 347, 348, 349, 350, 351]),
            AnnotationInfo(pdb_start=352, pdb_stop=373, seq_start=352, seq_end=373, type='H', chain='A', index=8,
                           length=21,
                           pdb_region=[352, 353, 354, 355, 356, 357, 358, 359, 360, 361, 362, 363, 364, 365, 366, 367,
                                       368, 369, 370, 371, 372, 373]),
            AnnotationInfo(pdb_start=374, pdb_stop=384, seq_start=374, seq_end=384, type='1', chain='A', index=9,
                           length=10, pdb_region=[374, 375, 376, 377, 378, 379, 380, 381, 382, 383, 384]),
            AnnotationInfo(pdb_start=385, pdb_stop=408, seq_start=385, seq_end=408, type='H', chain='A', index=10,
                           length=23,
                           pdb_region=[385, 386, 387, 388, 389, 390, 391, 392, 393, 394, 395, 396, 397, 398, 399, 400,
                                       401, 402, 403, 404, 405, 406, 407, 408]),
            AnnotationInfo(pdb_start=409, pdb_stop=418, seq_start=409, seq_end=418, type='2', chain='A', index=11,
                           length=9, pdb_region=[409, 410, 411, 412, 413, 414, 415, 416, 417, 418]),
            AnnotationInfo(pdb_start=419, pdb_stop=447, seq_start=419, seq_end=447, type='H', chain='A', index=12,
                           length=28,
                           pdb_region=[419, 420, 421, 422, 423, 424, 425, 426, 427, 428, 429, 430, 431, 432, 433, 434,
                                       435, 436, 437, 438, 439, 440, 441, 442, 443, 444, 445, 446, 447]),
            AnnotationInfo(pdb_start=448, pdb_stop=456, seq_start=448, seq_end=456, type='1', chain='A', index=13,
                           length=8, pdb_region=[448, 449, 450, 451, 452, 453, 454, 455, 456]),
            AnnotationInfo(pdb_start=457, pdb_stop=477, seq_start=457, seq_end=477, type='L', chain='A', index=14,
                           length=20,
                           pdb_region=[457, 458, 459, 460, 461, 462, 463, 464, 465, 466, 467, 468, 469, 470, 471, 472,
                                       473, 474, 475, 476, 477]),
            AnnotationInfo(pdb_start=478, pdb_stop=489, seq_start=478, seq_end=489, type='1', chain='A', index=15,
                           length=11, pdb_region=[478, 479, 480, 481, 482, 483, 484, 485, 486, 487, 488, 489]),
            AnnotationInfo(pdb_start=490, pdb_stop=518, seq_start=490, seq_end=518, type='H', chain='A', index=16,
                           length=28,
                           pdb_region=[490, 491, 492, 493, 494, 495, 496, 497, 498, 499, 500, 501, 502, 503, 504, 505,
                                       506, 507, 508, 509, 510, 511, 512, 513, 514, 515, 516, 517, 518]),
            AnnotationInfo(pdb_start=519, pdb_stop=526, seq_start=519, seq_end=526, type='2', chain='A', index=17,
                           length=7, pdb_region=[519, 520, 521, 522, 523, 524, 525, 526]),
            AnnotationInfo(pdb_start=527, pdb_stop=580, seq_start=527, seq_end=580, type='U', chain='A', index=18,
                           length=53,
                           pdb_region=[527, 528, 529, 530, 531, 532, 533, 534, 535, 536, 537, 538, 539, 540, 541, 542,
                                       543, 544, 545, 546, 547, 548, 549, 550, 551, 552, 553, 554, 555, 556, 557, 558,
                                       559, 560, 561, 562, 563, 564, 565, 566, 567, 568, 569, 570, 571, 572, 573, 574,
                                       575, 576, 577, 578, 579, 580]),
            AnnotationInfo(pdb_start=1, pdb_stop=37, seq_start=1, seq_end=37, type='U', chain='B', index=0, length=36,
                           pdb_region=[1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22,
                                       23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37]),
            AnnotationInfo(pdb_start=38, pdb_stop=62, seq_start=38, seq_end=62, type='2', chain='B', index=1, length=24,
                           pdb_region=[38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54, 55, 56, 57,
                                       58, 59, 60, 61, 62]),
            AnnotationInfo(pdb_start=63, pdb_stop=84, seq_start=63, seq_end=84, type='H', chain='B', index=2, length=21,
                           pdb_region=[63, 64, 65, 66, 67, 68, 69, 70, 71, 72, 73, 74, 75, 76, 77, 78, 79, 80, 81, 82,
                                       83, 84]),
            AnnotationInfo(pdb_start=85, pdb_stop=201, seq_start=85, seq_end=201, type='1', chain='B', index=3,
                           length=116,
                           pdb_region=[85, 86, 87, 88, 89, 90, 91, 92, 93, 94, 95, 96, 97, 98, 99, 100, 101, 102, 103,
                                       104, 105, 106, 107, 108, 109, 110, 111, 112, 113, 114, 115, 116, 117, 118, 119,
                                       120, 121, 122, 123, 124, 125, 126, 127, 128, 129, 130, 131, 132, 133, 134, 135,
                                       136, 137, 138, 139, 140, 141, 142, 143, 144, 145, 146, 147, 148, 149, 150, 151,
                                       152, 153, 154, 155, 156, 157, 158, 159, 160, 161, 162, 163, 164, 165, 166, 167,
                                       168, 169, 170, 171, 172, 173, 174, 175, 176, 177, 178, 179, 180, 181, 182, 183,
                                       184, 185, 186, 187, 188, 189, 190, 191, 192, 193, 194, 195, 196, 197, 198, 199,
                                       200, 201]),
            AnnotationInfo(pdb_start=202, pdb_stop=215, seq_start=202, seq_end=215, type='U', chain='B', index=4,
                           length=13,
                           pdb_region=[202, 203, 204, 205, 206, 207, 208, 209, 210, 211, 212, 213, 214, 215]),
            AnnotationInfo(pdb_start=216, pdb_stop=297, seq_start=216, seq_end=297, type='1', chain='B', index=5,
                           length=81,
                           pdb_region=[216, 217, 218, 219, 220, 221, 222, 223, 224, 225, 226, 227, 228, 229, 230, 231,
                                       232, 233, 234, 235, 236, 237, 238, 239, 240, 241, 242, 243, 244, 245, 246, 247,
                                       248, 249, 250, 251, 252, 253, 254, 255, 256, 257, 258, 259, 260, 261, 262, 263,
                                       264, 265, 266, 267, 268, 269, 270, 271, 272, 273, 274, 275, 276, 277, 278, 279,
                                       280, 281, 282, 283, 284, 285, 286, 287, 288, 289, 290, 291, 292, 293, 294, 295,
                                       296, 297]),
            AnnotationInfo(pdb_start=298, pdb_stop=324, seq_start=298, seq_end=324, type='H', chain='B', index=6,
                           length=26,
                           pdb_region=[298, 299, 300, 301, 302, 303, 304, 305, 306, 307, 308, 309, 310, 311, 312, 313,
                                       314, 315, 316, 317, 318, 319, 320, 321, 322, 323, 324]),
            AnnotationInfo(pdb_start=325, pdb_stop=351, seq_start=325, seq_end=351, type='2', chain='B', index=7,
                           length=26,
                           pdb_region=[325, 326, 327, 328, 329, 330, 331, 332, 333, 334, 335, 336, 337, 338, 339, 340,
                                       341, 342, 343, 344, 345, 346, 347, 348, 349, 350, 351]),
            AnnotationInfo(pdb_start=352, pdb_stop=373, seq_start=352, seq_end=373, type='H', chain='B', index=8,
                           length=21,
                           pdb_region=[352, 353, 354, 355, 356, 357, 358, 359, 360, 361, 362, 363, 364, 365, 366, 367,
                                       368, 369, 370, 371, 372, 373]),
            AnnotationInfo(pdb_start=374, pdb_stop=384, seq_start=374, seq_end=384, type='1', chain='B', index=9,
                           length=10, pdb_region=[374, 375, 376, 377, 378, 379, 380, 381, 382, 383, 384]),
            AnnotationInfo(pdb_start=385, pdb_stop=408, seq_start=385, seq_end=408, type='H', chain='B', index=10,
                           length=23,
                           pdb_region=[385, 386, 387, 388, 389, 390, 391, 392, 393, 394, 395, 396, 397, 398, 399, 400,
                                       401, 402, 403, 404, 405, 406, 407, 408]),
            AnnotationInfo(pdb_start=409, pdb_stop=418, seq_start=409, seq_end=418, type='2', chain='B', index=11,
                           length=9, pdb_region=[409, 410, 411, 412, 413, 414, 415, 416, 417, 418]),
            AnnotationInfo(pdb_start=419, pdb_stop=448, seq_start=419, seq_end=448, type='H', chain='B', index=12,
                           length=29,
                           pdb_region=[419, 420, 421, 422, 423, 424, 425, 426, 427, 428, 429, 430, 431, 432, 433, 434,
                                       435, 436, 437, 438, 439, 440, 441, 442, 443, 444, 445, 446, 447, 448]),
            AnnotationInfo(pdb_start=449, pdb_stop=456, seq_start=449, seq_end=456, type='1', chain='B', index=13,
                           length=7, pdb_region=[449, 450, 451, 452, 453, 454, 455, 456]),
            AnnotationInfo(pdb_start=457, pdb_stop=477, seq_start=457, seq_end=477, type='L', chain='B', index=14,
                           length=20,
                           pdb_region=[457, 458, 459, 460, 461, 462, 463, 464, 465, 466, 467, 468, 469, 470, 471, 472,
                                       473, 474, 475, 476, 477]),
            AnnotationInfo(pdb_start=478, pdb_stop=489, seq_start=478, seq_end=489, type='1', chain='B', index=15,
                           length=11, pdb_region=[478, 479, 480, 481, 482, 483, 484, 485, 486, 487, 488, 489]),
            AnnotationInfo(pdb_start=490, pdb_stop=518, seq_start=490, seq_end=518, type='H', chain='B', index=16,
                           length=28,
                           pdb_region=[490, 491, 492, 493, 494, 495, 496, 497, 498, 499, 500, 501, 502, 503, 504, 505,
                                       506, 507, 508, 509, 510, 511, 512, 513, 514, 515, 516, 517, 518]),
            AnnotationInfo(pdb_start=519, pdb_stop=526, seq_start=519, seq_end=526, type='2', chain='B', index=17,
                           length=7, pdb_region=[519, 520, 521, 522, 523, 524, 525, 526]),
            AnnotationInfo(pdb_start=527, pdb_stop=580, seq_start=527, seq_end=580, type='U', chain='B', index=18,
                           length=53,
                           pdb_region=[527, 528, 529, 530, 531, 532, 533, 534, 535, 536, 537, 538, 539, 540, 541, 542,
                                       543, 544, 545, 546, 547, 548, 549, 550, 551, 552, 553, 554, 555, 556, 557, 558,
                                       559, 560, 561, 562, 563, 564, 565, 566, 567, 568, 569, 570, 571, 572, 573, 574,
                                       575, 576, 577, 578, 579, 580]),
            AnnotationInfo(pdb_start=1, pdb_stop=37, seq_start=1, seq_end=37, type='U', chain='C', index=0, length=36,
                           pdb_region=[1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22,
                                       23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37]),
            AnnotationInfo(pdb_start=38, pdb_stop=62, seq_start=38, seq_end=62, type='2', chain='C', index=1, length=24,
                           pdb_region=[38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54, 55, 56, 57,
                                       58, 59, 60, 61, 62]),
            AnnotationInfo(pdb_start=63, pdb_stop=84, seq_start=63, seq_end=84, type='H', chain='C', index=2, length=21,
                           pdb_region=[63, 64, 65, 66, 67, 68, 69, 70, 71, 72, 73, 74, 75, 76, 77, 78, 79, 80, 81, 82,
                                       83, 84]),
            AnnotationInfo(pdb_start=85, pdb_stop=201, seq_start=85, seq_end=201, type='1', chain='C', index=3,
                           length=116,
                           pdb_region=[85, 86, 87, 88, 89, 90, 91, 92, 93, 94, 95, 96, 97, 98, 99, 100, 101, 102, 103,
                                       104, 105, 106, 107, 108, 109, 110, 111, 112, 113, 114, 115, 116, 117, 118, 119,
                                       120, 121, 122, 123, 124, 125, 126, 127, 128, 129, 130, 131, 132, 133, 134, 135,
                                       136, 137, 138, 139, 140, 141, 142, 143, 144, 145, 146, 147, 148, 149, 150, 151,
                                       152, 153, 154, 155, 156, 157, 158, 159, 160, 161, 162, 163, 164, 165, 166, 167,
                                       168, 169, 170, 171, 172, 173, 174, 175, 176, 177, 178, 179, 180, 181, 182, 183,
                                       184, 185, 186, 187, 188, 189, 190, 191, 192, 193, 194, 195, 196, 197, 198, 199,
                                       200, 201]),
            AnnotationInfo(pdb_start=202, pdb_stop=215, seq_start=202, seq_end=215, type='U', chain='C', index=4,
                           length=13,
                           pdb_region=[202, 203, 204, 205, 206, 207, 208, 209, 210, 211, 212, 213, 214, 215]),
            AnnotationInfo(pdb_start=216, pdb_stop=297, seq_start=216, seq_end=297, type='1', chain='C', index=5,
                           length=81,
                           pdb_region=[216, 217, 218, 219, 220, 221, 222, 223, 224, 225, 226, 227, 228, 229, 230, 231,
                                       232, 233, 234, 235, 236, 237, 238, 239, 240, 241, 242, 243, 244, 245, 246, 247,
                                       248, 249, 250, 251, 252, 253, 254, 255, 256, 257, 258, 259, 260, 261, 262, 263,
                                       264, 265, 266, 267, 268, 269, 270, 271, 272, 273, 274, 275, 276, 277, 278, 279,
                                       280, 281, 282, 283, 284, 285, 286, 287, 288, 289, 290, 291, 292, 293, 294, 295,
                                       296, 297]),
            AnnotationInfo(pdb_start=298, pdb_stop=324, seq_start=298, seq_end=324, type='H', chain='C', index=6,
                           length=26,
                           pdb_region=[298, 299, 300, 301, 302, 303, 304, 305, 306, 307, 308, 309, 310, 311, 312, 313,
                                       314, 315, 316, 317, 318, 319, 320, 321, 322, 323, 324]),
            AnnotationInfo(pdb_start=325, pdb_stop=351, seq_start=325, seq_end=351, type='2', chain='C', index=7,
                           length=26,
                           pdb_region=[325, 326, 327, 328, 329, 330, 331, 332, 333, 334, 335, 336, 337, 338, 339, 340,
                                       341, 342, 343, 344, 345, 346, 347, 348, 349, 350, 351]),
            AnnotationInfo(pdb_start=352, pdb_stop=373, seq_start=352, seq_end=373, type='H', chain='C', index=8,
                           length=21,
                           pdb_region=[352, 353, 354, 355, 356, 357, 358, 359, 360, 361, 362, 363, 364, 365, 366, 367,
                                       368, 369, 370, 371, 372, 373]),
            AnnotationInfo(pdb_start=374, pdb_stop=384, seq_start=374, seq_end=384, type='1', chain='C', index=9,
                           length=10, pdb_region=[374, 375, 376, 377, 378, 379, 380, 381, 382, 383, 384]),
            AnnotationInfo(pdb_start=385, pdb_stop=408, seq_start=385, seq_end=408, type='H', chain='C', index=10,
                           length=23,
                           pdb_region=[385, 386, 387, 388, 389, 390, 391, 392, 393, 394, 395, 396, 397, 398, 399, 400,
                                       401, 402, 403, 404, 405, 406, 407, 408]),
            AnnotationInfo(pdb_start=409, pdb_stop=418, seq_start=409, seq_end=418, type='2', chain='C', index=11,
                           length=9, pdb_region=[409, 410, 411, 412, 413, 414, 415, 416, 417, 418]),
            AnnotationInfo(pdb_start=419, pdb_stop=448, seq_start=419, seq_end=448, type='H', chain='C', index=12,
                           length=29,
                           pdb_region=[419, 420, 421, 422, 423, 424, 425, 426, 427, 428, 429, 430, 431, 432, 433, 434,
                                       435, 436, 437, 438, 439, 440, 441, 442, 443, 444, 445, 446, 447, 448]),
            AnnotationInfo(pdb_start=449, pdb_stop=456, seq_start=449, seq_end=456, type='1', chain='C', index=13,
                           length=7, pdb_region=[449, 450, 451, 452, 453, 454, 455, 456]),
            AnnotationInfo(pdb_start=457, pdb_stop=477, seq_start=457, seq_end=477, type='L', chain='C', index=14,
                           length=20,
                           pdb_region=[457, 458, 459, 460, 461, 462, 463, 464, 465, 466, 467, 468, 469, 470, 471, 472,
                                       473, 474, 475, 476, 477]),
            AnnotationInfo(pdb_start=478, pdb_stop=489, seq_start=478, seq_end=489, type='1', chain='C', index=15,
                           length=11, pdb_region=[478, 479, 480, 481, 482, 483, 484, 485, 486, 487, 488, 489]),
            AnnotationInfo(pdb_start=490, pdb_stop=518, seq_start=490, seq_end=518, type='H', chain='C', index=16,
                           length=28,
                           pdb_region=[490, 491, 492, 493, 494, 495, 496, 497, 498, 499, 500, 501, 502, 503, 504, 505,
                                       506, 507, 508, 509, 510, 511, 512, 513, 514, 515, 516, 517, 518]),
            AnnotationInfo(pdb_start=519, pdb_stop=526, seq_start=519, seq_end=526, type='2', chain='C', index=17,
                           length=7, pdb_region=[519, 520, 521, 522, 523, 524, 525, 526]),
            AnnotationInfo(pdb_start=527, pdb_stop=580, seq_start=527, seq_end=580, type='U', chain='C', index=18,
                           length=53,
                           pdb_region=[527, 528, 529, 530, 531, 532, 533, 534, 535, 536, 537, 538, 539, 540, 541, 542,
                                       543, 544, 545, 546, 547, 548, 549, 550, 551, 552, 553, 554, 555, 556, 557, 558,
                                       559, 560, 561, 562, 563, 564, 565, 566, 567, 568, 569, 570, 571, 572, 573, 574,
                                       575, 576, 577, 578, 579, 580]),
            AnnotationInfo(pdb_start=1, pdb_stop=37, seq_start=1, seq_end=37, type='U', chain='D', index=0, length=36,
                           pdb_region=[1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22,
                                       23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37]),
            AnnotationInfo(pdb_start=38, pdb_stop=62, seq_start=38, seq_end=62, type='2', chain='D', index=1, length=24,
                           pdb_region=[38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54, 55, 56, 57,
                                       58, 59, 60, 61, 62]),
            AnnotationInfo(pdb_start=63, pdb_stop=84, seq_start=63, seq_end=84, type='H', chain='D', index=2, length=21,
                           pdb_region=[63, 64, 65, 66, 67, 68, 69, 70, 71, 72, 73, 74, 75, 76, 77, 78, 79, 80, 81, 82,
                                       83, 84]),
            AnnotationInfo(pdb_start=85, pdb_stop=201, seq_start=85, seq_end=201, type='1', chain='D', index=3,
                           length=116,
                           pdb_region=[85, 86, 87, 88, 89, 90, 91, 92, 93, 94, 95, 96, 97, 98, 99, 100, 101, 102, 103,
                                       104, 105, 106, 107, 108, 109, 110, 111, 112, 113, 114, 115, 116, 117, 118, 119,
                                       120, 121, 122, 123, 124, 125, 126, 127, 128, 129, 130, 131, 132, 133, 134, 135,
                                       136, 137, 138, 139, 140, 141, 142, 143, 144, 145, 146, 147, 148, 149, 150, 151,
                                       152, 153, 154, 155, 156, 157, 158, 159, 160, 161, 162, 163, 164, 165, 166, 167,
                                       168, 169, 170, 171, 172, 173, 174, 175, 176, 177, 178, 179, 180, 181, 182, 183,
                                       184, 185, 186, 187, 188, 189, 190, 191, 192, 193, 194, 195, 196, 197, 198, 199,
                                       200, 201]),
            AnnotationInfo(pdb_start=202, pdb_stop=215, seq_start=202, seq_end=215, type='U', chain='D', index=4,
                           length=13,
                           pdb_region=[202, 203, 204, 205, 206, 207, 208, 209, 210, 211, 212, 213, 214, 215]),
            AnnotationInfo(pdb_start=216, pdb_stop=297, seq_start=216, seq_end=297, type='1', chain='D', index=5,
                           length=81,
                           pdb_region=[216, 217, 218, 219, 220, 221, 222, 223, 224, 225, 226, 227, 228, 229, 230, 231,
                                       232, 233, 234, 235, 236, 237, 238, 239, 240, 241, 242, 243, 244, 245, 246, 247,
                                       248, 249, 250, 251, 252, 253, 254, 255, 256, 257, 258, 259, 260, 261, 262, 263,
                                       264, 265, 266, 267, 268, 269, 270, 271, 272, 273, 274, 275, 276, 277, 278, 279,
                                       280, 281, 282, 283, 284, 285, 286, 287, 288, 289, 290, 291, 292, 293, 294, 295,
                                       296, 297]),
            AnnotationInfo(pdb_start=298, pdb_stop=324, seq_start=298, seq_end=324, type='H', chain='D', index=6,
                           length=26,
                           pdb_region=[298, 299, 300, 301, 302, 303, 304, 305, 306, 307, 308, 309, 310, 311, 312, 313,
                                       314, 315, 316, 317, 318, 319, 320, 321, 322, 323, 324]),
            AnnotationInfo(pdb_start=325, pdb_stop=351, seq_start=325, seq_end=351, type='2', chain='D', index=7,
                           length=26,
                           pdb_region=[325, 326, 327, 328, 329, 330, 331, 332, 333, 334, 335, 336, 337, 338, 339, 340,
                                       341, 342, 343, 344, 345, 346, 347, 348, 349, 350, 351]),
            AnnotationInfo(pdb_start=352, pdb_stop=373, seq_start=352, seq_end=373, type='H', chain='D', index=8,
                           length=21,
                           pdb_region=[352, 353, 354, 355, 356, 357, 358, 359, 360, 361, 362, 363, 364, 365, 366, 367,
                                       368, 369, 370, 371, 372, 373]),
            AnnotationInfo(pdb_start=374, pdb_stop=384, seq_start=374, seq_end=384, type='1', chain='D', index=9,
                           length=10, pdb_region=[374, 375, 376, 377, 378, 379, 380, 381, 382, 383, 384]),
            AnnotationInfo(pdb_start=385, pdb_stop=408, seq_start=385, seq_end=408, type='H', chain='D', index=10,
                           length=23,
                           pdb_region=[385, 386, 387, 388, 389, 390, 391, 392, 393, 394, 395, 396, 397, 398, 399, 400,
                                       401, 402, 403, 404, 405, 406, 407, 408]),
            AnnotationInfo(pdb_start=409, pdb_stop=418, seq_start=409, seq_end=418, type='2', chain='D', index=11,
                           length=9, pdb_region=[409, 410, 411, 412, 413, 414, 415, 416, 417, 418]),
            AnnotationInfo(pdb_start=419, pdb_stop=448, seq_start=419, seq_end=448, type='H', chain='D', index=12,
                           length=29,
                           pdb_region=[419, 420, 421, 422, 423, 424, 425, 426, 427, 428, 429, 430, 431, 432, 433, 434,
                                       435, 436, 437, 438, 439, 440, 441, 442, 443, 444, 445, 446, 447, 448]),
            AnnotationInfo(pdb_start=449, pdb_stop=456, seq_start=449, seq_end=456, type='1', chain='D', index=13,
                           length=7, pdb_region=[449, 450, 451, 452, 453, 454, 455, 456]),
            AnnotationInfo(pdb_start=457, pdb_stop=477, seq_start=457, seq_end=477, type='L', chain='D', index=14,
                           length=20,
                           pdb_region=[457, 458, 459, 460, 461, 462, 463, 464, 465, 466, 467, 468, 469, 470, 471, 472,
                                       473, 474, 475, 476, 477]),
            AnnotationInfo(pdb_start=478, pdb_stop=489, seq_start=478, seq_end=489, type='1', chain='D', index=15,
                           length=11, pdb_region=[478, 479, 480, 481, 482, 483, 484, 485, 486, 487, 488, 489]),
            AnnotationInfo(pdb_start=490, pdb_stop=518, seq_start=490, seq_end=518, type='H', chain='D', index=16,
                           length=28,
                           pdb_region=[490, 491, 492, 493, 494, 495, 496, 497, 498, 499, 500, 501, 502, 503, 504, 505,
                                       506, 507, 508, 509, 510, 511, 512, 513, 514, 515, 516, 517, 518]),
            AnnotationInfo(pdb_start=519, pdb_stop=526, seq_start=519, seq_end=526, type='2', chain='D', index=17,
                           length=7, pdb_region=[519, 520, 521, 522, 523, 524, 525, 526]),
            AnnotationInfo(pdb_start=527, pdb_stop=580, seq_start=527, seq_end=580, type='U', chain='D', index=18,
                           length=53,
                           pdb_region=[527, 528, 529, 530, 531, 532, 533, 534, 535, 536, 537, 538, 539, 540, 541, 542,
                                       543, 544, 545, 546, 547, 548, 549, 550, 551, 552, 553, 554, 555, 556, 557, 558,
                                       559, 560, 561, 562, 563, 564, 565, 566, 567, 568, 569, 570, 571, 572, 573, 574,
                                       575, 576, 577, 578, 579, 580])
        ]

        self.assertListEqual(ss2_annot, parser.ss2_annotation)


if __name__ == '__main__':
    unittest.main()
