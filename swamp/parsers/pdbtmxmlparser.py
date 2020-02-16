import collections
import xml.etree.ElementTree as ET
from swamp.parsers.parser import Parser


class PdbtmXmlParser(Parser):
    """Class to parse and store pdbtm data contained in xml format files

    :ivar ss2_annotation: List of named tuples containing secondary structure annotations found in the input file
    :type ss2_annotation: list, None

    :example:

    >>> from swamp.parsers import PdbtmXmlParser
    >>> my_parser = PdbtmXmlParser('<fname>')
    >>> my_parser.parse()
    """

    def __init__(self, fname, logger=None):
        self._ss2_annotation = None
        super(PdbtmXmlParser, self).__init__(fname, logger=logger)

    @property
    def summary(self):
        """Abstract property to store a summary of the parsed figures of merit"""
        return None

    @property
    def ss2_annotation(self):
        return self._ss2_annotation

    @ss2_annotation.setter
    def ss2_annotation(self, value):
        self._ss2_annotation = value

    @property
    def _annotation_template(self):
        """Template with the named tuple to contain information about the secondary structure information"""

        return collections.namedtuple("AnnotationInfo",
                                      ["pdb_start", "pdb_stop", "seq_start", "seq_end", "type", "chain", "index",
                                       "length", "pdb_region"])

    def parse(self):
        """Method to parse the input xml file and retrieve the ss2 annotation"""

        if self.error:
            self.logger.warning("Previous errors prevent parsing PDBTM file!")
            return

        tree = ET.parse(self.fname)
        root = tree.getroot()

        self.ss2_annotation = []
        for chain in root.iter("{http://pdbtm.enzim.hu}CHAIN"):
            for idx, region in enumerate(chain.iter("{http://pdbtm.enzim.hu}REGION")):
                current_region = self._annotation_template(index=idx,
                                                           pdb_start=int(region.attrib['pdb_beg']),
                                                           pdb_stop=int(region.attrib['pdb_end']),
                                                           seq_start=int(region.attrib['seq_beg']),
                                                           seq_end=int(region.attrib['seq_end']),
                                                           type=region.attrib['type'],
                                                           chain=chain.attrib["CHAINID"],
                                                           pdb_region=[x for x in range(int(region.attrib['pdb_beg']),
                                                                                        int(region.attrib[
                                                                                                'pdb_end']) + 1)],
                                                           length=(int(region.attrib['seq_end']) - int(
                                                               region.attrib['seq_beg'])))
                self.ss2_annotation.append(current_region)
