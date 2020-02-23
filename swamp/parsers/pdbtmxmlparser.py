import collections
import xml.etree.ElementTree as ET
from swamp.parsers.parser import Parser


class PdbtmXmlParser(Parser):
    """Class to parse and store pdbtm data contained in xml format files

    :param str fname: the file name to be parsed (default None)
    :param `~swamp.logger.swamplogger.SwampLogger` logger: logging interface for the parser (default None)
    :ivar list ss2_annotation: List of named tuples containing secondary structure annotations found in the input file

    :example:

    >>> from swamp.parsers import PdbtmXmlParser
    >>> my_parser = PdbtmXmlParser('<fname>')
    >>> my_parser.parse()
    """

    def __init__(self, fname, logger=None):
        self.ss2_annotation = None
        super(PdbtmXmlParser, self).__init__(fname, logger=logger)

    @property
    def summary(self):
        """Abstract property to store a summary of the parsed figures of merit (not in use)"""
        return None

    @property
    def _annotation_template(self):
        """Named tuple to contain information about the secondary structure information extracted from \
        :py:attr:`~swamp.parsers.parser.fname`"""

        return collections.namedtuple("AnnotationInfo",
                                      ["pdb_start", "pdb_stop", "seq_start", "seq_end", "type", "chain", "index",
                                       "length", "pdb_region"])

    def parse(self):
        """Method to parse :py:attr:`~swamp.parsers.parser.Parser.fname` and store \
        :py:attr:`~swamp.parsers.pdbtmxmlparser.PdbtmXmlParser.ss2 annotation`"""

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
