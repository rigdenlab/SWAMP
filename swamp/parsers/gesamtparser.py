import numpy as np
import pandas as pd
from swamp.parsers.parser import Parser


class GesamtParser(Parser):
    """Gesamt output parser

    :param str mode: corresponds with :py:attr:`~swamp.wrappers.gesamt.Gesamt.mode` used to create the output to be \
    parsed
    :param str stdout: the stdout to be parsed (default None)
    :param str fname: the file name to be parsed (default None)
    :param `~swamp.logger.swamplogger.SwampLogger` logger: logging interface for the parser (default None)
    :ivar bool error: if True an error has occurred along the process
    :ivar float qscore: qscore as reported by gesamt
    :ivar float rmsd: the obtained rmsd as reported by gesamt
    :ivar float seq_id: sequence identity between the input structures
    :ivar int n_align: number of aligned residues

    :example:

    >>> from swamp.parsers import GesamtParser
    >>> my_parser = GesamtParser('<mode>', '<stdout>', '<fname>')
    >>> my_parser.parse()
    """

    def __init__(self, mode, stdout=None, fname=None, logger=None):
        self.mode = mode
        self.qscore = None
        self.rmsd = None
        self.seq_id = None
        self.n_align = None
        self.hits_df = None
        super(GesamtParser, self).__init__(stdout=stdout, fname=fname, logger=logger)

    @property
    def summary(self):
        """Dataframe with hits found in the archive if :py:attr:`~swmap.parsers.gesamtparser.GesamtParser.mode` is
        'search-archive' otherwise a tuple with all the parsed figures of merit"""

        if self.mode == 'search-archive':
            return self.hits_df
        else:
            return self.qscore, self.rmsd, self.seq_id, self.n_align

    def parse(self):
        """Method to parse :py:attr:`~swamp.parsers.parser.Parser.fname` and store figures of merit"""
        if self.mode == 'search-archive':
            self.parse_hitfile()
        else:
            self.parse_stdout()

    def parse_stdout(self):
        """Method to retrieve qscore, rmsd, sequence identity and no. of aligned residues from \
        :py:attr:`~swamp.parsers.parser.gesamtparser.GesamtParser.stdout`

        :param str stdout: gesamt stdout to be parsed
        :param int n_models: number of models that were used in the structural alignment to generate the provided stdout
        :returns: qscore, rmsd, sequence identity and no. of aligned residues (tuple)
        """

        n_models = 0
        for line in self.stdout.split('\n'):
            if '... reading ' in line:
                n_models += 1

        if n_models == 2:
            qscore_mark = "Q-score"
            rmsd_mark = "RMSD"
            n_align_mark = "Aligned residues"
            seqid_mark = "Sequence Id"
        else:
            qscore_mark = "quality Q"
            rmsd_mark = "r.m.s.d"
            n_align_mark = "Nalign"
            seqid_mark = "SEQ_ID IS NOT FOUND IN MULTIPLE STRCUT. ALIGNMENT"

        self.qscore = np.nan
        self.rmsd = np.nan
        self.n_align = np.nan
        self.seq_id = np.nan
        for line in self.stdout.split("\n"):
            if len(line.split()) != 0 and line.split()[0] != "#":
                if qscore_mark in line and self.qscore is np.nan:
                    self.qscore = float(line.rstrip().lstrip().split(":")[-1].split()[0].rstrip().lstrip())
                elif rmsd_mark in line and self.rmsd is np.nan:
                    self.rmsd = float(line.rstrip().lstrip().split(":")[-1].split()[0].rstrip().lstrip())
                elif n_align_mark in line and self.n_align is np.nan:
                    self.n_align = int(line.rstrip().lstrip().split(":")[-1].split()[0].rstrip().lstrip())
                elif seqid_mark in line and self.seq_id is np.nan:
                    self.seq_id = float(line.rstrip().lstrip().split(":")[-1].split()[0].rstrip().lstrip())

    def parse_hitfile(self):
        """Method to parse a gesamt .hit output file

        :param str fname: file name of the .hit output file
        :returns: a dataframe with the results contained in the hit file (`pandas.Dataframe`)
        """

        self.hits_df = []
        with open(self.fname, "r") as fhandle:
            for line in fhandle:
                if line[0] != "#":
                    line = line[20:].split()
                    self.hits_df.append([line[-6], line[-5], line[-4], line[-3], line[-2], line[-1]])
        self.hits_df = pd.DataFrame(self.hits_df)
        self.hits_df.columns = ["qscore", "rmsd", "seq_id", "n_align", "n_res", "fname"]

    @staticmethod
    def get_pairwise_qscores(stdout):
        """Method to get the pairwise qscores of a given alignmnet between several models in an ensemble

        :param str stdout: gesamt stdout for the command
        :returns: qscores_dict: a dictionary with the pairwise qscores for each of the models in the alignment (dict)
        """

        qscores_dict = {}
        structure_id_dict = {}
        qscores_mark = "(o) pairwise Q-scores"
        file_mark = "... reading file"
        rmsd_mark = "(o) pairwise r.m.s.d."
        is_qscores = False

        for line in stdout.split("\n"):

            # Store file names and structure ids
            if file_mark in line:
                fname = line.split("'")[1]
                structure_id = "S%s" % str(len(qscores_dict.keys()) + 1).zfill(3)
                qscores_dict[fname] = None
                structure_id_dict[structure_id] = fname
            # Qscores will start appearing now
            elif qscores_mark in line:
                is_qscores = True
            # Store the qscore in the dictionary
            elif is_qscores and line.split("|")[0].rstrip().lstrip() in structure_id_dict.keys():
                structure_id = line.split("|")[0].rstrip().lstrip()
                idx = int(structure_id[1:])
                qscores_dict[structure_id_dict[structure_id]] = float(line.split()[idx].rstrip().lstrip())
            # If we reach the rmsd mark, break the loop
            elif rmsd_mark in line:
                break

        return qscores_dict
