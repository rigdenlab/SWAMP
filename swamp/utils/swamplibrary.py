import os
import gemmi
import conkit.io
import itertools
import pandas as pd
import swamp.utils as utils
from swamp.wrappers import Gesamt
from swamp.parsers import PdbtmXmlParser
from swamp.logger import SwampLogger


class SwampLibrary(object):
    """Class that implements certain methods to create a SWAMP fragment library and data structures to manage the data
    of interest in the library.

    :param str workdir: the working directory for this instance. Only used if a library will be created.
    :param `~swamp.logger.swamplogger.SwampLogger` logger: logging instance
    :ivar `pandas.DataFrame` rmsd_matrix: square dataframe with the rmsd distance across framgents in the library
    :ivar `pandas.DataFrame` qscore_matrix: square dataframe with the similarity across framgents in the library
    :ivar `pandas.DataFrame` nalign_matrix: square dataframe with the no. of aligned residues between framgents in the \
    library
    :ivar str pdb_library: location of the directory with the pdb files contained in the SWAMP library
    :ivar str pdbtm_svn: location of the pdbtm svn repository
    :ivar str outdir: an output directory for any operation of the :py:obj:`~swamp.utils.swamplibrary.SwampLibrary` \
    instance

    :example:

    >>> from swamp.utils.swamplibrary import SwampLibrary
    >>> my_library = SwampLibrary('<workdir>')
    >>> pdb_code_list = my_library.parse_nr_listfile("/path/to/nr_list")
    >>> my_library.pdbtm_svn = "/path/to/pdbtm_svn"
    >>> my_library.pdb_library = "/path.to/pdb_library"
    >>> my_library.make_library(outdir="/path/to/outdir", pdb_codes=pdb_code_list)
    >>> my_library.all_vs_all_gesamt(outdir="/path/to/outdir", inputdir="/path/to/library", nthreads=1)
    >>> my_library.create_distance_mtx(gesamt_dir="/path/to/gesamt_dir")

    """

    def __init__(self, workdir, logger=None):
        self.workdir = workdir
        self._make_workdir()
        if logger is None:
            self.logger = SwampLogger(__name__)
            self.logger.init(logfile=None, use_console=True, console_level='info')
        else:
            self.logger = logger
        self.qscore_matrix = None
        self.nalign_matrix = None
        self.rmsd_matrix = None
        self.pdb_library = None
        self.pdbtm_svn = None
        self.outdir = None

    # ------------------ Properties ------------------

    @property
    def pdbfiles_list(self):
        """A list of file names in :py:attr:`~swamp.utils.swamplibrary.SwampLibrary.pdb_library`"""
        return [os.path.join(self.pdb_library, fname) for fname in os.listdir(self.pdb_library)]

    @property
    def _pdbfname_template(self):
        """A template file name for pdb files"""
        return os.path.join(self.pdb_library, "{}", "pdb{}.ent.gz")

    @property
    def _xmlfname_template(self):
        """A template file name for xml files"""
        return os.path.join(self.pdbtm_svn, "{}", "{}.xml")

    @property
    def _fragfile_template(self):
        """A template file name for fragment pdb files"""
        return os.path.join(self.pdb_library, '{}.pdb')

    @property
    def _library_out_template(self):
        """A template file name for pdb files with helical pairs"""
        return os.path.join(self.workdir, '{}_{}{}_{}{}.{}')

    @property
    def _ensemble_pdbout_template(self):
        """A template pdb file name for ensemble output files"""
        if self.outdir is None:
            return None
        else:
            return os.path.join(self.outdir, 'ensemble_{}.pdb')

    @property
    def _centroid_template(self):
        """A centroid pdb file name template"""
        if self.outdir is None:
            return None
        else:
            return os.path.join(self.outdir, 'centroid_{}.pdb')

    # ------------------ Hidden methods ------------------

    def _is_valid_entry(self, pdbcode):
        """For a given pdb code, check if there is a PDB and a XML file in the \
        :py:attr:`~swamp.utils.swamplibrary.pdb_library` and :py:attr:`~swamp.utils.swamplibrary.pdbtm_svn` respectively

        :param str pdbcode: the pdb code of interest
        :returns: True if all the files are present (bool)
        """

        if not os.path.isfile(self._pdbfname_template.format(pdbcode[1:3], pdbcode)):
            self.logger.warning("Entry %s not found in input dir %s" % (pdbcode, self.pdb_library))
            return False
        if not os.path.isfile(self._xmlfname_template.format(pdbcode[1:3], pdbcode)):
            self.logger.warning("Entry %s not found in input dir %s" % (pdbcode, self.pdbtm_svn))
            return False
        return True

    def _make_workdir(self):
        """Create the :py:attr:`~swamp.utils.swamplibrary.workdir`

        :raises ValueError: if :py:attr:`~swamp.utils.swamplibrary.workdir` is None
        """

        if self.workdir is None:
            raise ValueError("Impossible to create workdir, please set workdir value first!")

        if not os.path.isdir(self.workdir):
            os.mkdir(self.workdir)

    def _determine_orientation(self, frag_ids):
        """For a given set of fragment ids, determine the optimal orientation to ensemble them and return the tuple of \
        file names

        :param frag_ids: a list with the fragment ids of interest
        :type frag_ids: list, tuple
        :returns: a tuple with the file names of the alignment that scored the highest qscore (tuple)
        :raises ValueError: if there are less than 2 fragments in the input list
        """

        if len(frag_ids) < 2:
            raise ValueError("Impossible to determine the orientation of less than two fragments!")

        qscores = []
        frag_list = [(frag, SwampLibrary._get_reciprocal_id(frag)) for frag in frag_ids]

        all_combinations = list(itertools.product(*frag_list))

        for combination in all_combinations:
            gesamt = Gesamt(pdbin=[self._fragfile_template.format(x) for x in combination], workdir=None,
                            mode="alignment")
            gesamt.run()
            qscores.append(float(gesamt.summary_results["qscore"]))

        return all_combinations[qscores.index(max(qscores))]

    # ------------------ Public methods ------------------

    def remove_homologs(self, pdb_ids_to_remove):
        """Remove fragments originating from a set of pdb structures out of \
        :py:attr:`~swamp.utils.swamplibrary.SwampLibrary.qscore_matrix`, \
        :py:attr:`~swamp.utils.swamplibrary.SwampLibrary.rmsd_matrix`, \
        :py:attr:`~swamp.utils.swamplibrary.SwampLibrary.nalign_matrix`

        :argument tuple pdb_ids_to_remove: tuple with the pdb codes of the structures to be removed
        """

        # Detect the fragments comming from homolog structures (convert everything to lower case)
        pdb_ids_to_remove = [pdb.lower() for pdb in pdb_ids_to_remove]
        frag_ids_to_remove = []
        for frag_id in self.qscore_matrix.columns:
            if frag_id.split('_')[0].lower() in pdb_ids_to_remove:
                frag_ids_to_remove.append(frag_id)
        # Remove the fragments
        self.qscore_matrix.drop(frag_ids_to_remove, 0, inplace=True)
        self.qscore_matrix.drop(frag_ids_to_remove, 1, inplace=True)
        self.rmsd_matrix.drop(frag_ids_to_remove, 0, inplace=True)
        self.rmsd_matrix.drop(frag_ids_to_remove, 1, inplace=True)
        self.nalign_matrix.drop(frag_ids_to_remove, 0, inplace=True)
        self.nalign_matrix.drop(frag_ids_to_remove, 1, inplace=True)

    def create_distance_mtx(self, gesamt_dir):
        """Create the square distance matrices for the library: \
        :py:attr:`~swamp.utils.swamplibrary.SwampLibrary.qscore_matrix`, \
        :py:attr:`~swamp.utils.swamplibrary.SwampLibrary.rmsd_matrix` and \
        :py:attr:`~swamp.utils.swamplibrary.SwampLibrary.nalign_matrix`

        Requires the :py:func:`~swamp.utils.swamplibrary.SwampLibrary.all_vs_all_gesamt` results. The distance \
        matrices contain the optimal structural alignment between every set of fragments present in the library.

        :param str gesamt_dir: directory containing the .hit files resulting from the all vs all gesamt search results
        """

        frag_dict = self._get_frag_id_dict(gesamt_dir)
        self.qscore_matrix = pd.DataFrame()
        self.qscore_matrix["frag_id"] = list(frag_dict.keys())
        self.rmsd_matrix = pd.DataFrame()
        self.rmsd_matrix["frag_id"] = list(frag_dict.keys())
        self.nalign_matrix = pd.DataFrame()
        self.nalign_matrix["frag_id"] = list(frag_dict.keys())
        self.logger.info("Creating distance matrices now...")
        n_frags = len(frag_dict.keys())

        for idx, unique_id in enumerate(frag_dict.keys()):
            self.logger.info("Working on entry %s (%s/%s)" % (unique_id, idx + 1, n_frags))
            fragment_distances = None

            for hitfile in frag_dict[unique_id]:
                # Get the current distances
                current_hits = Gesamt.parse_hitfile(hitfile)
                current_hits.drop("n_res", 1, inplace=True)
                current_hits.drop("seq_id", 1, inplace=True)
                current_hits.drop("rmsd", 1, inplace=True)
                current_hits.fname = current_hits.fname.str.replace('.pdb', '')
                current_hits.rename(columns={"fname": "frag_id"}, inplace=True)
                current_hits.frag_id = current_hits.frag_id.apply(lambda x: self._get_unique_frag_id(x))
                current_hits["max_qscore"] = current_hits.groupby(["frag_id"], sort=False)["qscore"].transform(max)
                current_hits = current_hits[current_hits.qscore == current_hits.max_qscore]
                current_hits.drop("qscore", 1, inplace=True)
                current_hits.drop_duplicates(inplace=True)
                # Append results to the current fragment distances
                fragment_distances = pd.concat([fragment_distances, current_hits]).reset_index(drop=True)

            # Get the final distances for this fragment
            fragment_distances.rename(columns={'max_qscore': 'qscore'}, inplace=True)
            fragment_distances['max_qscore'] = fragment_distances.groupby(["frag_id"], sort=False)["qscore"].transform(
                max)
            fragment_distances = fragment_distances[fragment_distances.qscore == fragment_distances.max_qscore]
            fragment_distances.drop("max_qscore", 1, inplace=True)
            fragment_distances.drop_duplicates(subset='frag_id', inplace=True)
            # Store it in the final matrix
            self.qscore_matrix = self.qscore_matrix.merge(fragment_distances.loc[:, ['frag_id', 'qscore']], how="left",
                                                          on=["frag_id"])
            self.qscore_matrix.rename(columns={'qscore': unique_id}, inplace=True)
            self.nalign_matrix = self.nalign_matrix.merge(fragment_distances.loc[:, ['frag_id', 'n_align']], how="left",
                                                          on=["frag_id"])
            self.nalign_matrix.rename(columns={'n_align': unique_id}, inplace=True)
            self.rmsd_matrix = self.rmsd_matrix.merge(fragment_distances.loc[:, ['frag_id', 'rmsd']], how="left",
                                                      on=["frag_id"])
            self.rmsd_matrix.rename(columns={'rmsd': unique_id}, inplace=True)

        self.rmsd_matrix = self.rename_axis(self.rmsd_matrix)
        self.nalign_matrix = self.rename_axis(self.nalign_matrix)
        self.qscore_matrix = self.rename_axis(self.qscore_matrix)

    def make_library(self, pdb_codes):
        """Create the pdb files for each contacting TM helical pair in detected with the information at \
        :py:attr:`~swamp.utils.swamplibrary.pdb_library` and  :py:attr:`~swamp.utils.swamplibrary.pdbtm_svn`. Files \
        will be created at :py:attr:`~swamp.utils.swamplibrary.workdir`

        :param list pdb_codes: a list with the pdb codes that will be included to the library
        """

        for idx, entry in enumerate(pdb_codes):

            pdbcode = entry[0]
            chain = entry[1]
            if not self._is_valid_entry(pdbcode):
                self.logger.warning("Skipping invalid entry %s" % pdbcode)
                continue
            self.logger.info(
                "Processing %s:%s entry to the library (%s / %s)" % (pdbcode, chain, idx + 1, len(pdb_codes)))

            # TM helices
            pdbtm_parser = PdbtmXmlParser(self._xmlfname_template.format(pdbcode[1:3], pdbcode))
            pdbtm_parser.parse()
            tmhelices = [ss_annot for ss_annot in pdbtm_parser.ss2_annotation if
                         ss_annot.type == "H" and ss_annot.chain == chain]

            # Extract pdb hierarchy
            full_hierarchy = gemmi.read_structure(self._pdbfname_template.format(pdbcode[1:3], pdbcode))
            if full_hierarchy.info.__getitem__('_exptl.method') != 'X-RAY DIFFRACTION':
                self.logger.info('Not a X-ray structure, skipping...')
                continue
            full_hierarchy.remove_waters()

            # Check helical pairs individually
            for idx, helix_a in enumerate(tmhelices):
                for helix_b in tmhelices[idx + 1:]:

                    helix_a_hierarchy = utils.extract_hierarchy(to_extract=helix_a.pdb_region,
                                                                chainID=helix_a.chain,
                                                                full_hierarchy=full_hierarchy)
                    helix_b_hierarchy = utils.extract_hierarchy(to_extract=helix_b.pdb_region,
                                                                chainID=helix_b.chain,
                                                                full_hierarchy=full_hierarchy)
                    fragment_hierarchy = utils.merge_hierarchies((helix_a_hierarchy, helix_b_hierarchy),
                                                                 renumber=False)
                    fragment_cmap = utils.extract_fragment_cmap(fragment_hierarchy,
                                                                (helix_a.pdb_region, helix_b.pdb_region))

                    if fragment_cmap is None:
                        self.logger.warning(
                            "No contacts loaded from %s:%s %s - %s" % (pdbcode, chain, helix_a.index, helix_b.index))
                        continue

                    if len(fragment_cmap) >= 2:
                        self.logger.info(
                            "Found contacting helical pair! %s %s %s" % (pdbcode, helix_a.index, helix_b.index))

                        # Write pdb files
                        fragment_hierarchy.cell = full_hierarchy.cell
                        utils.renumber_hierarchy(fragment_hierarchy)
                        inverted_fragment = utils.invert_hiearchy(fragment_hierarchy)
                        inverted_fragment.cell = full_hierarchy.cell
                        pdbout = self._library_out_template.format(pdbcode, helix_a.index, helix_a.chain, helix_b.index,
                                                                   helix_b.chain, "pdb")
                        fragment_hierarchy.write_pdb(pdbout)
                        pdbout = self._library_out_template.format(pdbcode, helix_b.index, helix_b.chain, helix_a.index,
                                                                   helix_a.chain, "pdb")
                        inverted_fragment.write_pdb(pdbout)

                        # Write contact maps
                        conkit.io.write(
                            self._library_out_template.format(pdbcode, helix_a.index, helix_a.chain, helix_b.index,
                                                              helix_b.chain, "mapalign"), "mapalign", fragment_cmap)
                        conkit.io.write(
                            self._library_out_template.format(pdbcode, helix_a.index, helix_a.chain, helix_b.index,
                                                              helix_b.chain, "aleigen"), "aleigen", fragment_cmap)
                        inverted_cmap = utils.invert_contactmap(fragment_cmap)
                        conkit.io.write(
                            self._library_out_template.format(pdbcode, helix_b.index, helix_b.chain, helix_a.index,
                                                              helix_a.chain, "mapalign"), "mapalign", inverted_cmap)
                        conkit.io.write(
                            self._library_out_template.format(pdbcode, helix_b.index, helix_b.chain, helix_a.index,
                                                              helix_a.chain, "aleigen"), "aleigen", inverted_cmap)

    def all_vs_all_gesamt(self, inputdir, outdir, nthreads=1):
        """For each the members of the library, obtain the distance with all the others. This step is required to \
        obtain the distance matrices: :py:attr:`~swamp.utils.swamplibrary.SwampLibrary.qscore_matrix`, \
        :py:attr:`~swamp.utils.swamplibrary.SwampLibrary.rmsd_matrix` and \
        :py:attr:`~swamp.utils.swamplibrary.SwampLibrary.nalign_matrix`

        :param str inputdir: the input directory with the pdb files created by \
        :py:func:`~swamp.utils.swamplibrary.SwampLibrary.make_library`
        :param str outdir: the output directory where the .hit files will be created
        :param int nthreads: number of threads to be used in the gesamt archive scan (default 1)
        """

        # Make the archive
        self.logger.info("Creating gesamt archive at %s" % os.path.join(outdir, "gesamt_archive"))
        gesamt_makearchive = Gesamt(workdir=None, mode="make-archive", pdb_archive=inputdir, pdbin=None,
                                    gesamt_archive=os.path.join(outdir, "gesamt_archive"))
        gesamt_makearchive.run()

        # Scan the archive with all the fragments in the library
        self.logger.info("Now scanning archive with all fragments in the library...")
        fragment_list = [os.path.join(inputdir, fname) for fname in os.listdir(inputdir) if
                         fname[-4:] == ".pdb" and os.path.isfile(os.path.join(inputdir, fname))]
        for pdbfile in fragment_list:
            id = os.path.basename(pdbfile)[:-4]
            gesamt = Gesamt(mode="search-archive", pdbin=pdbfile, gesamt_archive=os.path.join(outdir, "gesamt_archive"),
                            min2="0", nthreads=str(nthreads), min1="0", workdir=None,
                            hits_out=os.path.join(outdir, "%s_hits.txt" % id))
            gesamt.run()

    # ------------------ Static methods ------------------

    @staticmethod
    def rename_axis(df):
        """Rename the axis of a `pandas.DataFrame` so that the row names correspond with the column names

        :param `pandas.DataFrame` df: the dataframe to be renamed
        :returns: renamed dataframe (`pandas.DataFrame`)
        """

        new_rownames = {}
        for rowidx in [x for x in list(df.index)]:
            new_rownames[rowidx] = df.columns[rowidx + 1]

        df = df.rename(index=new_rownames, inplace=False)
        return df.drop("frag_id", 1)

    @staticmethod
    def _get_unique_frag_id(frag_id):
        """For a given fragment id (pdb_tmhelix_tmhelix) obtain the unique fragment id (tmhelix fields are sorted)

         The unique id corresponds with the pdb code of the structure where the fragment was found, plus the sorted
         indeces of the two TM helices that form the fragment. This fragment id serves as an unique pointer for each
         component of the library.

         :param str frag_id: the fragment id of interest
         :returns: unique fragment id (str)
         :raises ValueError: if the fragment id has more than 3 components
         """

        frag_id = frag_id.split("_")
        if len(frag_id) > 3:
            raise ValueError("The frag ID cannot have more than 3 components! %s" % "_".join(frag_id))
        frag_id[1:] = sorted(frag_id[1:])
        return "_".join(frag_id)

    @staticmethod
    def _get_reciprocal_id(frag_id):
        """Method to get the reciprocal fragment id.

        A reciprocal fragment corresponds with the same fragment id, but the order at which the TM helices appear
        at the sequence the fragment has been inverted.

        :param str frag_id: fragment id of interest
        :returns: reciprocal id where the helical order is inverted (str)
        """

        frag_id = frag_id.split("_")
        if len(frag_id) > 3:
            raise ValueError("The frag ID cannot have more than 3 components! %s" % "_".join(frag_id))
        frag_id[1:] = reversed(frag_id[1:])
        return "_".join(frag_id)

    @staticmethod
    def _get_frag_id_dict(gesamt_dir):
        """In a directory full of gesamt .hit files, parse all the fragment ids and their hit files into a dictionary

        It will compute both orientations into a dictionary where the key is the unique frag id

        :param str gesamt_dir: directory where the .hit files are located
        :returns: a dictionary with the unique fragment id as key and a list with the hit filenames as values
        """

        result = {}

        for hitfile in os.listdir(gesamt_dir):
            hitfile = os.path.join(gesamt_dir, hitfile)
            if not os.path.isfile(hitfile) or not os.path.basename(hitfile)[-4:] == ".txt":
                continue
            frag_id = os.path.basename(hitfile)[:-9]
            unique_id = SwampLibrary._get_unique_frag_id(frag_id)
            if unique_id in result.keys():
                result[unique_id].append(hitfile)
            else:
                result[unique_id] = [hitfile]

        return result

    @staticmethod
    def parse_nr_listfile(fname):
        """Method to parse a file with pdb structures listed as PDB:CHAIN into a nested list

        :param str fname: file name with the list to be parsed
        :returns: a nested list where each element contains the non-redundant pdb code and chain name (list)
        """

        nr_pdbIDsChains = []
        with open(fname, "r") as nr_list_file:
            for line in nr_list_file:
                nr_pdbIDsChains.append((line[:4], line.rstrip()[-1]))
        return nr_pdbIDsChains
