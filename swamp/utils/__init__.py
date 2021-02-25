"""This is SWAMP: Solving structures With Alpha Membrane Pairs

This module implements useful classes and methods used across all modules in SWAMP
"""

__author__ = "Filomeno Sanchez Rodriguez"
__credits__ = "Daniel Rigden, & Ronan Keegan"
__email__ = "filomeno.sanchez-rodriguez@liv.ac.uk"

import sys
import os
import gzip
import shutil
import tempfile
import logging
from swamp import version

__version__ = version.__version__

if 'DISABLE_DEPENDENCY_CHECKS' not in os.environ:

    if "CCP4" not in os.environ:
        raise RuntimeError("Cannot find CCP4 root directory")

    import gemmi
    import conkit.io
    from conkit.core import Contact, ContactMap, Sequence
    from Bio.PDB.parse_pdb_header import _parse_pdb_header_list


def SwampLibrary(*args, **kwargs):
    """:py:obj:`~swamp.utils.swamplibrary.SwampLibrary` instance"""
    from swamp.utils.swamplibrary import SwampLibrary

    return SwampLibrary(*args, **kwargs)


def ThreadResults(*args, **kwargs):
    """:py:obj:`~swamp.utils.threadresults.ThreadResults` instance"""
    from swamp.utils.threadresults import ThreadResults

    return ThreadResults(*args, **kwargs)


def TargetSplit(*args, **kwargs):
    """:py:obj:`~swamp.utils.targetsplit.TargetSplit` instance"""
    from swamp.utils.targetsplit import TargetSplit

    return TargetSplit(*args, **kwargs)


def compress(fname, out=None):
    """Compress a text file into .gz

    :param str fname: the file name to be compressed
    :param str out: specify an output file name, otherwise default is fname.gz
    :returns: compressed file name (str)
    """

    if out is None:
        out = '%s.gz' % fname

    with open(fname, 'rb') as f_in, gzip.open(out, 'wb') as f_out:
        data = f_in.read()
        if sys.version_info[0] < 3:
            bindata = data
        else:
            bindata = bytearray(data)
        f_out.write(bindata)

    return out


def decompress(fname, out=None):
    """Decompress a .gz file into text file

    :param str fname: the file name to be decompressed
    :param str out: specify an output file name, otherwise default is fname without .gz
    :returns: the decompressed file name (str)
    """

    if out is None:
        out = fname.replace('.gz', '')

    with open(out, "wb") as f_out, gzip.open(fname, "rb") as f_in:
        bindata = f_in.read()
        f_out.write(bindata)

    return out


def touch(fname, content='', mode='w'):
    """Create a file with the specified contents

    :param str fname: file name to be created
    :param str content: content to write into the file (default '')
    :param str mode: mode to open the file handler (default: 'w')
    """

    with open(fname, mode) as fhandle:
        fhandle.write(content)
    fhandle.close()


def get_tempfile():
    """Method to get a temporary file name

    :returns: temporary file name (str)
    """

    temp_name = next(tempfile._get_candidate_names())
    return os.path.join(os.environ['CCP4_SCR'], '%s.pdb' % temp_name)


def remove(path):
    if os.path.exists(path):
        if os.path.isdir(path):
            shutil.rmtree(path)
        else:
            os.remove(path)


def create_tempfile(content, mode="w"):
    """Create a temporary file with a given set of contents

    :param str content: content to dump into the temporary file
    :param str mode: mode to open the file handler (default: 'w')
    :returns: the path to the temporary file name (str)
    """

    fname = get_tempfile()
    touch(fname, content, mode)
    return fname


def invert_contactmap(cmap):
    """Method to invert a contact map

    :param :py:obj:`~conkit.core.ContactMap` cmap: the contact map of interest
    :returns: and inverted_cmap: the contact map corresponding with the inverted sequence (1-res_seq) \
    (:py:obj:`~conkit.core.ContactMap`)
    """

    inverted_cmap = ContactMap('inverted')

    highest_residue_number = max([max(contact.id) for contact in cmap])
    for contact in cmap:
        new_contact = Contact(highest_residue_number + 1 - contact.res1_seq,
                              highest_residue_number + 1 - contact.res2_seq,
                              contact.raw_score)
        inverted_cmap.add(new_contact)

    inverted_cmap.sequence = cmap.sequence
    return inverted_cmap


def extract_interhelical_cmap(cmap, helices, residues, new_id, score_threshold=0.0):
    """Method to extract the interhelical contacts

    :param :py:obj:`~conkit.core.ContactMap` cmap: the contact map of interest
    :param tuple helices: a nested list with the listed residue numbers of the two helices of interest
    :param str new_id: the new identifies given to the resulting contact map
    :param float score_threshold: the raw score threshold at which contacts will be included (default 0.0)
    :returns: inverted_cmap: the contact map containing only the inter-helical contacts between the pair of helices \
    (:py:obj:`~conkit.core.ContactMap`)
    """

    result_cmap = ContactMap(new_id)
    dummy_sequence = Sequence('id', 'A' * len(residues))

    for contact in cmap:
        if contact.raw_score >= score_threshold:
            helix_1 = [(contact.res1_seq in x) for x in helices]
            helix_2 = [(contact.res2_seq in x) for x in helices]
            if helix_1 != helix_2 and any(helix_1) and any(helix_2):
                # Create a new contact (IMPORTANT: Renumber the contact position within the new map)
                new_contact = Contact(residues.index(contact.res1_seq) + 1,
                                      residues.index(contact.res2_seq) + 1, contact.raw_score)
                result_cmap.add(new_contact)

    result_cmap.sequence = dummy_sequence
    return result_cmap


def extract_fragment_cmap(pdb_hierarchy, helices):
    """Method to extract the interhelical contact map of a given pdb file

    :param :py:obj:`~gemmi.Structure` pdb_hierarchy: the pdb hierarchy of the fragment
    :param tuple helices: a nested list with the listed residues numbers of the helices of interest
    :returns: the contact map with the interhelical contacts (:py:obj:`~conkit.core.ContactMap`)
    """

    residues = [i for sub in helices for i in sub]
    temp_pdb = get_tempfile()
    pdb_hierarchy.write_pdb(temp_pdb)
    fragment_cmap = conkit.io.read(temp_pdb, "pdb").top_map
    os.remove(temp_pdb)

    if fragment_cmap is None:
        return None

    return extract_interhelical_cmap(fragment_cmap, helices, residues, "InterhelicalContacts")


def merge_hierarchies(hiearchies, new_chain_id="A", new_model_id="1", renumber=False):
    """Method to merge two given hierarchies into one (same chain and model)

    :param tuple hiearchies: a list with the pdb hierarchies to be merged
    :param str new_chain_id: the new chain id for the result hierarchy
    :param str new_model_id: the new model name for the result hierarchy
    :param bool renumber: if True the residues of the resulting hierarchy will be renumbered starting at 1
    :returns: a new :py:obj:`~gemmi.Structure` hierarchy corresponding to the merged input hierarchies
    """

    if not isinstance(hiearchies, list) and not isinstance(hiearchies, tuple):
        raise ValueError("Please provide hierarchies to be merged as lists!")
    if len(hiearchies) < 2:
        raise ValueError("Please provide at least two hierarchies to merge!")

    new_model = gemmi.Model(new_model_id)
    new_chain = gemmi.Chain(new_chain_id)
    new_hierarchy = gemmi.Structure()

    for hierarchy in hiearchies:
        for res in hierarchy[0][0]:
            new_chain.add_residue(res)

    new_model.add_chain(new_chain)
    new_hierarchy.add_model(new_model)
    if renumber:
        renumber_hierarchy(new_hierarchy)

    return new_hierarchy


def renumber_hierarchy(hierarchy, start=1):
    """Method to renumber a given hierarchy to start in a given value. Renumbered inplace

    :param :py:obj:`~gemmi.Structure` hierarchy: pdb hierarchy to be renumbered
    :param int start: first residue to start renumbering of the hierarchy
    """

    atom_idx = 1
    for model in hierarchy:
        for chain in model:
            for idx, residue in enumerate(chain):
                residue.seqid.num = idx + start
                for atom in residue:
                    atom.serial = atom_idx
                    atom_idx += 1


def extract_hierarchy(full_hierarchy, to_extract, chainID=None):
    """Method to extract a given set of residues from a pdbfile in form of a gemmi structure hierarchy

    :param :py:obj:`~gemmi.Structure` full_hierarchy: pdb hierarchy of interest
    :param tuple to_extract: list with the residue numbers to be extracted
    :param str chainID: the chain id where the residues to be extracted are located (default None)
    :returns: a new :py:obj:`~gemmi.Structure` hierarchy containing the extracted residues
    """

    new_hierarchy = gemmi.Structure()
    new_chain = gemmi.Chain("A")
    new_model = gemmi.Model("1")

    # Check model number
    if len(full_hierarchy) > 1:
        logging.debug("pdb {0} has > 1r model - only first model will be kept".format(full_hierarchy.name))

    # Check chain number
    if len(full_hierarchy[0]) > 1:
        if chainID is None:
            logging.debug(
                "pdb {0} has > 1 chain - only first chain will be kept".format(full_hierarchy.name))
            old_chain = full_hierarchy[0][0]
        else:
            try:
                old_chain = full_hierarchy[0][chainID]
            except ValueError as e:
                raise ValueError("Chain %s not found in %s!" % (chainID, full_hierarchy.name))
    else:
        old_chain = full_hierarchy[0][0]

    # Extract region
    for residue in old_chain:
        if residue.seqid.num in to_extract:
            new_chain.add_residue(residue)

    # Append the model to the new hierarchy and exit
    new_model.add_chain(new_chain)
    new_hierarchy.add_model(new_model)
    return new_hierarchy


def invert_hiearchy(hierarchy):
    """Method to return the inverted hierarchy (1-res_seq)

    :param :py:obj:`~gemmi.Structure` hierarchy: pdb hierarchy to be inverted
    :returns: the :py:obj:`~gemmi.Structure` hierarchy corresponding with the inverted sequence (1-res_seq)
    """

    inverted_model = gemmi.Model("1")
    inverted_chain = gemmi.Chain("A")
    inverted_hierarchy = gemmi.Structure()

    tmp_list = []
    for residue in hierarchy[0][0]:
        tmp_list.append(residue)

    for idx, residue in enumerate(tmp_list[::-1]):
        inverted_chain.add_residue(residue)
        inverted_chain[-1].seqid.num = idx + 1

    inverted_model.add_chain(inverted_chain)
    inverted_hierarchy.add_model(inverted_model)
    renumber_hierarchy(inverted_hierarchy)

    return inverted_hierarchy


def get_missing_residues(header_list):
    """Get a dictionary with the missing residues described in the REMARK section of a pdb file

    :param tuple header_list: a list with the lines of the header section of the pdb file
    :returns: a dictionary with the missing residues present in each chain (chain ids are used as keys)
    """

    head = _parse_pdb_header_list(header_list)
    rslt = {}

    for residue in head['missing_residues']:
        if residue['chain'] in rslt.keys():
            rslt[residue['chain']].append(residue['ssseq'])
        else:
            rslt[residue['chain']] = [residue['ssseq']]

    return rslt


def extract_hierarchy_seqnumber(hierarchy, seq_numbers, chain_id='A'):
    """ Extract the hierarchy corresponding with a given set of residue sequence numbers. Considers missing residues.

    :argument :py:obj:`~gemmi.Structure` hierarchy: original hierarchy to trim
    :argument tuple seq_numbers: residue sequence number to extract
    :argument str chain_id: chain where the residues should be extracted from
    :returns: a new :py:obj:`~gemmi.Structure` with the residues of interest

    :example

        >>> import gemmi
        >>> from swamp.utils import extract_hierarchy_seqnumber
        >>> hierarchy = gemmi.read_structure('/mnt/sda1/MR_edge_cases/3txt_MR/3txt.pdb')
        >>> subtrgt_1 = extract_hierarchy_seqnumber(hierarchy, [x for x in range(104 ,125)] + [x for x in range(158, 179)])
        >>> subtrgt_1.write_minimal_pdb('/home/filo/test.pdb')
    """

    header = hierarchy.make_pdb_headers().split('\n')
    try:
        missing_res = get_missing_residues(header)[chain_id]
    except KeyError:
        missing_res = []
    new_hierarchy = gemmi.Structure()
    new_model = gemmi.Model("1")
    new_chain = gemmi.Chain(chain_id)

    idx = 1
    for residue in hierarchy[0][chain_id]:
        factor = sum([1 for x in missing_res if x < residue.seqid.num])
        residue.seqid.num = idx + factor
        idx += 1
        if residue.seqid.num in seq_numbers:
            new_chain.add_residue(residue)
    new_model.add_chain(new_chain)
    new_hierarchy.add_model(new_model)
    return new_hierarchy


def merge_into_ensemble(hierarchies):
    """Method to merge a series of hierarchies into an ensemble where each of the original hierarchies is
     represented as a model

     :argument tuple hierarchies: the hierarchies that will merged to form an ensemble
     :returns: a new :py:obj:`~gemmi.Structure` hierarchy containing the ensemble
     """

    new_hierarchy = gemmi.Structure()

    for hierarchy in hierarchies:
        new_model = gemmi.Model(str(len(new_hierarchy) + 1))
        new_model.add_chain(hierarchy[0][0])
        new_hierarchy.add_model(new_model)

    return new_hierarchy


def split_ensemble_into_models(hierarchy):
    """Method to split a ensemble into its constituent models

    :argument :py:obj:`~gemmi.Structure` hierarchy: the input ensemble to be splited
    :returns: a tuple containing :py:obj:`~gemmi.Structure`, each formed by a single model originating from the input \
    ensemble
    """

    result = []

    for model in hierarchy:
        new_hierarchy = gemmi.Structure()
        new_hierarchy.add_model(model)
        result.append(new_hierarchy)

    return tuple(result)
