import gemmi
import logging
import tempfile
import os
from Bio.PDB.parse_pdb_header import _parse_pdb_header_list


def get_tempfile():
    """Method to get a temporary file name

    :returns temporary file name
    :rtype str
    """

    temp_name = next(tempfile._get_candidate_names())
    return os.path.join(os.environ['CCP4_SCR'], '%s.pdb' % temp_name)


def merge_hierarchies(hiearchies, new_chain_id="A", new_model_id="1", renumber=False):
    """Method to merge two given hierarchies into one (same chain and model)

    :param hiearchies: a list with the pdb hierarchies to be merged
    :type hiearchies: list, tuple
    :param new_chain_id: the new chain id for the result hierarchy
    :type new_model_id: str
    :param new_model_id: the new model name for the result hierarchy
    :type new_model_id:str
    :param renumber: if True the residues of the resulting hierarchy will be renumbered starting at 1
    :type renumber: bool
    :returns new_hierarchy: a new pdb hierarchy corresponding to the merged input hierarchies
    :rtype :obj:`gemmi.Structure`
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

    :param hierarchy: pdb hierarchy to be renumbered
    :type hierarchy: :obj:`gemmi.Structure`
    :param start: first residue to start renumbering of the hierarchy
    :type start: int
    :returns nothing
    :rtype None
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

    :param full_hierarchy: pdb hierarchy of interest
    :type full_hierarchy: :obj:`gemmi.Structure`
    :param to_extract: list with the residue numbers to be extracted
    :type to_extract: list, tuple
    :param chainID: the chain id where the residues to be extracted are located (default None)
    :type chainID: None, str
    :returns new_hierarchy: the pdb hierarchy containing the extracted residues
    :rtype :obj:`gemmi.Structure`
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

    :param hierarchy: pdb hierarchy to be inverted
    :type hierarchy: :obj:`gemmi.Structure`
    :returns inverted_hierarchy: the pdb hierarchy corresponding with the inverted sequence (1-res_seq)
    :rtype :obj:`gemmi.Structure`
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

    :param header_list: a list with the lines of the header section of the pdb file
    :type header_list: list, tuple
    :returns rslt: a dictionary with the missing residues present in each chain (chain ids are used as keys)
    :rtype dict
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

    :argument hierarchy: original hierarchy to trim
    :type hierarchy: :obj:`gemmi.Structure`
    :argument seq_numbers: residue sequence number to extract
    :type seq_numbers: list
    :argument chain_id: chain where the residues should be extracted from
    :type chain_id: str
    :returns new_hierarchy: a new hierarchy with the residues of interest
    :rtype :object:`gemmi.Structure`

    :example

        >>> import gemmi
        >>> from swamp.tools.pdb_tools import extract_hierarchy_seqnumber
        >>> hierarchy = gemmi.read_structure('/mnt/sda1/MR_edge_cases/3txt_MR/3txt.pdb')
        >>> subtrgt_1 = extract_hierarchy_seqnumber(hierarchy, [x for x in range(104 ,125)] + [x for x in range(158, 179)])
        >>> subtrgt_1.write_minimal_pdb('/home/filo/test.pdb')
    """

    header = hierarchy.make_pdb_headers().split('\n')
    missing_res = get_missing_residues(header)[chain_id]
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

     :argument hierarchies
     :type hierarchies: tuple, list
     :returns new_hierarchy: a new hierarchy containing the ensemble
     :rtype :object:`gemmi.Structure`
     """

    new_hierarchy = gemmi.Structure()

    for hierarchy in hierarchies:
        new_model = gemmi.Model(len(new_hierarchy) + 1)
        new_model.add_chain(hierarchy[0][0])
        new_hierarchy.add_model(new_model)

    return new_hierarchy


def split_ensemble_into_models(hierarchy):
    """Method to split a ensemble into its constituent models

    :argument hierarchy: the input ensemble to be splited
    :type hierarchy :object:`gemmi.Structure`
    :return a tuple containing hierarchies, each formed by a single model originating from the input ensemble
    :rtype tuple
    """

    result = []

    for model in hierarchy:
        new_hierarchy = gemmi.Structure()
        new_hierarchy.add_model(model)
        result.append(new_hierarchy)

    return tuple(result)
