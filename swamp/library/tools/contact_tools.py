import os
import conkit.io
from swamp.tools.pdb_tools import get_tempfile
from conkit.core import Contact, ContactMap, Sequence


def invert_contactmap(cmap):
    """Method to invert a contact map

    :param cmap: the contact map of interest
    :type cmap: :object:`conkit.core.ContactMap`
    :returns inverted_cmap: the contact map corresponding with the inverted sequence (1-res_seq)
    :rtype :object:`conkit.core.ContactMap`
    """

    inverted_cmap = ContactMap('inverted')

    for contact in cmap:
        new_contact = Contact(cmap.highest_residue_number + 1 - contact.res1_seq,
                              cmap.highest_residue_number + 1 - contact.res2_seq,
                              contact.raw_score)
        inverted_cmap.add(new_contact)

    inverted_cmap.sequence = cmap.sequence
    return inverted_cmap


def extract_interhelical_cmap(cmap, helices, residues, new_id, score_threshold=0.0):
    """Method to extract the interhelical contacts

    :param cmap: the contact map of interest
    :type cmap: :object:`conkit.core.ContactMap`
    :param helices: a nested list with the listed residue numbers of the two helices of interest
    :type helices: tuple, list
    :param new_id: the new identifies given to the resulting contact map
    :type new_id: str, int
    :param score_threshold: the raw score threshold at which contacts will be included (default 0.0)
    :type score_threshold: float
    :returns inverted_cmap: the contact map containing only the inter-helical contacts between the pair of helices
    :rtype :object:`conkit.core.ContactMap`
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

    :param pdb_hierarchy: the pdb hierarchy of the fragment
    :type pdb_hierarchy: :obj:`gemmi.Structure`
    :param helices: a nested list with the listed residues numbers of the helices of interest
    :type helices: tuple, list
    :returns the contact map with the interhelical contacts
    :rtype :object:`conkit.core.ContactMap`
    """

    residues = [i for sub in helices for i in sub]
    temp_pdb = get_tempfile()
    pdb_hierarchy.write_pdb(temp_pdb)
    fragment_cmap = conkit.io.read(temp_pdb, "pdb").top_map
    os.remove(temp_pdb)

    if fragment_cmap is None:
        return None

    return extract_interhelical_cmap(fragment_cmap, helices, residues, "InterhelicalContacts")
