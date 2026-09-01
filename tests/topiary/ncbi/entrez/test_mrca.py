import pytest

from topiary.ncbi.entrez.mrca import get_mrca_taxid


def _record(lineage_taxids, taxid):
    """
    Build the record shape Entrez.read returns for a taxonomy efetch: a
    LineageEx list of ancestors, plus the taxon's own TaxId.
    """

    return {"LineageEx": [{"TaxId": str(t)} for t in lineage_taxids],
            "TaxId": str(taxid)}


def test_get_mrca_taxid_no_species_is_root_of_life(mocker):

    mocker.patch("topiary.ncbi.get_taxid", return_value=[])

    # Nothing to find a common ancestor of, so fall all the way back to root
    assert get_mrca_taxid([]) == 1


def test_get_mrca_taxid_single_species_is_itself(mocker):

    mocker.patch("topiary.ncbi.get_taxid", return_value=["9606"])

    # The MRCA of one taxon is that taxon. No efetch needed.
    efetch = mocker.patch("topiary.ncbi.entrez.mrca.Entrez.efetch")

    assert get_mrca_taxid(["Homo sapiens"]) == 9606
    efetch.assert_not_called()


def test_get_mrca_taxid_finds_shared_ancestor(mocker):
    """
    Human and mouse share the lineage up through Mammalia, then diverge.
    """

    mocker.patch("topiary.ncbi.get_taxid", return_value=["9606", "10090"])
    mocker.patch("topiary.ncbi.entrez.mrca.Entrez.efetch")

    # root(1) -> Chordata(7711) -> Mammalia(40674) -> then diverge
    human = _record([7711, 40674, 9443], 9606)
    mouse = _record([7711, 40674, 9989], 10090)

    mocker.patch("topiary.ncbi.entrez.mrca.Entrez.read",
                 return_value=[human, mouse])

    # Lineages are [1,7711,40674,9443,9606] and [1,7711,40674,9989,10090];
    # they agree through index 2, so the MRCA is Mammalia.
    assert get_mrca_taxid(["Homo sapiens", "Mus musculus"]) == 40674


def test_get_mrca_taxid_unrelated_species_fall_back_to_root(mocker):
    """
    Two lineages that diverge immediately below root have root as their MRCA.
    """

    mocker.patch("topiary.ncbi.get_taxid", return_value=["9606", "3702"])
    mocker.patch("topiary.ncbi.entrez.mrca.Entrez.efetch")

    animal = _record([33208], 9606)
    plant = _record([33090], 3702)

    mocker.patch("topiary.ncbi.entrez.mrca.Entrez.read",
                 return_value=[animal, plant])

    assert get_mrca_taxid(["Homo sapiens", "Arabidopsis thaliana"]) == 1


def test_get_mrca_taxid_identical_lineages(mocker):
    """
    Two records for the same taxon: the MRCA is that taxon, and the loop must
    run to the end of the shortest lineage without breaking.
    """

    mocker.patch("topiary.ncbi.get_taxid", return_value=["9606", "9606"])
    mocker.patch("topiary.ncbi.entrez.mrca.Entrez.efetch")

    rec = _record([7711, 40674], 9606)
    mocker.patch("topiary.ncbi.entrez.mrca.Entrez.read",
                 return_value=[rec, rec])

    assert get_mrca_taxid(["Homo sapiens", "Homo sapiens"]) == 9606


def test_get_mrca_taxid_nested_lineages(mocker):
    """
    One lineage is a prefix of the other (a genus and a species within it).
    The loop is bounded by the shorter lineage, so the MRCA is the ancestor.
    """

    mocker.patch("topiary.ncbi.get_taxid", return_value=["9605", "9606"])
    mocker.patch("topiary.ncbi.entrez.mrca.Entrez.efetch")

    genus = _record([7711, 40674], 9605)
    species = _record([7711, 40674, 9605], 9606)

    mocker.patch("topiary.ncbi.entrez.mrca.Entrez.read",
                 return_value=[genus, species])

    # [1,7711,40674,9605] vs [1,7711,40674,9605,9606] -- agree for all four
    # positions of the shorter one
    assert get_mrca_taxid(["Homo", "Homo sapiens"]) == 9605


def test_get_mrca_taxid_empty_records_is_root(mocker):
    """
    efetch came back with nothing to work from.
    """

    mocker.patch("topiary.ncbi.get_taxid", return_value=["9606", "10090"])
    mocker.patch("topiary.ncbi.entrez.mrca.Entrez.efetch")
    mocker.patch("topiary.ncbi.entrez.mrca.Entrez.read", return_value=[])

    assert get_mrca_taxid(["Homo sapiens", "Mus musculus"]) == 1


def test_get_mrca_taxid_queries_every_taxid(mocker):
    """
    All taxids must reach efetch in one comma-joined request -- dropping one
    would silently widen the MRCA.
    """

    mocker.patch("topiary.ncbi.get_taxid", return_value=["9606", "10090", "7955"])
    efetch = mocker.patch("topiary.ncbi.entrez.mrca.Entrez.efetch")
    mocker.patch("topiary.ncbi.entrez.mrca.Entrez.read",
                 return_value=[_record([7711], 9606),
                               _record([7711], 10090),
                               _record([7711], 7955)])

    get_mrca_taxid(["Homo sapiens", "Mus musculus", "Danio rerio"])

    efetch.assert_called_once()
    assert efetch.call_args.kwargs["id"] == "9606,10090,7955"
    assert efetch.call_args.kwargs["db"] == "taxonomy"
