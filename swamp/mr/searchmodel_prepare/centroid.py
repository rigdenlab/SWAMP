import os
from swamp.mr.searchmodel_prepare.prepare import PrepareSearchModel
from swamp.mr.searchmodel_prepare.polyala import PolyALA


class Centroid(PrepareSearchModel):
    """Centroid wrapper to prepare search model

    Examples
    --------
    work in progress...

    """

    # ------------------ Class specific properties ------------------

    @property
    def cmd(self):

        """Property to store the command to be executed"""

        return None

    @property
    def modification(self):

        """Property to store the modification to be applied"""

        return "centroid"

    # ------------------ Class specific methods ------------------

    def prepare(self):

        """Method to prepare the search model using molrep"""

        self.make_workdir()
        os.chdir(self.workdir)
        PolyALA.truncate_polyALA(self.model_list[0], self.pdbout)
        self.check_output()


