import os
from pyjob import cexec
from simbad.util import mtz_util
from swamp.wrappers.wrapper import Wrapper


# Class to hold refmac run
class Crank2(Wrapper):

    def __init__(self, workdir, pdbin, mtzin, fastafile, solvent='0.45', nchains_asu=1, anomalous_scatterer="S",
                 mode="SAD", wavelength_type="peak", logger=None):

        super(Crank2, self).__init__(workdir=workdir, logger=logger)

        self._solution = "NO"
        self._rfactor = "NA"
        self._rfree = "NA"
        self._FOM = "NA"
        self._wavelength_type = wavelength_type
        self._mode = mode
        self._anomalous_scatterer = anomalous_scatterer
        self._pdbin = pdbin
        self._fastafile = fastafile
        self._solvent = solvent
        self._nchains_asu = nchains_asu
        self._mtzin = mtzin

    @property
    def wrapper_name(self):
        return "crank2"

    @property
    def cmd(self):
        return ["ccp4-python", "%s/share/ccp4i/crank2/crank2.py" % os.environ["CCP4"], '--xyzout', self.pdbout,
                '--hklout', self.mtzout, '--dirout', self.workdir, '--logout', self.logout, '--keyin', self.keyin,
                '--graceful-preliminary-stop', '--disable-rvapi']

    @property
    def keywords(self):
        mtz_head = mtz_util.GetLabels(self.mtzin)
        return [
            "fsigf plus dname=%s file=%s i=%s sigi=%s" % (
                self.wavelength_type, self.mtzin, mtz_head.iplus, mtz_head.sigiplus),
            "fsigf minus dname=%s file=%s i=%s sigi=%s" % (
                self.wavelength_type, self.mtzin, mtz_head.iminus, mtz_head.sigiminus), "target::%s" % self.mode,
            "model substr atomtype=%s d_name=%s" % (self.anomalous_scatterer, self.wavelength_type), "refatompick",
            'model unknown "file=%s" atomtype=%s' % (self.pdbin, self.anomalous_scatterer),
            "sequence monomers_asym=%s solvent_content=%s file=%s" % (self.nchains_asu, self.solvent, self.fastafile),
            'comb_phdmmb exclude free=%s "file=%s" mb buccaneer' % (mtz_head.free, self.mtzin),
            'ref target::MLHL exclude free=%s "file=%s"' % (mtz_head.free, self.mtzin)
        ]

    @property
    def wavelength_type(self):
        """wavelength_type"""
        return self._wavelength_type

    @wavelength_type.setter
    def wavelength_type(self, value):
        self._wavelength_type = value

    @property
    def fastafile(self):
        """fastafile"""
        return self._fastafile

    @fastafile.setter
    def fastafile(self, value):
        self._fastafile = value

    @property
    def solvent(self):
        """solvent"""
        return self._solvent

    @solvent.setter
    def solvent(self, value):
        self._solvent = value

    @property
    def nchains_asu(self):
        """nchains_asu"""
        return self._nchains_asu

    @nchains_asu.setter
    def nchains_asu(self, value):
        self._nchains_asu = value

    @property
    def mtzin(self):
        """mtzin"""
        return self._mtzin

    @mtzin.setter
    def mtzin(self, value):
        self._mtzin = value

    @property
    def mode(self):
        """mode"""
        return self._mode

    @mode.setter
    def mode(self, value):
        self._mode = value

    @property
    def anomalous_scatterer(self):
        """anomalous_scatterer"""
        return self._anomalous_scatterer

    @anomalous_scatterer.setter
    def anomalous_scatterer(self, value):
        self._anomalous_scatterer = value

    @property
    def pdbin(self):
        """pdbin"""
        return self._pdbin

    @pdbin.setter
    def pdbin(self, value):
        self._pdbin = value

    @property
    def rfactor(self):
        """Rfactor"""
        return self._rfactor

    @rfactor.setter
    def rfactor(self, value):
        self._rfactor = value

    @property
    def rfree(self):
        """Rfree"""
        return self._rfree

    @rfree.setter
    def rfree(self, value):
        self._rfree = value

    @property
    def solution(self):
        """solution"""
        return self._solution

    @solution.setter
    def solution(self, value):
        self._solution = value

    @property
    def FOM(self):
        """FOM"""
        return self._FOM

    @FOM.setter
    def FOM(self, value):
        self._FOM = value

    @property
    def keyin(self):
        return os.path.join(self.run_info.workdir, "crank2_keyin.txt")

    @property
    def logout(self):
        return os.path.join(self.run_info.workdir, "crank2_logout.txt")

    def _create_keyin(self):
        """Create the kwyword input file"""

        with open(self.run_info.keyin, "w") as keyin_file:
            keyin_file.write("\n".join(self.keywrds))

    def _check_valid_mtz(self):
        """ check the sanity of the mtz labels"""

        mtz_head = mtz_util.GetLabels(self.mtzin)
        if mtz_head is None:
            self.logger.error('No labels found, impossible to proceed with crank2 !!')
            self.error = True
            return

        missing_labels = []
        if mtz_head.iplus is None:
            missing_labels.append('I (+)')
        if mtz_head.sigiplus is None:
            missing_labels.append('SIGI (+)')
        if mtz_head.iminus is None:
            missing_labels.append('I (-)')
        if mtz_head.sigiminus is None:
            missing_labels.append('SIGI (-)')
        if mtz_head.free is None:
            missing_labels.append('FREE')

        if len(missing_labels) != 0:
            self.logger.error(
                'Following labels are missing in the mtzfile, cannot proceed with crank2 !!\n%s\n') % (
                ', '.join(missing_labels))
            self.error = True

    def get_scores(self, logfile=None):
        """Parse scores from the logfile"""

        for line in self.logcontents.split("\n"):
            if "R factor after refinement is" in line:
                self.rfactor = line.split(" ")[-1].rstrip()
            elif "R-free factor after refinement is" in line:
                self.rfree = line.split()[-1].rstrip()
            elif "FOM after" in line:
                self.FOM = line.split()[-1].rstrip()
        if self.FOM == "NA" or self.rfactor == "NA" or self.rfree == "NA":
            self.logger.warning("Parsing crank2 log file did not found all scores!")
            self.logger.warning("FOM: %s / Rfactor: %s / Rfree: %s" % (self.FOM, self.rfactor, self.rfree))
        elif float(self.rfree) < 0.45:
            self.solution = "YES"

    def run(self):
        """Run crank2"""

        # Check that the mtz has correct labels
        self._check_valid_mtz()
        if self.error:
            return

        # Manage dirs
        self.make_workdir()
        tmp_dir = os.getcwd()
        os.chdir(self.run_info.workdir)

        # Run crank2
        self._create_keyin()
        self.logger.debug(" ".join(self.cmd))
        self.logcontents = cexec(self.cmd, permit_nonzero=True)

        # Parse the output and exit
        self.get_scores()
        os.chdir(tmp_dir)
