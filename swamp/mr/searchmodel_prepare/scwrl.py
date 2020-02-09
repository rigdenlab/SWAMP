import os
import subprocess
from pyjob import cexec
from swamp.wrappers.wrapper import Wrapper


class Scwrl(Wrapper):
    """Scwrl wrapper

    Examples
    --------
        Working on it...
    """

    def __init__(self, work_dir, pdbin, pdbout, seqin):
        self.work_dir = work_dir
        self.pdbin = pdbin
        self.seqin = seqin
        self.pdbout = pdbout
        self.logfile = os.path.join(self.work_dir, "scwrl.log")
        package_path = os.path.dirname(os.path.dirname(os.path.dirname(__file__)))
        self.scwrl4_source = os.path.join(package_path, "lib", "scwrl4", "Scwrl4")

        super(Scwrl, self).__init__(run_info, target_info, searchmodel_info)

    # Create a method to run scwrl
    def run(self):
        # Get in the workdir
        if not os.path.isdir(self.work_dir):
            os.mkdir(self.work_dir)
        os.chdir(self.work_dir)

        # Strip the side chains from the original pdbin
        pdbin_resseq = os.path.join(self.work_dir, "%s_resseq.pdb" % os.path.basename(self.pdbin)[:-4])
        cmd = "pdbset xyzin %s xyzout %s" % (self.pdbin, pdbin_resseq)
        print(cmd)
        p = subprocess.Popen(cmd, stdin=subprocess.PIPE, stdout=subprocess.PIPE, shell=True)
        stdout, stderr = p.communicate(input="exclude side" + os.linesep + "chain A" + os.linesep + "end")
        # Create the log file
        with open(self.work_dir + "/pdbset.log", 'w') as f_out:
            f_out.write(stdout)
        # Resequence the file
        resequence_pdb_file(pdbin=pdbin_resseq, resseq_fasta=self.seqin)
        # Renumber the file
        renumber_residues(pdbin=pdbin_resseq, in_situ=True)

        # Load the command and run it
        cmd = [self.scwrl4_source, "-i", pdbin_resseq, "-o", self.pdbout, "-t", "-h"]
        print(" ".join(cmd))
        stdout = cexec(cmd, permit_nonzero=False)
        # Create the log file
        with open(self.logfile, 'w') as f_out:
            f_out.write(stdout)
