import os
import dill
from pyjob import Script
import swamp.mr.mrarray


class MrJob(object):
    """Class to manage the creation and exection of a MrRun in the context of a parent :obj:`swamp.mr.mrarray.MrArray`

     This class implements methods to create an executable python script that can be submitted as single HPC job.
     It also implements utilities to retrieve results for the resulting :obj:`swamp.mr.mrrun.MrRun` instance.

     :param str id: the id given to the MrRun and identifying this job
     :param str workdir: the working directory for the MrRun
     :ivar list searchmodel_list: a list with the searchmodels to be used in the MrRun
     :ivar :obj:`~swamp.mr.mrarray.MrArray` parent_array: acquires the value of the parent MrArray instance
     :ivar str target_mtz: target's mtz filename
     :ivar str target_fa: target's fasta filename
     :ivar str phased_mtz: target's mtz filename containing phases (default None)
     :ivar str python_interpreter: location of the python interpreter for the script
     """

    def __init__(self, id, workdir, python_interpreter=os.path.join(os.environ['CCP4'], 'bin', 'ccp4-python')):

        self._id = id
        self._workdir = workdir
        self._python_interpreter = python_interpreter
        self._target_mtz = None
        self._target_fa = None
        self._phased_mtz = None
        self._parent_array = None
        self._searchmodel_list = []

    def __repr__(self):
        return '{}(id={}, workdir="{}")'.format(self.__class__.__name__, self.id, self.workdir)

    # ------------------ General properties ------------------

    @property
    def id(self):
        return self._id

    @id.setter
    def id(self, value):
        self._id = value

    @property
    def python_interpreter(self):
        return self._python_interpreter

    @python_interpreter.setter
    def python_interpreter(self, value):
        self._python_interpreter = value

    @property
    def workdir(self):
        return self._workdir

    @workdir.setter
    def workdir(self, value):
        self._workdir = value

    @property
    def phased_mtz(self):
        return self._phased_mtz

    @phased_mtz.setter
    def phased_mtz(self, value):
        self._phased_mtz = value

    @property
    def target_mtz(self):
        return self._target_mtz

    @target_mtz.setter
    def target_mtz(self, value):
        self._target_mtz = value

    @property
    def target_fa(self):
        return self._target_fa

    @target_fa.setter
    def target_fa(self, value):
        self._target_fa = value

    @property
    def searchmodel_list(self):
        return self._searchmodel_list

    @searchmodel_list.setter
    def searchmodel_list(self, value):
        self._searchmodel_list = value

    @property
    def parent_array(self):
        return self._parent_array

    @parent_array.setter
    def parent_array(self, value):
        """Setter of parent array value

        :param value: MrArray to be set
        :type value: :obj:`~swamp.mr.mr.mrarray.MrArray
        :raises TypeError: value is not an instance of :obj:`~swamp.mr.mr.mrarray.MrArray
        """

        if value is None:
            pass
        elif not isinstance(value, swamp.mr.mrarray.MrArray):
            raise TypeError('Parent array must be a swamp.mr.mrarray.MrArray instance!')
        else:
            self._parent_array = value
            self.target_mtz = self.parent_array.target_mtz
            self.target_fa = self.parent_array.target_fa
            self.phased_mtz = self.parent_array.phased_mtz

    @property
    def results(self):
        """Results of the MrRun resulting from the execution of the python script

        :returns results: a nested list with the results of the MrRun resulting from the execution of the python script
        :rtype list, None
        """

        pickle_fname = os.path.join(self.workdir, "results.pckl")

        if os.path.isfile(pickle_fname):
            results = []
            with open(pickle_fname, 'rb') as pickle_fhandle:
                mr_run = dill.load(pickle_fhandle)
            pickle_fhandle.close()
            results += mr_run.results
            return results

        else:
            return None

    @property
    def _python_script(self):
        """Python script to create and execute the corresponding :obj:`swamp.mr.mrrun.MrRun` instance"""

        script = """cd {_workdir}
{_python_interpreter} << EOF
from swamp.mr.mrrun import MrRun
mr_run = MrRun(id='{_id}', workdir='{_workdir}', target_fa='{_target_fa}', target_mtz='{_target_mtz}')\n""".format(
            **self.__dict__)

        if self.phased_mtz is not None:
            script += 'mr_run.phased_mtz = "%s"\n' % self.phased_mtz

        for searchmodel in self.searchmodel_list:
            args_list = []
            for arg in searchmodel.keys():
                args_list.append('%s=%s' % (arg, searchmodel[arg]))
            script += 'mr_run.add_searchmodel(%s)\n' % ', '.join(args_list)

        script += """if not mr_run.error:
    mr_run.run()
    mr_run.create_result_table_outfile()
    mr_run.store_pickle()
EOF
"""
        return script

    @property
    def script(self):
        """Instance of :object:`pyjob.Script` that corresponds with the job to be executed"""

        script = Script(directory=os.path.join(os.path.abspath(os.path.join(self.workdir, os.pardir))),
                        prefix=self.id.lower(), stem='', suffix='.sh')
        script.append(self._python_script)
        return script

    # ------------------ Methods ------------------

    def add_searchmodel(self, **kwargs):
        """Add the necessary information to add a given search model to the MrRun

        :param kwargs: the arguments will be passed directly to :obj:`~swamp.mr.mrrun.add_searchmodel()`
        :returns nothing
        :rtype None
        """

        self.searchmodel_list.append(kwargs)
