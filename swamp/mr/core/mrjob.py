import os
import dill
from pyjob import Script
import swamp.mr.core.mrarray


class MrJob(object):
    """Class to manage the creation and execution of a :py:obj:`~swamp.mr.core.mrrun.MrRun` in the context of a parent
     container :py:obj:`~swamp.mr.core.mrarray.MrArray` instance.

     This class implements methods to create a python script that can be executed as single independent job.
     It also implements utilities and data structures to retrieve and store the results obtained with the resulting
     instance of :py:obj:`~swamp.mr.core.mrrun.MrRun`.

     :param str id: unique identifier of this :py:obj:`~swamp.mr.core.mrjob.MrJob` isntance
     :param str workdir: working directory for :py:obj:`~swamp.mr.core.mrjob.MrJob` instance
     :param str python_interpreter: python interpreter for the :py:obj:`pyjob.Script`
     :ivar str id: Unique identifier of this :py:obj:`~swamp.mr.core.mrjob.MrJob` isntance
     :ivar str phased_mtz: target's mtz filename containing phases (default None)
     :ivar str target_mtz: target's mtz filename (default None)
     :ivar str target_fa: target's fasta filename (default None)
     :ivar list searchmodel_list: A list with the search models to be used in the \
     :py:obj:`~swamp.mr.core.mrrun.MrRun` instance
     """

    def __init__(self, id, workdir, python_interpreter=os.path.join(os.environ['CCP4'], 'bin', 'ccp4-python')):

        self.id = id
        self.workdir = workdir
        self.python_interpreter = python_interpreter
        self.target_mtz = None
        self.target_fa = None
        self.phased_mtz = None
        self.parent_array = None
        self.searchmodel_list = []

        if not os.path.isdir(self.workdir):
            os.makedirs(self.workdir)

    def __repr__(self):
        return '{}(id={}, workdir="{}")'.format(self.__class__.__name__, self.id, self.workdir)

    # ------------------ General properties ------------------

    @property
    def parent_array(self):
        """The parent :py:obj:`~swamp.mr.core.mrarray.MrArray` instance"""
        return self._parent_array

    @parent_array.setter
    def parent_array(self, value):
        """Property setter for :py:attr:`~swamp.mr.core.mrjob.MrJob.parent_array`

        :param value: MrArray to be set
        :type value: :py:obj:`~swamp.mr.core.mrarray.MrArray`
        :raises TypeError: value is not an instance of :py:obj:`~swamp.mr.core.mrarray.MrArray`
        """

        if value is None:
            pass
        elif not isinstance(value, swamp.mr.core.mrarray.MrArray):
            raise TypeError('Parent array must be a swamp.mr.mrarray.MrArray instance!')
        else:
            self._parent_array = value
            self.target_mtz = self.parent_array.target_mtz
            self.target_fa = self.parent_array.target_fa
            self.phased_mtz = self.parent_array.phased_mtz

    @property
    def results(self):
        """A nested list with the results of the :py:obj:`~swamp.mr.core.mrrun.MrRun` instance created with the \
        execution of :py:attr:`~swamp.mr.core.mrjob.MrJob.python_script`"""

        pickle_fname = os.path.join(self.workdir, "results.pckl")

        if os.path.isfile(pickle_fname):
            with open(pickle_fname, 'rb') as pickle_fhandle:
                mr_run = dill.load(pickle_fhandle)
            pickle_fhandle.close()
            return mr_run.results

        else:
            return None

    @property
    def python_script(self):
        """String with the python script to create and execute the :py:obj:`~swamp.mr.core.mrrun.MrRun` instance \
        associated with this :py:obj:`~swamp.mr.core.mrjob.MrJob`"""

        script = """cd {_workdir}
{_python_interpreter} << EOF
from swamp.mr.core.mrrun import MrRun
mr_run = MrRun(id='{_id}', workdir='{_workdir}', target_fa='{_target_fa}', target_mtz='{_target_mtz}')\n""".format(
            **self.__dict__)

        if self.phased_mtz is not None:
            script += 'mr_run.phased_mtz = "%s"\n' % self.phased_mtz

        for searchmodel in self.searchmodel_list:
            args_list = []
            for arg in searchmodel.keys():
                args_list.append('%s="%s"' % (arg, searchmodel[arg]))
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
        """A :py:obj:`pyjob.Script` instance that will be executed on this :py:obj:`~swamp.mr.core.mrjob.MrJob`"""

        script = Script(directory=os.path.join(os.path.abspath(os.path.join(self.workdir, os.pardir))),
                        prefix=self.id.lower(), stem='', suffix='.sh')
        script.append(self.python_script)
        return script

    # ------------------ Methods ------------------

    def add_searchmodel(self, **kwargs):
        """Provide necessary information to add a given search model to the :py:obj:`~swamp.mr.core.mrrun.MrRun` \
        instance.

        :param kwargs: directly used into :py:func:`~swamp.mr.core.mrrun.MrRun.add_searchmodel`
        """

        self.searchmodel_list.append(kwargs)
