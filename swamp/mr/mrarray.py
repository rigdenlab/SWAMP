from pyjob import TaskFactory
from swamp.mr.mr import Mr
from swamp.mr.mrjob import MrJob


class MrArray(Mr):
    """An array of molecular replacement tasks to solve a given structure.

    This class implements data structures to hold all the MR tasks to be executed on a target. It implements functions
    to run and store results of these tasks, contained as instances of :py:obj:`~swamp.mr.mrjob.MrJob` instances.

    :param str id: unique identifier for the :py:obj:`~swamp.mr.mrarray.MrArray` instance
    :param str workdir: working directory where the :py:obj:`~swamp.mr.mrjob.MrJob` instances will be executed
    :param str target_mtz: target's mtz filename
    :param str target_fa: target's fasta filename
    :param str platform: platform where the array of tasks will be executed (default 'sge')
    :param str queue_name: name of the queue where the tasks should be submitted (default None)
    :param str queue_environment: queue environment where the tasks should be submitted (default None)
    :param str phased_mtz: target's mtz filename containing phase information (default None)
    :param int max_concurrent_nprocs: maximum number of concurrent tasks to be executed at any given time (default 1)
    :param int job_kill_time: kill time assigned to :py:obj:`~swamp.mr.mrjob.MrJob` instances (default None)
    :param `~swamp.logger.swamplogger.SwampLogger` logger: logging interface for the MR pipeline (default None)
    :param bool silent: if set to True the logger will not print messages
    :param int max_array_size: set the maximum permitted number of :py:obj:`pyjob.Scripts` instances in a submitted \
    :py:obj:`pyjob.ClusterTask` (default None)
    :ivar list results: A list with the figures of merit obtained after the completion of the pipeline
    :ivar bool error: True if errors have occurred at some point on the pipeline
    :ivar list job_list: A list of the :py:obj:`~swamp.mr.mrjob.MrJob` instances contained on this \
    :py:obj:`~swamp.mr.mrarray.MrArray` instance.
    :ivar dict job_dict: A dictionary of the :py:obj:`~swamp.mr.mrjob.MrJob` instances contained on this \
    :py:obj:`~swamp.mr.mrarray.MrArray` instance. Key corresponds with :py:attr:`swamp.mr.mrjob.MrJob.id`
    :ivar list scripts: List of :py:obj:`pyjob.Scripts` instances to be executed on this \
    :py:obj:`~swamp.mr.mrarray.MrArray` instance
    :ivar str shell_interpreter: Indicates shell interpreter to execute \
    :py:obj:`~swamp.mr.mrjob.MrJob` (default '/bin/bash')

    :example:

    >>> from swamp.mr import MrArray, MrJob
    >>> mr_array = MrArray('<id>', '<workdir>', '<target_mtz>', 'target_fasta>')
    >>> mr_array.add(MrJob('<id>', '<workdir>'))
    >>> print(mr_array)
    MrArray(id="<id>", njobs=1)
    >>> mr_array.run()
    """

    def __init__(self, id, workdir, target_mtz, target_fa, platform="sge", queue_name=None, logger=None,
                 max_array_size=None, queue_environment=None, phased_mtz=None, max_concurrent_nprocs=1,
                 job_kill_time=None, silent=False):

        super(MrArray, self).__init__(id, target_fa, target_mtz, workdir, phased_mtz=phased_mtz,
                                      logger=logger, silent=silent)

        self.init_params = locals()
        self.logger.info(self.pipeline_header.format('MR-ARRAY'))
        self.logger.info(self._inform_args(**self.init_params))
        self.max_concurrent_nprocs = max_concurrent_nprocs
        self.platform = platform
        self.job_kill_time = job_kill_time
        self.queue_name = queue_name
        self.queue_environment = queue_environment
        self.job_list = []
        self.job_dict = {}
        self.scripts = []
        self.shell_interpreter = "/bin/bash"
        self.max_array_size = max_array_size

    def __repr__(self):
        return '{}(id={}, njobs={})'.format(self.__class__.__name__, self.id, len(self.job_list))

    def __contains__(self, id):
        """True if there is a job with the given id"""
        return id in self.job_dict

    def __delitem__(self, id):
        """Remove a job with given id"""
        job = self[id]
        job.parent_array = None
        self.job_dict.pop(id)
        self.job_list.remove(job)

    def __getitem__(self, id):
        """Return the job with the given id"""
        if isinstance(id, slice):
            indexes_to_keep = set(range(*id.indices(len(self))))
            copy_to_return = self.copy()
            for i, job in enumerate(self):
                if i not in indexes_to_keep:
                    copy_to_return.remove(job.id)
            return copy_to_return
        elif isinstance(id, int):
            return self.job_list[id]
        else:
            return self.job_dict[id]

    def __iter__(self):
        """Iterate over the job list"""
        for job in self.job_list:
            yield job

    def __len__(self):
        """Return the number of jobs"""
        return len(self.job_list)

    def __reversed__(self):
        """Reversed list of jobs"""
        for job in reversed(self.job_list):
            yield job

    # ------------------ General properties ------------------

    @property
    def cleanup_dir_list(self):
        """List of directories to cleanup after completion of the pipeline (None)"""
        return None

    @property
    def _other_task_info(self):
        """A dictionary with the extra kwargs for :py:obj:`pyjob.TaskFactory`"""

        info = {'directory': self.workdir, 'shell': self.shell_interpreter}

        if self.platform == 'local':
            info['processes'] = self.max_concurrent_nprocs
        else:
            info['max_array_size'] = self.max_concurrent_nprocs
        if self.queue_environment is not None:
            info['environment'] = self.queue_environment
        if self.queue_name is not None:
            info['queue'] = self.queue_name
        if self.job_kill_time is not None:
            info['runtime'] = self.job_kill_time

        return info

    # ------------------ Methods ------------------

    def add(self, value):
        """Add an instance of :py:obj:`~swamp.mr.mrjob.MrJob` to the array. This includes both the MrJob object \
        and its :py:obj:`pyjob.Script` attribute.

        :argument value: :py:obj:`~swamp.mr.mrjob.MrJob` instance to be added \
        to the array for execution
        :raises TypeError: value is not an instance of :py:obj:`~swamp.mr.mrjob.MrJob`
        :raises ValueError: a :py:obj:`~swamp.mr.mrjob.MrJob` instance with the same \
        :py:attr:`swamp.mr.mrjob.MrJob.id` is already contained in the array
        """

        if not isinstance(value, MrJob):
            raise TypeError('Can only add MrJob instances to an MrArray!')

        if value.id in self:
            raise ValueError("MrJob %s defined twice!" % value)

        self.logger.debug('Registering job %s into the array' % value)

        value.parent_array = self

        self.job_list.append(value)
        self.job_dict[value.id] = value
        self.scripts.append(value.script)

    def run(self, store_results=False):
        """Send the array for execution in the HPC using :py:obj:`pyjob.TaskFactory`

        :argument bool store_results: Not implemented
        """

        self.logger.info("Sending the MR task array to the HPC for execution")

        if self.max_array_size is not None:
            all_scripts = [tuple(self.scripts[x:x + self.max_array_size]) for x in
                           range(0, len(self.scripts), self.max_array_size)]
        else:
            all_scripts = (self.scripts,)

        for idx, scripts in enumerate(all_scripts, 1):
            self.logger.info('Sending task array %s / %s' % (idx, len(all_scripts)))
            with TaskFactory(self.platform, scripts, **self._other_task_info) as task:
                task.name = 'swamp_%s' % idx
                task.run()

        self.logger.info('All tasks in the array have been completed!')
        self.logger.info('Retrieving results')

    def append_results(self):
        """Append the results obtained in each :py:obj:`~swamp.mr.mrjob.MrJob` instance listed at \
        :py:attr:`~swamp.mr.mrarray.MrArray.job_list` into :py:attr:`~swamp.mr.mr.Mr.results`"""

        for job in self.job_list:

            if job.results is not None:
                self.logger.debug('Recovering results of job %s' % job.id)
                self.results += job.results

            else:
                self.logger.debug('Cannot find any results for job %s' % job.id)
