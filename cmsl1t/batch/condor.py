import logging
import os
import socket
import re

import htcondor
from plumbum import local

from .common import Status, Batch

logger = logging.getLogger(__name__)

CONDOR_STATUS = [
    Status.CREATED,
    Status.PENDING,
    Status.RUNNING,
    Status.FINISHED,
    Status.HELD,
]


def submit(config_files, batch_directory, batch_log_dir, run_script):
    logger.info("Will submit {0} jobs".format(len(config_files)))
    schedd = htcondor.Schedd()
    results = []

    job_cfgs = []
    for i, cfg in enumerate(config_files):
        job_cfgs.append(__create_job_cfg(i, cfg, batch_directory, batch_log_dir, run_script))

    if 'cern.ch' in socket.gethostname():
        results = list(_submit_via_command_line(job_cfgs, config_files, batch_directory))
    else:
        with schedd.transaction() as txn:
            for job_cfg, cfg in zip(job_cfgs, config_files):
                result = __submit_one(txn, job_cfgs, cfg)
                results.append(result)

    return results


def __create_job_cfg(index, config_file, batch_directory, batch_log_dir, run_script):
    cfg = os.path.realpath(config_file)
    stderr_log = os.path.join(batch_log_dir, 'job_{0}.out'.format(index))
    stdout_log = os.path.join(batch_log_dir, 'job_{0}.out'.format(index))
    job_log = os.path.join(batch_log_dir, 'job_{0}.log'.format(index))
    environment = 'HOME={}'.format(os.environ["HOME"])
    return dict(
        executable=run_script,
        arguments=cfg,
        output=stdout_log,
        error=stderr_log,
        log=job_log,
        environment=environment,
    )


def __submit_one(txn, job_cfg, cfg):
    sub = htcondor.Submit(job_cfg)
    out = sub.queue(txn)
    return dict(
        batch_id=int(out),
        batch=Batch.condor,
        config_file=cfg,
        stderr_log=job_cfg['error'],
        stdout_log=job_cfg['output'],
        job_log=job_cfg['log'],
        status=Status.CREATED,
    )


def _submit_via_command_line(job_cfgs, config_files, batch_directory):
    condor_submit_file = os.path.join(batch_directory, 'job.submit')
    with open(condor_submit_file, 'w+') as f:
        content = '\n'
        for job_cfg in job_cfgs:
            content += '\n'.join([k + ' = ' + v for k, v in job_cfg.items()])
            content += '\n' + 'queue' + '\n'
        f.write(content)
    condor_submit = local['condor_submit']
    out = condor_submit(condor_submit_file)
    batch_ids = list(_parse_condor_submit_output(out))
    for i, job_cfg in enumerate(job_cfgs):
        yield dict(
            # batch_id=int(out),
            batch_id=batch_ids[i],
            batch=Batch.condor,
            config_file=config_files[i],
            stderr_log=job_cfg['error'],
            stdout_log=job_cfg['output'],
            job_log=job_cfg['log'],
            status=Status.CREATED,
        )


def _parse_condor_submit_output(out):
    # last_line = out.split('\n')[-1]
    important_line = ''
    for line in out.split('\n'):
        if 'submitted to cluster' in line:
            important_line = line
            break
    m = re.search(r"(\d+)( .* )(\d+)\.", important_line)
    n_jobs, _, condor_id = m.groups()
    for i in range(int(n_jobs)):
        yield condor_id + '.' + str(i)


def get_status(batch_id):
    schedd = htcondor.Schedd()
    status, exit_code = __status_from_schedd(batch_id, schedd)
    if status == Status.UNKNOWN:
        from_history = __status_from_history(batch_id, schedd)
        if from_history is None:
            return Status.UNKNOWN
        status, exit_code = from_history

    if exit_code is None or exit_code == 0:
        return status
    else:
        return Status.FAILED


def __status_from_schedd(batch_id, schedd):
    cluster, process = batch_id.split(".")
    query = 'ClusterId=={} && ProcId=={}'.format(cluster, process)
    query = schedd.query(query, ['JobStatus', 'ExitCode'])
    if not query or query is None:
        return Status.UNKNOWN, None
    for job in query:
        exit_code = job['ExitCode'] if 'ExitCode' in job else None
        status = CONDOR_STATUS[job['JobStatus']]
        return status, exit_code


def __status_from_history(batch_id, schedd):
    query = 'ClusterId=={0}'.format(batch_id)
    for job in schedd.history(query, ['JobStatus', 'ExitCode'], 1):
        exit_code = job['ExitCode'] if 'ExitCode' in job else None
        status = CONDOR_STATUS[job['JobStatus']]
        return status, exit_code


def resubmit(config_files, local_ids, batch_directory, batch_log_dir, run_script):
    logger.info("Will resubmit {0} jobs".format(len(config_files)))
    schedd = htcondor.Schedd()
    results = []
    with schedd.transaction() as txn:
        for i, cfg in zip(local_ids, config_files):
            result = __submit_one(
                txn, i, cfg, batch_directory, batch_log_dir, run_script
            )
            results.append(result)

    return results
