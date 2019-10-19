import six
import math
import logging
import os
import pandas as pd
import stat
from cmsl1t.config import get_unique_out_dir
logger = logging.getLogger(__name__)


class Status:
    CREATED = 'CREATED'
    PENDING = 'PENDING'
    RUNNING = 'RUNNING'
    FAILED = 'FAILED'
    FINISHED = 'FINISHED'
    UNKNOWN = 'UNKNOWN'
    HELD = 'HELD'


class Batch:
    lsf = 'LSF'
    condor = 'HTCondor'


class PATHS:
    RUN_SCRIPT = 'run.sh'
    LOG = 'logs'
    CONFIGS = '_configs'
    INFO_FILE = 'info.csv'


def _get_info_file_path(batch_dir):
    return os.path.join(batch_dir, 'info.csv')


def get_config_name_template(config_file, batch_config_dir):
    batch_filename = os.path.basename(config_file.name)
    batch_filename = list(os.path.splitext(batch_filename))
    batch_filename.insert(1, "_{index}")
    batch_filename = "".join(batch_filename)
    batch_filename = os.path.join(batch_config_dir, batch_filename)
    return batch_filename


def create_run_script(setup_script, project_root, batch_dir):
    run_script_content = _get_run_script(setup_script, project_root)
    run_script = os.path.join(batch_dir, PATHS.RUN_SCRIPT)
    with open(run_script, 'w') as f:
        f.write(run_script_content)
    # make executable
    st = os.stat(run_script)
    os.chmod(run_script, st.st_mode | stat.S_IEXEC)
    return run_script


def create_info_file(info, batch_dir):
    info_file = os.path.join(batch_dir, PATHS.INFO_FILE)
    df = pd.DataFrame(info)
    df['run_script'] = PATHS.RUN_SCRIPT
    df['batch_dir'] = batch_dir
    df['config_file'] = df['config_file'].str.replace(batch_dir, '')
    df['job_log'] = df['job_log'].str.replace(batch_dir, '')
    df['output_folder'] = df['output_folder'].str.replace(batch_dir, '')
    df['stderr_log'] = df['stderr_log'].str.replace(batch_dir, '')
    df['stdout_log'] = df['stdout_log'].str.replace(batch_dir, '')
    df.to_csv(info_file, index=False)
    return info_file


def _get_run_script(setup_script, project_root, shared_fs=True):
    run_script_contents = [
        '#!/usr/bin/env bash',
        'pushd {project_root}',
        'source {setup_script}',
        'echo == start ENV ==',
        'env',
        'echo == END ENV ==',
        'cmsl1t -c "$1"',
        '',
    ]
    if not shared_fs:
        # infer project root
        pass
    run_script_contents = '\n'.join(run_script_contents)
    run_script_contents = run_script_contents.format(
        project_root=project_root,
        setup_script=setup_script,
    )
    return run_script_contents


def _prepare_input_file_groups(input_ntuples, files_per_job):
    file_lists = []
    current_list = []
    for infile in input_ntuples:
        if not infile.startswith("root:"):
            infile = os.path.realpath(infile)
        current_list.append(infile)

        # Is the current list full?
        if len(current_list) >= files_per_job:
            file_lists.append(current_list)
            current_list = []

    # Even if the last list had fewer files than needed, make sure to use this
    # too
    if current_list:
        file_lists.append(current_list)

    return file_lists


def prepare_output_folders(output_folder):
    batch_dir = os.path.join(output_folder, "batch")
    batch_dir = get_unique_out_dir(batch_dir)
    batch_config_dir = os.path.join(batch_dir, PATHS.CONFIGS)
    batch_log_dir = os.path.join(batch_dir, PATHS.LOG)
    logger.info("Batch config files will be placed under: " + batch_config_dir)
    os.makedirs(batch_config_dir)
    os.makedirs(batch_log_dir)
    return batch_dir, batch_config_dir, batch_log_dir


def prepare_jobs(config, batch_filename_template, outdir, files_per_job):
    job_generator = _prepare_jobs(
        config, batch_filename_template, outdir, files_per_job)
    job_configs, job_ids, output_folders = six.moves.zip(*job_generator)
    return job_configs, job_ids, output_folders


def _prepare_jobs(config, batch_filename_template, outdir, files_per_job):
    # Get the list of input files
    input_ntuples = config.get('input', 'files')
    input_ntuples = _prepare_input_file_groups(input_ntuples, files_per_job)

    n_jobs = len(input_ntuples)
    n_jobs_pad_width = int(math.log10(n_jobs)) + 1
    padding = "{{:0{}}}".format(n_jobs_pad_width)

    for i, in_files in enumerate(input_ntuples):
        padded_index = padding.format(i)

        # Reset the input file list
        config.config['input']['files'] = in_files

        # Reset the output directory
        # TODO: assumes shared_fs
        output_folder = outdir.format(index=padded_index)
        config.config['output']['folder'] = output_folder
        os.makedirs(output_folder)

        # Dump the config file
        batch_file = batch_filename_template.format(index=padded_index)
        config.dump(batch_file)

        yield (batch_file, padded_index, output_folder)
