'''
YAML Config parser for cmsl1t. The details of the configrations are described
in http://cms-l1t-analysis.readthedocs.io/en/latest/tutorial/configuration.html
'''
from __future__ import print_function
import yaml
import os
from datetime import datetime
import logging
from cmsl1t.utils import module
from copy import deepcopy
import re


logger = logging.getLogger(__name__)
TODAY = datetime.now().timetuple()


def get_unique_out_dir(outdir=None, revision=1):
    full_outdir = outdir + "-v{rev}".format(rev=revision)
    if os.path.isdir(full_outdir):
        return get_unique_out_dir(outdir, revision + 1)
    return full_outdir


def get_last_version_of(outdir):
    paths = resolve_file_paths([outdir + '*'])
    max_version = -1
    last_version_path = None
    version_re = re.compile(r".*-v(\d+)$")
    for path in paths:
        v_match = version_re.match(path)
        if v_match:
            version = v_match.group(1)
            if version > max_version:
                max_version = version
                last_version_path = path
    return last_version_path


def resolve_file_paths(paths):
    from cmsl1t.utils.root_glob import glob
    all_files = []
    for p in paths:
        all_files.extend(glob(p))
    return all_files


class ConfigParser(object):
    SECTIONS = ['general', 'input', 'analysis', 'output']

    def __init__(self):
        self.config = {}
        self.config_errors = []

    def read(self, input_file, reload_histograms=False, hist_files=None):
        cfg = yaml.load(input_file)
        self._read_config(cfg, reload_histograms, hist_files)

    def _read_config(self, cfg, reload_histograms=False, hist_files=None):
        cfg['general'] = dict(version=cfg['version'], name=cfg['name'])
        del cfg['version'], cfg['name']

        input_files = cfg['input']['files']
        try:
            input_files = resolve_file_paths(input_files)
        except Exception as e:
            msg = 'Could not resolve file paths:' + str(e)
            logger.exception(msg)
            raise IOError(msg)
        cfg['input']['files'] = input_files
        ntuple_map_file = 'config/ntuple_content.yaml'
        ntuple_map_file = cfg['input'].get("ntuple_map_file", ntuple_map_file)
        if ntuple_map_file:
            cfg['input']['ntuple_map_file'] = os.path.realpath(ntuple_map_file)

        self.config = cfg

        if not self.is_valid():
            msg = '\n'.join(self.config_errors)
            logger.exception(msg)
            raise IOError(msg)

        try:
            self.__fill_outdir_and_reload_files(reload_histograms, hist_files)
        except Exception as e:
            msg = 'Could not fill out output template: ' + str(e)
            logger.exception(msg)
            raise IOError(msg)

        self._reduce_scopes()

    def sections(self):
        return self.config.keys()

    def options(self, section):
        return self.config[section].keys()

    def get(self, *args):
        opt = self.config
        for arg in args:
            opt = opt[arg]
        return opt

    def try_get(self, *args, **kwargs):
        default = kwargs.pop("default", None)
        opt = self.config
        for arg in args:
            opt = opt.get(arg, None)
            if opt is None:
                return default
        return opt

    def is_valid(self):
        results = [self.validate_sections()]
        results += [self.validate_input_files()]
        results += [self.validate_analyzers()]
        results += [self.validate_producers()]
        return all(results)

    def validate_sections(self):
        expectedSections = sorted(ConfigParser.SECTIONS)
        sections = sorted(self.config.keys())
        hasValidSections = expectedSections == sections
        if hasValidSections is not True:
            msg = 'Configuration has missing or invalid sections.\n'
            msg += self.__section_format(
                self.__compare_sections(expectedSections, sections))
            self.config_errors += [msg]
        return hasValidSections

    def __section_format(self, sections):
        section_template = ['  - {}' for _ in sections]
        section_template = '\n'.join(section_template)
        return section_template.format(*sections)

    def __compare_sections(self, expected, current):
        expected = set(expected)
        current = set(current)
        missing = expected.difference(current)
        invalid = current.difference(expected)
        decorated_sections = []
        for e in expected:
            if e in missing:
                decorated_sections.append(e + ' (missing)')
            else:
                decorated_sections.append(e)
        for i in invalid:
            decorated_sections.append(i + ' (invalid)')
        return sorted(decorated_sections)

    def validate_input_files(self):
        input_files = self.config['input']['files']
        if not input_files:
            msg = "Could not find any existing files.\n"
            msg += "Given files:\n"
            msg += '\n'.join(self.config['input']['files'])
            self.config_errors += [msg]
        return input_files != []

    def validate_analyzers(self):
        return self.__validate_module_imports(['analysis', 'analyzers'])

    def validate_producers(self):
        return self.__validate_module_imports(['analysis', 'producers'])

    def __validate_module_imports(self, config_keys):
        modules = deepcopy(self.config)
        for key in config_keys:
            if key in modules:
                modules = modules[key]
            else:
                msg = "Cannot find required config section: "
                msg += "::".join(config_keys)
                self.config_errors.append(msg)
                return False
        if isinstance(modules, dict):
            msg, results = self.__validate_module_setup_dict(modules)
        elif isinstance(modules, list):
            msg, results = self.__validate_module_setup_list(modules)
        if msg:
            self.config_errors.append('\n'.join(msg))
        return all(results)

    def __validate_module_setup_dict(self, modules):
        msg = []
        results = []
        for name in modules.keys():
            m = modules[name]['module']
            if isinstance(m, dict):
                m = m.keys()[0]
            if not module.exists(m):
                msg += ['Module {0} does not exist!'.format(m)]
                results += [False]
            else:
                results += [True]
        return msg, results

    def __validate_module_setup_list(self, modules):
        msg = []
        results = []
        for m in modules:
            if isinstance(m, dict):
                m = m['module']
            if not module.exists(m):
                msg += ['Module {0} does not exist!'.format(m)]
                results += [False]
            else:
                results += [True]
        return msg, results

    def __repr__(self):
        return self.config.__repr__()

    def __fill_outdir_and_reload_files(self, reload_histograms, hist_files):
        cfg = self.config

        # Deduce what sort of reload we want:
        if reload_histograms:
            if hist_files:
                hist_files = resolve_file_paths(hist_files.split())
                cfg['input']['hist_files'] = hist_files
                if len(hist_files) > 1:
                    reload_histograms = "merge"
                else:
                    reload_histograms = "plot specific"
            else:
                reload_histograms = "plot last"
        cfg["input"]["reload_histograms"] = reload_histograms

        output_folder = self.try_get('output', 'folder')
        if not output_folder:
            if reload_histograms == "plot specific":
                output_folder = os.path.dirname(hist_files[0])
            else:
                template = os.path.join(*cfg['output']['template'])

                date = '{y}{m:02d}{d:02d}'.format(
                    y=TODAY.tm_year, m=TODAY.tm_mon, d=TODAY.tm_mday)
                sample_name = cfg['input']['sample']['name']
                trigger_name = cfg['input']['trigger']['name']
                run_number = cfg['input']['run_number']

                output_folder = template.format(
                    date=date, sample_name=sample_name, trigger_name=trigger_name,
                    run_number=run_number)

                # Find the version of this output dir to use
                if cfg['input']['reload_histograms'] == "plot last":
                    latest_version = get_last_version_of(output_folder)
                    if not latest_version:
                        msg = "Cannot find valid input histogram-file directory."
                        msg += " Looking for: " + output_folder
                        logger.error(msg)
                        raise IOError(msg)
                    search_path = os.path.join(latest_version, "*.root")
                    self.config['input']['hist_files'] = resolve_file_paths([search_path])
                else:
                    # Either merging multiple hists, or we're reading trees
                    # Essentially, this is a new analysis output
                    output_folder = get_unique_out_dir(output_folder)

        output_folder = os.path.realpath(output_folder)
        plots_folder = os.path.join(output_folder, "plots")
        cfg['output']['folder'] = output_folder
        cfg['output']['plots_folder'] = get_unique_out_dir(plots_folder)

    def describe(self):
        return __doc__

    def dump(self, out_filename):
        # Remove contents that weren't in the actual input config file
        config = deepcopy(self.config)

        del config['output']['plots_folder']
        del config['input']['reload_histograms']
        if 'hist_files' in config['input']:
            del config['input']['hist_files']

        # copy the "general" section back to the top-level
        config.update(config.get('general', {}))

        with open(out_filename, "w") as out_file:
            out_file.write(yaml.dump(config))

    def _reduce_scopes(self):
        '''
            Reduces config scopes for analyzers and producers
        '''
        cfg = self.config
        analyzers = self.get('analysis', 'analyzers')
        analyzers = [self.reduce_scope_for_analyzer(a) for a in analyzers]
        cfg['analysis']['analyzers'] = analyzers

        producers = self.get('analysis', 'producers')
        producers = [self.reduce_scope_for_producer(p) for p in producers]
        cfg['analysis']['producers'] = producers

    def reduce_scope_for_analyzer(self, analyzer_name):
        analyzer_spec = self.get('analysis', 'analyzers')
        if isinstance(analyzer_spec, list):
            # Already reduced to a list
            return [a for a in analyzer_spec if a == analyzer_name][-1]

        analyzer = analyzer_spec[analyzer_name]
        forbidden_local_settings = ['name', 'input_files']
        for s in forbidden_local_settings:
            if s in analyzer:
                logger.warn('Setting {0} is forbidden in analysis::analyzers::{1}'.format(s, analyzer_name))
                analyzer.pop(s)

        global_settings = dict(
            triggerName=self.get('input', 'trigger')['name'],
            lumiJson=self.try_get('input', 'lumi_json', default=''),
            # TODO: do better for legacy analyzer
            input_files=self.get('input', 'files'),
        )
        # move analysis section minus producers and analyzers into analyzer scope
        analysis = deepcopy(self.config['analysis'])
        analysis.pop('analyzers')
        analysis.pop('producers')
        global_settings.update(analysis)

        reduced_scope = {'name': analyzer_name}
        reduced_scope.update(global_settings)
        reduced_scope.update(analyzer)

        return reduced_scope

    def reduce_scope_for_producer(self, producer):
        producers_spec = self.get('analysis', 'producers')
        if isinstance(producers_spec, list):
            # Already reduced to a list
            return producer

        producer_dict = producers_spec[producer]
        reduced_scope = {'name': producer}
        reduced_scope.update(producer_dict)

        return reduced_scope


if __name__ == '__main__':
    config = ConfigParser()
    config.read('config/demo.yaml')
    print(config)
    print()
    print('Sections', config.sections())
    for section in config.sections():
        print('Section', section, 'has options', config.options(section))
        for option in config.options(section):
            print('>> {} = {}'.format(option, config.get(section, option)))
    print()
    print(config.is_valid())
