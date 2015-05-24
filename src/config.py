'''
Configuration file reading and access functions.

The configuration file is written in YAML and is supplied
by the user.

TODO: validation of config file input.
'''

import yaml


class Config(object):
    def __init__(self, config_filename):
        # Try to open and parse the YAML formatted config file
        with open(config_filename) as config_file:
            try:
                config = yaml.load(config_file)
            except yaml.YAMLError, exc:
                print("Error in configuration file:", exc)
                raise exc
        self.config = config
        self.config_filename = config_filename


    def get_option(self, option):
        '''Retrieve a global option from the configuration'''
        if option in self.config:
            return self.config[option]
        else:
            raise Exception("Unknown option: {}, not in configuration " \
                "file: {}".format(option, self.config_filename))


    def get_stage_option(self, stage, option):
        '''Try to retrieve a configuration option for a particular stage.
        If the stage does not define the option then look for it in the
        default options.
        '''
        stages = self.config['stages']
        if stage in stages:
            this_stage = stages[stage]
            # Look for the option in the specified stage
            if option in this_stage:
                return this_stage[option]
            else:
                # Look for the option in the defaults
                defaults = self.config['defaults']
                if option in defaults:
                    return defaults[option]
                else:
                    # Option does not have a default value
                    raise Exception("Option: {} not defined in config for " \
                        "stage: {} nor in defaults in configuration " \
                        "file {}".format(option, stage, self.config_filename))
        else:
            # Stage does not exist in the config file
            raise Exception("Unknown stage: {}, not in configuration " \
                "file: {}".format(stage, self.config_filename))


    def validate(self):
        '''Check that the configuration is valid.'''
        config = self.config
        filename = self.config_filename
        # Test for required fields: defaults, stages, fastqs, read_groups
        check_required_field(config, filename, 'defaults')
        check_required_field(config, filename, 'stages')
        check_required_field(config, filename, 'fastqs')
        check_required_field(config, filename, 'read_groups')


    def get_read_group(self, sample):
        '''Get the read group information for a given sample'''
        config = self.config
        filename = self.config_filename
        # This is already validated to exist
        read_groups = config['read_groups']
        if sample in read_groups:
            return read_groups[sample]
        else:
            raise Exception("Configuration file {} does not have read group " \
                "for sample {}".format(filename, sample))


def check_required_field(config, filename, field):
    '''Utility to check whether a field exists in the config dictionary'''
    if field  not in config:
        raise Exception("Configuration file {} does not have '{}' " \
            "field".format(filename, field))
