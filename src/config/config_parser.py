'''
Config parser

This file serves as the loader and reader of the yaml config file

Functions
---------
parse_config
'''

import yaml

def parse_config() -> dict:
    '''
    Parses config.yml file and returns its content as a dictionary.

    Returns
    -------
    config : dict
        config dictionary 
    '''

    with open('src/config/config.yml', 'r') as yml:
        config = yaml.safe_load(yml)
    return config

config = parse_config()