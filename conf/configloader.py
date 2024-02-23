import json
import dict_digger

import os
import os.path as pth

_config_path = __file__ + os.sep + os.pardir + os.sep + "configuration.json"

with open(_config_path) as config_file:
    _configuration = json.load(config_file)

def get_configuration_map():
    return _configuration

def get_configuration_item(*args):
    return dict_digger.dig(_configuration, *args)
