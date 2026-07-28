import numpy as np
from scipy.constants import e, h, k, pi

class PhysicsConstants:
    ME = 9.1093837015e-31
    E = e
    H = h  
    HBAR = h / (2 * pi)
    K = k / e
    BCS_INTEGRAL = 1.764
    NUMBER_OF_SINS_IN_CEB = 2.0
    INTEGRATION_SCALE = 0.2e-3
    ESSENTIALLY_INFINITY_SCALE = 20.0
    
    @staticmethod
    def get_amplifier_constants(amp_type='AD745'):
        amplifiers = {
            'AD745': {'voltage_noise': 1.0e-8, 'current_noise': 6.9e-15},
            'OPA111': {'voltage_noise': 8.0e-9, 'current_noise': 0.8e-15},
            'AD797': {'voltage_noise': 0.9e-9, 'current_noise': 2000.0e-15},
            'IFN146': {'voltage_noise': 1.1e-9, 'current_noise': 0.3e-15},
            'OPA1641': {'voltage_noise': 5.1e-9, 'current_noise': 0.8e-15},
        }
        return amplifiers.get(amp_type, amplifiers['AD745'])

class Configuration:
    def __init__(self, config_file=None):
        self.data_file = "SPC-CEB_300mK_Triton11-2026.txt"
        self.amp_type = 'AD745'
        self.sep = '\t'
        self.threads = 56
        
        if config_file:
            self.load_from_json(config_file)
    
    def load_from_json(self, config_file):
        import json
        with open(config_file, 'r') as f:
            config = json.load(f)
        
        for key, value in config.items():
            setattr(self, key, value)
    
    def to_json(self, config_file):
        import json
        config_dict = {key: value for key, value in vars(self).items() 
                      if not key.startswith('_')}
        with open(config_file, 'w') as f:
            json.dump(config_dict, f, indent=4)