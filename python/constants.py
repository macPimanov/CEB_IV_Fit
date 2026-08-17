from scipy.constants import e, h, k, pi
from typing import Dict

class PhysicsConstants:
    ME: float = 9.1093837015e-31
    E: float = e
    H: float = h  
    HBAR: float = h / (2 * pi)
    K: float = k / e
    BCS_INTEGRAL: float = 1.764
    NUMBER_OF_SINS_IN_CEB: float = 2.0
    INTEGRATION_SCALE: float = 0.2e-3
    ESSENTIALLY_INFINITY_SCALE: float = 20.0
    
    @staticmethod
    def get_amplifier_constants(amp_type: str = 'AD745') -> Dict[str, float]:
        amplifiers: Dict[str, Dict[str, float]] = {
            'AD745': {'voltage_noise': 1.0e-8, 'current_noise': 6.9e-15},
            'OPA111': {'voltage_noise': 8.0e-9, 'current_noise': 0.8e-15},
            'AD797': {'voltage_noise': 0.9e-9, 'current_noise': 2000.0e-15},
            'IFN146': {'voltage_noise': 1.1e-9, 'current_noise': 0.3e-15},
            'OPA1641': {'voltage_noise': 5.1e-9, 'current_noise': 0.8e-15},
        }
        return amplifiers.get(amp_type, amplifiers['AD745'])
