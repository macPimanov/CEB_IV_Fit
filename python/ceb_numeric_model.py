import numpy as np
from scipy import integrate
from .constants import PhysicsConstants

class CEBNumericModel:
    def __init__(self):
        self.constants = PhysicsConstants()
    
    def current(self, v, tau):
        a0 = 1.0 + 0.375 * tau - 0.1171875 * np.power(tau, 2)
        a1 = tau * np.power(a0, 2)
        a2 = 1.0 + np.exp((np.abs(v) - 1.0) / tau - (1.15 + tau))
        a3 = (np.power(v, 2) - 1.0) / (1.0 + np.exp(-(np.power(v, 2) - 1.0) / tau))
        
        s = (v >= 0.0) * 2.0 - 1.0
        return s * np.sqrt(2.0 * np.pi * a1 / a2 + a3) * (1.0 / (2.0 * np.exp((1.0 - v) / tau) + 1.0) - 1.0 / (2.0 * np.exp((1.0 + v) / tau) + 1.0))
    
    def current_integral(self, DT, v, tau, tauE):
        dx = self.constants.INTEGRATION_SCALE * DT
        inf = self.constants.ESSENTIALLY_INFINITY_SCALE * DT
        DT2 = np.power(DT, 2)
        v1 = np.abs(v)
        
        i = 0.0
        x = DT + dx
        
        while x <= inf:
            accuracy = 1e-15
            a2 = np.sqrt(np.power(x, 2) - DT2)
            a0 = 1.0 / (np.exp((x - v1) / tauE) + 1.0)
            a1 = 1.0 / (np.exp(x / tau) + 1.0)
            i_step = x / a2 * (a0 - a1)
            
            if np.abs(i_step) < accuracy:
                break
            i += i_step
            
            a0 = 1.0 / (np.exp((-x - v1) / tauE) + 1.0)
            a1 = 1.0 / (np.exp(-x / tau) + 1.0)
            i_step = x / a2 * (a0 - a1)
            
            if np.abs(i_step) < accuracy:
                break
            i += i_step
            
            x += dx
        
        s = (v >= 0.0) * 2.0 - 1.0
        return i * dx * s
    
    def power_cool_piece(self, v, tau):
        PItau = np.pi * tau
        a0 = (1.0 - v) / tau
        tmp1 = 2.0 * np.exp(a0)
        tmp2 = np.exp(-2.5 * (a0 + 2.0))
        a1 = np.sqrt(2.0 * PItau) * ((1.0 - v) / (tmp1 + 1.28) + 0.5 * tau / (tmp1 + 0.64)) / (tmp2 + 1.0)
        
        if v - 1.0 > tau:
            a2 = np.sqrt(np.power(v, 2) - 1.0)
            return a1 + 0.5 * (-v * a2 + np.log(v + a2) + np.power(PItau, 2) / 3.0 * v / a2) / (1.0 / tmp2 + 1.0)
        return a1
    
    def power_cool(self, v, tau, tauE):
        a1 = self.power_cool_piece(v, tauE)
        b1 = self.power_cool_piece(-v, tauE)
        c1 = self.power_cool_piece(0.0, tau)
        return a1 + b1 - c1
    
    def and_current(self, DT, v, tauE, Wt, tm):
        dx = self.constants.INTEGRATION_SCALE * DT
        DT2 = np.power(DT, 2)
        twoWt = 2.0 * Wt
        twoTauE = 2.0 * tauE
        
        i = 0.0
        x = 0.0
        
        while x <= DT:
            DT2x2 = DT2 - np.power(x, 2)
            sqrtDT2x2 = np.sqrt(DT2x2)
            da0a1 = np.tanh((x + v) / twoTauE) - np.tanh((x - v) / twoTauE)
            a2 = twoWt * sqrtDT2x2 / tm / (np.power(x * (twoWt - sqrtDT2x2 / DT), 2) + DT2x2 / np.power(DT * tm, 2))
            i += DT / sqrtDT2x2 * da0a1 * a2
            x += dx
        
        return i * dx
    
    def power_cool_integral(self, DT, v, tau, tauE):
        dx = self.constants.INTEGRATION_SCALE * DT
        inf = self.constants.ESSENTIALLY_INFINITY_SCALE * DT
        v1 = np.abs(v)
        DT2 = np.power(DT, 2)
        
        po = 0.0
        ps = 0.0
        x = DT + dx
        
        while x <= inf:
            x2 = np.power(x, 2)
            a = np.sqrt(x2 - DT2)
            
            da2a1 = (1.0 / (np.exp((x - v1) / tauE) + 1.0) - 1.0 / (np.exp(x / tau) + 1.0)) / a
            po += x * (x - v1) * da2a1
            ps += x2 * da2a1
            
            da2a1 = (1.0 / (np.exp((-x - v1) / tauE) + 1.0) - 1.0 / (np.exp(-x / tau) + 1.0)) / a
            po += x * (-x - v1) * da2a1
            ps += -x2 * da2a1
            
            x += dx
        
        po *= dx
        ps *= dx
        return po, ps