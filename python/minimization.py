import numpy as np
from scipy.optimize import minimize, brent
from lmfit import Parameters, minimize as lmfit_minimize, report_fit
from typing import Tuple, Callable, Any, List, Optional
from typing import Sequence

class MinimizationAlgorithms:
    @staticmethod
    def golden_section_minimization(f: Callable[[float], float], a: float, b: float, tolerance: float = 1e-8) -> Tuple[float, float]:
        R = (np.sqrt(5.0) - 1.0) / 2.0
        C = 1.0 - R
        
        x1 = b - R * (b - a)
        x2 = a + R * (b - a)
        f1 = f(x1)
        f2 = f(x2)
        
        while abs(b - a) > tolerance * (abs(x1) + abs(x2)):
            if f2 < f1:
                a = x1
                x1 = x2
                f1 = f2
                x2 = a + R * (b - a)
                f2 = f(x2)
            else:
                b = x2
                x2 = x1
                f2 = f1
                x1 = b - R * (b - a)
                f1 = f(x1)
        
        x_min = x1 if f1 < f2 else x2
        f_min = min(f1, f2)
        return x_min, f_min
    
    @staticmethod
    def brent_minimization(f: Callable[[float], float], a: float, b: float, tolerance: float = 1e-8) -> Tuple[float, float]:
        result = brent(f, brack=(a, b), tol=tolerance, full_output=True)
        return result[0], result[1]
    
    @staticmethod
    def scipy_minimize(f: Callable[[np.ndarray], float], x0: np.ndarray, bounds: Optional[List[Tuple[float, float]]] = None, method: str = 'L-BFGS-B', tolerance: float = 1e-8) -> Tuple[np.ndarray, float]:
        result = minimize(f, x0, method=method, bounds=bounds, tol=tolerance)
        return result.x, result.fun
    
    @staticmethod
    def lmfit_minimize(params: Parameters, objective_func: Callable[[Parameters], Any], method: str = 'leastsq', tolerance: float = 1e-8) -> Tuple[Parameters, float]:
        result = lmfit_minimize(objective_func, params, method=method, tol=tolerance)
        return result.params, result.chisqr
    
    @staticmethod
    def differential_brent(f: Callable[[float], float], df: Callable[[float], float], a: float, b: float, tolerance: float = 1e-8) -> Tuple[float, float]:
        def combined_obj(x: float) -> float:
            fx = f(x)
            dfx = df(x)
            return fx
        
        result = brent(combined_obj, brack=(a, b), tol=tolerance, full_output=True)
        return result[0], result[1]
    
    @staticmethod
    def bracket_minimum(f: Callable[[float], float], a: float, b: float, factor: float = 2.0) -> Tuple[float, float, float, float, float, float]:
        fa = f(a)
        fb = f(b)
        
        if fb > fa:
            a, b = b, a
            fa, fb = fb, fa
        
        c = b + factor * (b - a)
        fc = f(c)
        
        while fb > fc:
            TINY = 1.0e-20
            GLIMIT = 100.0
            
            r = (b - a) * (fb - fc)
            q = (b - c) * (fb - fa)
            
            denom = max(abs(q - r), TINY)
            u = b - ((b - c) * q - (b - a) * r) / (2.0 * np.sign(q - r) * denom)
            
            ulim = b + GLIMIT * (c - b)
            
            if (b > u) == (u > c):
                fu = f(u)
                if fu < fc:
                    a = b
                    b = u
                    fa, fb = fb, fu
                    return a, b, c, fa, fb, fc
                elif fu > fb:
                    c = u
                    fc = fu
                    return a, b, c, fa, fb, fc
                u = c + factor * (c - b)
                fu = f(u)
            elif (c > u) == (u > ulim):
                fu = f(u)
                if fu < fc:
                    b, c, u = c, u, u + factor * (u - c)
                    fb, fc, fu = fc, fu, f(u)
            elif (u >= ulim) == (ulim >= c):
                u = ulim
                fu = f(u)
            else:
                u = c + factor * (c - b)
                fu = f(u)
            
            a, fa = b, fb
            b, fb = c, fc
            c, fc = u, fu
        
        return a, b, c, fa, fb, fc