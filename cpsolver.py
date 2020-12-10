# -*- coding: utf-8 -*-
import numpy as np

class cpsolver:
    def __init__(self, 
                 func, 
                 N = 64,
                 thr = 0.1, 
                 margin = 1e-7,
                 err = 1e-7,
                 ):
        self.F = func
        self.N = N
        self.thr = thr
        self.margin = margin
        self.err = err
        self.Res = set()
        
    def reset(self):
        self.Res = set()
        
    def proc(self, Rr , Ir):
        D0 = np.hstack((np.linspace(Rr[0] + Ir[0], Rr[1] + Ir[0], self.N),
                        np.linspace(Rr[1] + Ir[0], Rr[1] + Ir[1], self.N),
                        np.linspace(Rr[1] + Ir[1], Rr[0] + Ir[1], self.N),
                        np.linspace(Rr[0] + Ir[1], Rr[0] + Ir[0], self.N)))
        
        D1 = np.array(list(map(self.F,D0)))
        
        if 0j in D1:S = 1
        else:S = np.abs(np.sum(np.diff(D1) / D1[:-1]).imag / (2*np.pi))
        
        if S > self.thr:
            Rc = (Rr[0] + Rr[1]) / 2
            Ic = (Ir[0] + Ir[1]) / 2
            
            Val = Rc + Ic
            Span = np.abs(Rr[0] - Rr[1] + Ir[0] - Ir[1])
            Error = np.abs(self.F(Val))
            
            if  Span < self.margin and  Error < self.err:
                self.Res.add((Val, Span, Error))
            else:
                self.proc(Rr=(Rr[0],Rc),Ir=(Ir[0],Ic))
                self.proc(Rr=(Rr[0],Rc),Ir=(Ic,Ir[1]))
                self.proc(Rr=(Rc,Rr[1]),Ir=(Ir[0],Ic))
                self.proc(Rr=(Rc,Rr[1]),Ir=(Ic,Ir[1]))
            
    
    def solve(self, Rr = (-5, 5), Ir = (-5j, 5j)) :
        self.proc(Rr = (-5, 5),Ir = (-5j, 5j))
                
        for i in self.Res:print('Main Value={:.3},Margin={:.3},Error={:.3}'.format(*i))
