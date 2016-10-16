class Logspace:
    def __init__(self):
        self.LOGZERO = 'LOGZERO'
    def eexp(self,x):
        if x == self.LOGZERO:
            return 0
        else:
            return np.exp(x)
    def eln(self,x):
        if x == 0:
            return self.LOGZERO
        elif x>0:
            return np.log(x)
        else:
            print 'Wrong!!!\n\t negative input error'
            return np.nan
    def elnsum(self,elnx,elny):
        if elnx == self.LOGZERO:
            return elny
        elif elny == self.LOGZERO:
            return elnx
        elif elnx > elny:
            return elnx + self.eln(1+np.exp(elny-elnx))
        else:
            return elny + self.eln(1+np.exp(elnx-elny))
    def elnproduct(self,elnx,elny):
        if elnx == self.LOGZERO or elny == self.LOGZERO:
            return self.LOGZERO
        else:
            return elnx + elny
            
            
            
logspace = Logspace()

elnx = logspace.eln(np.exp(60))
elny = logspace.eln(np.exp(-10))
elnx_plus_y = logspace.elnsum(elnx,elny)
print elnx_plus_y
elnx_prod_y = logspace.elnproduct(elnx,elny)
print elnx_prod_y
logzero = logspace.eln(0)
elnx_plus_y = logspace.elnsum(elnx,logzero)
print elnx_plus_y
elnx_prod_y = logspace.elnproduct(elnx,logzero)
print elnx_prod_y
print logspace.eln(-1)
