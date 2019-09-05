from HiggsAnalysis.CombinedLimit.PhysicsModel import *
 
#this version measures AFB and gauge boson polarization (alpha)
class DY_AFB(PhysicsModel):
    def __init__(self):
        return

    def doParametersOfInterest(self):
        """Create POI and other parameters, and define the POI set."""

        self.modelBuilder.doVar("Afb[0.6,-0.75,0.75]");
        self.modelBuilder.doVar("A0[0.05, -1.0, 1.0]");
        self.modelBuilder.doSet("POI","Afb,A0")

        # ss templates
        self.modelBuilder.doVar("Rdy_mumu_ss[1.0,0.0,10.0]");
        self.modelBuilder.doVar("Rdy_ee_ss[1.0,0.0,10.0]");
      
        self.modelBuilder.factory_('expr::Alph("2.0*@0/(2.0-@0)",A0)')
        self.modelBuilder.factory_('expr::Norm("3.0/(4.0*(2.0+@0))",Alph)')
        self.modelBuilder.factory_('expr::RAlph("@0*@1",Alph,Norm)')
        self.modelBuilder.factory_('expr::Rpl("(@0+@1)",Norm,Afb)')
        self.modelBuilder.factory_('expr::Rmn("(@0-@1)",Norm,Afb)')


 
 
 
 
    def getYieldScale(self,bin,process):
        if 'alpha' in process: return "RAlph"
        if 'pl' in process: return "Rpl"
        if 'mn' in process : return "Rmn"
        if 'dy' in process and 'ee_ss' in bin: return "Rdy_ee_ss"
        if 'dy' in process and 'mumu_ss' in bin: return "Rdy_mumu_ss"
        else:
            return 1
 
### This is class measures AFB only
class DY_AFB_LO(PhysicsModel):
    def __init__(self):
        return

    def doParametersOfInterest(self):
        """Create POI and other parameters, and define the POI set."""

        self.modelBuilder.doVar("Afb[0.6,-0.75,0.75]");
        self.modelBuilder.doSet("POI","Afb")

        # ss templates
        self.modelBuilder.doVar("Rdy_mumu_ss[1.0,0.0,10.0]");
        self.modelBuilder.doVar("Rdy_ee_ss[1.0,0.0,10.0]");
      
        self.modelBuilder.factory_('expr::Rpl("(1.+@0)",Afb)')
        self.modelBuilder.factory_('expr::Rmn("(1.-@0)",Afb)')


 
 
 
 
    def getYieldScale(self,bin,process):
        if 'pl' in process: return "Rpl"
        if 'mn' in process : return "Rmn"
        if 'dy' in process and 'ee_ss' in bin: return "Rdy_ee_ss"
        if 'dy' in process and 'mumu_ss' in bin: return "Rdy_mumu_ss"
        else:
            return 1

class DY_AFB_dilu(PhysicsModel):
    def __init__(self):
        return

    def doParametersOfInterest(self):
        """Create POI and other parameters, and define the POI set."""

        self.modelBuilder.doVar("Afb[0.6,-0.75,0.75]");
        self.modelBuilder.doSet("POI","Afb")

        # ss templates
        self.modelBuilder.doVar("Dilu_ratio[1.0,0.0,10.0]");
        self.modelBuilder.doVar("Rdy_mumu_ss[1.0,0.0,10.0]");
        self.modelBuilder.doVar("Rdy_ee_ss[1.0,0.0,10.0]");
      
        self.modelBuilder.factory_('expr::Rpl("0.5*(1.+@0*@1)",Afb, Dilu_ratio)')
        self.modelBuilder.factory_('expr::Rmn("0.5*(1.-@0*@1)",Afb, Dilu_ratio)')


 
 
 
 
    def getYieldScale(self,bin,process):
        if 'pl' in process: return "Rpl"
        if 'mn' in process : return "Rmn"
        if 'dy' in process and 'ee_ss' in bin: return "Rdy_ee_ss"
        if 'dy' in process and 'mumu_ss' in bin: return "Rdy_mumu_ss"
        else:
            return 1

class Samesign(PhysicsModel):
    def __init__(self):
        return

    def doParametersOfInterest(self):
        """Create POI and other parameters, and define the POI set."""

        self.modelBuilder.doVar("Rdy[1.,0.0,10.0]");
        self.modelBuilder.doVar("Rqcd[1,0.0,10.0]");
        self.modelBuilder.doSet("POI","Rdy,Rqcd")
      
 
 
 
    def getYieldScale(self,bin,process):
        if 'dy' in process : return "Rdy"
        if 'qcd' in process : return "Rqcd"
        else:
            return 1

class EMu(PhysicsModel):
    def __init__(self):
        return

    def doParametersOfInterest(self):
        """Create POI and other parameters, and define the POI set."""

        self.modelBuilder.doVar("Rdy[1.,0.0,10.0]");
        self.modelBuilder.doVar("Rbk[1.,0.0,10.0]");
        self.modelBuilder.doVar("Rqcd_emu[1,0.0,10.0]");
        self.modelBuilder.doSet("POI","Rbk,Rdy,Rqcd_emu")
      
 
 
 
    def getYieldScale(self,bin,process):
        if 'bk' in process : return "Rbk"
        if 'dy' in process : return "Rdy"
        if 'qcd' in process : return "Rqcd_emu"
        else:
            return 1
samesign = Samesign()
emu = EMu()
dy_AFB = DY_AFB() 
dy_AFB_LO = DY_AFB_LO() 
dy_AFB_dilu = DY_AFB_dilu() 
