from HiggsAnalysis.CombinedLimit.PhysicsModel import *
 
 
### This is the base python class to measure DY AFB
 
class DY_AFB(PhysicsModel):
    def __init__(self):
        return

    def doParametersOfInterest(self):
        """Create POI and other parameters, and define the POI set."""

        self.modelBuilder.doVar("Afb[0.6,-0.75,0.75]");
        self.modelBuilder.doVar("Rbk_ee[0.15,0.0,1.0]");
        self.modelBuilder.doVar("RQCD_ee[0.05,0.0,1.0]");
        self.modelBuilder.doVar("Rbk_mumu[0.15,0.0,1.0]");
        self.modelBuilder.doVar("RQCD_mumu[0.05,0.0,1.0]");
        self.modelBuilder.doSet("POI","Afb,Rbk_ee,Rbk_mumu,RQCD_ee,RQCD_mumu")
      
        self.modelBuilder.factory_('expr::Rpl_ee("(1.-@1-@2)*(1.+@0)",Afb, Rbk_ee, RQCD_ee)')
        self.modelBuilder.factory_('expr::Rmn_ee("(1.-@1-@2)*(1.-@0)",Afb, Rbk_ee, RQCD_ee)')

        self.modelBuilder.factory_('expr::Rpl_mumu("(1.-@1-@2)*(1.+@0)",Afb, Rbk_mumu, RQCD_mumu)')
        self.modelBuilder.factory_('expr::Rmn_mumu("(1.-@1-@2)*(1.-@0)",Afb, Rbk_mumu, RQCD_mumu)')

 
 
 
 
    def getYieldScale(self,bin,process):
        if 'pl' in process and 'ee' in bin: return "Rpl_ee"
        if 'mn' in process and 'ee' in bin: return "Rmn_ee"
        if 'pl' in process and 'mumu' in bin: return "Rpl_mumu"
        if 'mn' in process and 'mumu' in bin: return "Rmn_mumu"

        if 'bk' in process and 'ee' in bin: return "Rbk_ee"
        if 'bk' in process and 'mumu' in bin: return "Rbk_mumu"
        if 'qcd' in process and 'ee' in bin: return "RQCD_ee"
        if 'qcd' in process and 'mumu' in bin: return "RQCD_mumu"
        else:
            return 1

class DY_AFB_noQCD(PhysicsModel):
    def __init__(self):
        return

    def doParametersOfInterest(self):
        """Create POI and other parameters, and define the POI set."""

        self.modelBuilder.doVar("Afb[0.6,-0.75,0.75]");
        self.modelBuilder.doVar("Rbk_ee[0.15,0.0,1.0]");
        self.modelBuilder.doVar("Rbk_mumu[0.15,0.0,1.0]");
        self.modelBuilder.doSet("POI","Afb,Rbk_ee,Rbk_mumu")
      
        self.modelBuilder.factory_('expr::Rpl_ee("(1.-@1)*(1.+@0)",Afb, Rbk_ee)')
        self.modelBuilder.factory_('expr::Rmn_ee("(1.-@1)*(1.-@0)",Afb, Rbk_ee)')

        self.modelBuilder.factory_('expr::Rpl_mumu("(1.-@1)*(1.+@0)",Afb, Rbk_mumu)')
        self.modelBuilder.factory_('expr::Rmn_mumu("(1.-@1)*(1.-@0)",Afb, Rbk_mumu)')

 
 
 
 
    def getYieldScale(self,bin,process):
        if 'pl' in process and 'ee' in bin: return "Rpl_ee"
        if 'mn' in process and 'ee' in bin: return "Rmn_ee"
        if 'pl' in process and 'mumu' in bin: return "Rpl_mumu"
        if 'mn' in process and 'mumu' in bin: return "Rmn_mumu"

        if 'bk' in process and 'ee' in bin: return "Rbk_ee"
        if 'bk' in process and 'mumu' in bin: return "Rbk_mumu"
        else:
            return 1

class ElEl_samesign(PhysicsModel):
    def __init__(self):
        return

    def doParametersOfInterest(self):
        """Create POI and other parameters, and define the POI set."""

        self.modelBuilder.doVar("Rdy_ee[1.,0.0,10.0]");
        self.modelBuilder.doVar("RQCD_ee[1,0.0,10.0]");
        self.modelBuilder.doSet("POI","Rdy_ee,RQCD_ee")
      
 
 
 
    def getYieldScale(self,bin,process):
        if 'dy' in process : return "Rdy_ee"
        if 'qcd' in process : return "RQCD_ee"
        else:
            return 1

elel_samesign = ElEl_samesign()
dy_AFB = DY_AFB() 
dy_AFB_noQCD = DY_AFB_noQCD() 
