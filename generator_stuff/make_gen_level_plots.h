
#include "TObject.h"
class TGenParticle : public TObject
{
    public:
        TGenParticle():
            parent(-1), pdgId(0),status(0),
            pt(0), eta(0), phi(0), mass(0), y(0)
    {}
        ~TGenParticle(){}

        int   parent;
        int   pdgId;
        int   status;
        float pt, eta, phi, mass, y;

        ClassDef(TGenParticle,1)
};

class TGenEventInfo : public TObject
{
    public:
        TGenEventInfo():
            id_1(0),  id_2(0),
            x_1(0),   x_2(0),
            scalePDF(0), xs(0), weight(1)
    {}
        ~TGenEventInfo(){}

        int   id_1,  id_2;   // parton flavor PDG ID
        float x_1,   x_2;    // parton momentum fraction
        float scalePDF;      // Q-scale used for PDF evaluation
        float xs;            // cross section from LHE file
        float weight;        // generator-level event weight

        ClassDef(TGenEventInfo,3)
};
