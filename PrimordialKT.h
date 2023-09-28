#ifndef PRIMORDIALKT_H
#define PRIMORDIALKT_H

#include "Pythia8/Pythia.h"
namespace Pythia8{
  
// Class to handle rotations for the determination of the string axis
// when the momenta of the string endpoints are not in the same plane.
class StringRotation{
  public:
  // Empty constructor.
  StringRotation(){};
  // Constructor which needs the struck quark momentum Pq, the remnant momentum Pr,
  // the normal to the lepton scattering plane and the flag of a struck antiquark.
  // The momenta are required in the GNS system.
  StringRotation(Vec4 Pq, Vec4 Pr, Vec4 LeptonPlane_in, bool SeaAntiQuarkHit){
    // Determination of the normal vector to Pq and Pr.
    QuarkPlane = Vec4();
    if(SeaAntiQuarkHit)
      QuarkPlane = cross3(Pq,Pr);
    else
      QuarkPlane = cross3(Pr,Pq);
    QuarkPlane /= QuarkPlane.pAbs();
    // Normal to the lepton scattering plane.
    LeptonPlane = LeptonPlane_in;
    LeptonPlane /= LeptonPlane.pAbs();
    // Normal to the plane formed by QuarkPlane and LeptonPlane,
    // and relative angle.
    Normal = cross3(QuarkPlane,LeptonPlane);
    ThetaNormal = acos(dot3(QuarkPlane,LeptonPlane));
  }

  // Rotate a vector around the normal to the plane
  // formed by the normal vector to the lepton scattering plane
  // and the normal vector to the struck quark - remnant in the
  // GNS system.
  void Rotate(Vec4& v){
    v.rotaxis(ThetaNormal, Normal);}
  // Inverse of Rotate.
  void RotateBack(Vec4& v){
    v.rotaxis(-ThetaNormal, Normal);}

  public:
  // Normal vector to the plane formed by the momenta of
  // the string end points in the GNS.
  Vec4 QuarkPlane;
  // Normal to the lepton scattering plane in the GNS.
  Vec4 LeptonPlane;
  // Normal vector to the plane given by QuarkPlane x LeptonPlane.
  Vec4 Normal;
  // Angle between QuarkPlane and LeptonPlane.
  double ThetaNormal;
  
};

// End of PrimordialKT.h.
}
#endif

