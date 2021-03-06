#ifndef BPHParticleMasses_H
#define BPHParticleMasses_H
/** \class BPHParticleMasses
 *
 *  Description: 
 *     Class to contain particle masses
 *
 *
 *  $Date: 2016-08-02 15:44:55 $
 *  $Revision: 1.1 $
 *  \author Paolo Ronchese INFN Padova
 *
 */

//----------------------
// Base Class Headers --
//----------------------


//------------------------------------
// Collaborating Class Declarations --
//------------------------------------


//---------------
// C++ Headers --
//---------------


//              ---------------------
//              -- Class Interface --
//              ---------------------

class BPHParticleMasses {

 public:

  static const double        elMass;
  static const double      muonMass;
  static const double       tauMass;
  static const double      pionMass;
  static const double      kaonMass;
  static const double       kx0Mass;
  static const double       phiMass;
  static const double      jPsiMass;
  static const double      psi2Mass;
  static const double      ups1Mass;
  static const double      ups2Mass;
  static const double      ups3Mass;
  static const double    protonMass;
  static const double   lambda0Mass;
  static const double   lambdaxMass;
  static const double      penQMass;
  static const double Lambda0_bMass;
  
  static const double        elMSigma;
  static const double      muonMSigma;
  static const double       tauMSigma;
  static const double      pionMSigma;
  static const double      kaonMSigma;
  static const double    protonMSigma;
  static const double      jPsiMSigma;
  static const double Lambda0_bMSigma;

  static const double       kx0MWidth;
  static const double       phiMWidth;
  static const double      jPsiMWidth;
  static const double      psi2MWidth;
  static const double      ups1MWidth;
  static const double      ups2MWidth;
  static const double      ups3MWidth;
  static const double   lambda0MWidth;
  static const double   lambdaxMWidth;
  static const double      penQMWidth;
  static const double Lambda0_bMWidth;

 private:

};


#endif // BPHParticleMasses_H

