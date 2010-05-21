///////////////////////////////////////////////////////////////////////////
//
//    Copyright 2010
//
//    This file is part of rootpwa
//
//    rootpwa is free software: you can redistribute it and/or modify
//    it under the terms of the GNU General Public License as published by
//    the Free Software Foundation, either version 3 of the License, or
//    (at your option) any later version.
//
//    rootpwa is distributed in the hope that it will be useful,
//    but WITHOUT ANY WARRANTY; without even the implied warranty of
//    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//    GNU General Public License for more details.
//
//    You should have received a copy of the GNU General Public License
//    along with rootpwa. If not, see <http://www.gnu.org/licenses/>.
//
///////////////////////////////////////////////////////////////////////////
//-------------------------------------------------------------------------
// File and Version Information:
// $Rev::                             $: revision of last commit
// $Author::                          $: author of last commit
// $Date::                            $: date of last commit
//
// Description:
//      class hierarchy for mass-dependent part of the amplitude
//
//
// Author List:
//      Boris Grube          TUM            (original author)
//
//
//-------------------------------------------------------------------------


#ifndef MASSDEPENDENCE_H
#define MASSDEPENDENCE_H


#include <iostream>
#include <complex>

#include <boost/shared_ptr.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/triangular.hpp>
#include <boost/numeric/ublas/lu.hpp>


namespace ublas = boost::numeric::ublas;


namespace rpwa {

  class isobarDecayVertex;
  typedef boost::shared_ptr<isobarDecayVertex> isobarDecayVertexPtr;


  class massDependence {

  public:
  
    massDependence()          { }
    virtual ~massDependence() { }

    virtual std::complex<double> amp(const isobarDecayVertex& v) = 0;

    virtual std::complex<double> operator ()(const isobarDecayVertex& v) { return amp(v); }

    virtual std::ostream& print(std::ostream& out) const;

    static bool debug() { return _debug; }                             ///< returns debug flag
    static void setDebug(const bool debug = true) { _debug = debug; }  ///< sets debug flag
    

  private:
    
    static bool _debug;  ///< if set to true, debug messages are printed
    
  };


  typedef boost::shared_ptr<massDependence> massDependencePtr;


  inline
  std::ostream&
  operator <<(std::ostream&         out,
	      const massDependence& massDep)
  {
    return massDep.print(out);
  }


  class flatMassDependence : public massDependence {

  public:
  
    flatMassDependence()          { }
    virtual ~flatMassDependence() { }

    virtual std::complex<double> amp(const isobarDecayVertex&);

    virtual std::ostream& print(std::ostream& out) const;

    static bool debug() { return _debug; }                             ///< returns debug flag
    static void setDebug(const bool debug = true) { _debug = debug; }  ///< sets debug flag
    

  private:
    
    static bool _debug;  ///< if set to true, debug messages are printed
    
  };


  typedef boost::shared_ptr<flatMassDependence> flatMassDependencePtr;


  inline
  flatMassDependencePtr
  createFlatMassDependence()
  {
    flatMassDependencePtr md(new flatMassDependence());
    return md;
  }


  class relativisticBreitWigner : public massDependence {

  public:

    relativisticBreitWigner()          { }
    virtual ~relativisticBreitWigner() { }

    virtual std::complex<double> amp(const isobarDecayVertex& v);

    virtual std::ostream& print(std::ostream& out) const;

    static bool debug() { return _debug; }                             ///< returns debug flag
    static void setDebug(const bool debug = true) { _debug = debug; }  ///< sets debug flag
    

  private:
    
    static bool _debug;  ///< if set to true, debug messages are printed
    
  };


  typedef boost::shared_ptr<relativisticBreitWigner> relativisticBreitWignerPtr;


  inline
  relativisticBreitWignerPtr
  createRelativisticBreitWigner()
  {
    relativisticBreitWignerPtr md(new relativisticBreitWigner());
    return md;
  }


  class piPiSWaveAuMorganPenningtonM : public massDependence {

  public:

    piPiSWaveAuMorganPenningtonM();
    virtual ~piPiSWaveAuMorganPenningtonM() { }

    virtual std::complex<double> amp(const isobarDecayVertex& v);

    virtual std::ostream& print(std::ostream& out) const;

    static bool debug() { return _debug; }                             ///< returns debug flag
    static void setDebug(const bool debug = true) { _debug = debug; }  ///< sets debug flag

    int vesSheet;


  protected:

    ublas::matrix<std::complex<double> >               _T;
    std::vector<ublas::matrix<std::complex<double> > > _a;
    std::vector<ublas::matrix<std::complex<double> > > _c;
    ublas::matrix<double>                              _sP;


  private:

    double _piChargedMass;
    double _piNeutralMass;
    double _kaonChargedMass;
    double _kaonNeutralMass;
    double _kaonMeanMass;
    
    static bool _debug;  ///< if set to true, debug messages are printed
    
  };


  typedef boost::shared_ptr<piPiSWaveAuMorganPenningtonM> piPiSWaveAuMorganPenningtonMPtr;


  inline
  piPiSWaveAuMorganPenningtonMPtr
  createPiPiSWaveAuMorganPenningtonM()
  {
    piPiSWaveAuMorganPenningtonMPtr md(new piPiSWaveAuMorganPenningtonM());
    return md;
  }


  // /** @brief AMP parameterization of pipi s-wave
  //  *  
  //  *  We have introduced a small modification by setting the off-diagonal 
  //  *  elements of the M-matrix to zero.
  //  */
  // class AMP_M : public massDependence {

  // public:

  //   AMP_M();
  //   virtual ~AMP_M() { }
  //   AMP_M(const AMP_M&) { }
  //   virtual massDependence* create() const { return new AMP_M(); }
  //   virtual massDependence* clone() const { return new AMP_M(*this); }

  //   virtual void print() { std::cout << "AMP_M"; }
  //   std::complex<double> val(const particle& p);
  // };


  // /** @brief old VES parameterization
  //  *
  //  *  Brute force subtraction of the f0(980)
  //  */ 
  // class AMP_ves : public AMP_M {
  // public:
  //   AMP_ves() : AMP_M() { ves_sheet = 1; }
  //   virtual ~AMP_ves() { }
  //   AMP_ves(const AMP_ves&) { }
  //   virtual massDependence* create() const { return new AMP_ves(); }
  //   virtual massDependence* clone() const { return new AMP_ves(*this); }

  //   virtual void print() { std::cout << "AMP_ves"; }
  //   std::complex<double> val(const particle& p);
  // };


  // /** @brief Kachaev's version of the AMP parameterization
  //  * 
  //  * From the original fortran code:
  //  * Source: K.L.Au et al, Phys.Rev. D35, P 1633. M solution.
  //  * 04-Mar-2003 See eps_k1.for for description.
  //  * Here matrix M=K^{-1} is parametrized with one pole.
  //  * Misprint in the article (other than in K1--K3 solutions)
  //  * was corrected.
  //  *
  //  * 14-Mar-2003 Nice amplitude for pi-pi S-wave without f0(975).
  //  * It is smooth and nicely tends to zero after approx 1.5 GeV.
  //  * f0(975) pole excluded; coupling to KK zeroed; set C411=C422=0.
  //  * The largest effect from C411, zeroing of C422 looks insignificant.
  //  */ 
  // class AMP_kach : public AMP_M {
  // public:
  //   AMP_kach();
  //   virtual ~AMP_kach() { }
  //   AMP_kach(const AMP_kach&) { }
  //   virtual massDependence* create() const { return new AMP_kach(); }
  //   virtual massDependence* clone() const { return new AMP_kach(*this); }

  //   virtual void print() { std::cout << "AMP_kach"; }
  //   //std::complex<double> val(const particle& p);
  // };


  // http://www.crystalclearsoftware.com/cgi-bin/boost_wiki/wiki.pl?LU_Matrix_Inversion
  // matrix inversion routine using lu_factorize and lu_substitute
  template<class T>
  bool
  invertMatrix(const ublas::matrix<T>& A,
	       ublas::matrix<T>&       inverseA)
  {
    // create working copy of input
    ublas::matrix<T> M(A);
    // create permutation matrix for LU-factorization
    ublas::permutation_matrix<std::size_t> pM(M.size1());
    // perform LU-factorization
    if (ublas::lu_factorize(M, pM) != 0)
      return false;
    // create identity matrix of "inverse"
    inverseA.assign(ublas::identity_matrix<T>(M.size1()));
    // backsubstitute to get the inverse
    ublas::lu_substitute(M, pM, inverseA);
    return true;
  }


  template<class T>
  ublas::matrix<T>
  invertMatrix(const ublas::matrix<T>& A,
	       bool&                   isSingular)
  { 
    ublas::matrix<T> inverseA (A.size1(), A.size2());
    isSingular = !invert (A, inverseA);
    return inverseA;
  }


}  // namespace rpwa


#endif  // MASSDEPENDENCE_H
