/* FILE: kl_uniformexchange.h            -*-Mode: c++-*-
 *
 * Uniform 6 neighbor exchange energy on rectangular mesh,
 * derived from Oxs_Energy class.
 * Based on file: uniformexchange.h
 * Extensions done by Kristof Lebecki, Universität Konstanz (Feb 2010):
 * - periodic boundary onditions in z-dimension
 * - taking into account non-uniform magnetization, caused by temperature
 *
 */

#ifndef _KLM_UNIFORMEXCHANGE
#define _KLM_UNIFORMEXCHANGE

#include "oc.h"
#include "director.h"
#include "chunkenergy.h"
#include "energy.h"
#include "meshvalue.h"
#include "simstate.h"
#include "threevector.h"
#include "util.h"

/* End includes */

class Klm_UniformExchange:public Oxs_ChunkEnergy {
private:
  enum ExchangeCoefType {
    A_UNKNOWN, A_TYPE, LEX_TYPE
  }  excoeftype;
  OC_REAL8m A;
  OC_REAL8m lex;

  enum ExchangeKernel { NGBR_UNKNOWN,
			NGBR_6_FREE,
                        NGBR_6_MIRROR, NGBR_6_MIRROR_STD,
                        NGBR_6_BIGANG_MIRROR, NGBR_6_ZD2,
// KL(m) Both are based on NGBR_6_MIRROR, i.e. CalcEnergy6NgbrMirror
// (eventually on CalcEnergy6NgbrMirror_lex)
			NGBR_6_Z_PERIOD, // KL(m) z-periodic conditions.
			NGBR_6_LLB,      // KL(m) LLB: influence of temperature
			NGBR_6_LLB_Z_P,  // KL(m) LLB + z-periodicity
			NGBR_12_FREE, NGBR_12_ZD1, NGBR_12_ZD1B,
			NGBR_12_MIRROR,	NGBR_26 } kernel;
  /// Exchange formulation to use.  "unknown" is invalid; it
  /// is defined for error detection.
  /// NOTE: "kernel" is set inside the constructor, and should
  ///  be fixed thereafter.
  
  // "?coef" and "?integ" are coefficient matrices used by some
  // of the CalcEnergy routines.  (Well, currently just 12NgbrZD1.)
  // They are "mutable" so they can be changed from inside the
  // (const) CalcEnergy routines.
  mutable Oxs_Mutex thread_mutex;
  mutable OC_UINT4m mesh_id;
  mutable Nb_2DArrayWrapper<OC_REAL8m> xcoef;
  mutable Nb_2DArrayWrapper<OC_REAL8m> ycoef;
  mutable Nb_2DArrayWrapper<OC_REAL8m> zcoef;
  mutable OC_REAL8m xinteg[3];
  mutable OC_REAL8m yinteg[3];
  mutable OC_REAL8m zinteg[3];

  mutable OC_REAL8m energy_density_error_estimate; // Cached value,
                                 /// initialized when mesh changes.

  void InitCoef_12NgbrZD1(OC_INDEX size,
			  OC_REAL8m wgt[3],
			  Nb_2DArrayWrapper<OC_REAL8m>& coef) const;
  /// Sets "coef" and "integ" arrays for use by NGBR_12_ZD1 kernel.

  // Support for threaded mindot calculations
  mutable vector<OC_REAL8m> mindot;

  // Utility routine for CalcEnergy6NgbrBigAngleMirror
  OC_REAL8m ComputeAngle(const ThreeVector& u1,const ThreeVector& u2) const;

  // Calculation routines for each of the
  // aforementioned energy formulations.
  void CalcEnergy6NgbrFree
  (const Oxs_MeshValue<ThreeVector>& spin,
   const Oxs_MeshValue<OC_REAL8m>& Ms_inverse,
   const Oxs_RectangularMesh* mesh,
   Oxs_ComputeEnergyDataThreaded& ocedt,
   Oxs_ComputeEnergyDataThreadedAux& ocedtaux,
   OC_INDEX node_start,OC_INDEX node_stop,
   int threadnumber) const;
  void CalcEnergy6NgbrMirror
  (const Oxs_MeshValue<ThreeVector>& spin,
   const Oxs_MeshValue<OC_REAL8m>& Ms_inverse,
   const Oxs_RectangularMesh* mesh,
   Oxs_ComputeEnergyDataThreaded& ocedt,
   Oxs_ComputeEnergyDataThreadedAux& ocedtaux,
   OC_INDEX node_start,OC_INDEX node_stop,
   int threadnumber) const;
  void CalcEnergy6NgbrMirror_lex
  (const Oxs_MeshValue<ThreeVector>& spin,
   const Oxs_MeshValue<OC_REAL8m>& Ms,
   const Oxs_RectangularMesh* mesh,
   Oxs_ComputeEnergyDataThreaded& ocedt,
   Oxs_ComputeEnergyDataThreadedAux& ocedtaux,
   OC_INDEX node_start,OC_INDEX node_stop,
   int threadnumber) const;
  void CalcEnergy6NgbrMirrorStd
  (const Oxs_MeshValue<ThreeVector>& spin,
   const Oxs_MeshValue<OC_REAL8m>& Ms_inverse,
   const Oxs_RectangularMesh* mesh,
   Oxs_ComputeEnergyDataThreaded& ocedt,
   Oxs_ComputeEnergyDataThreadedAux& ocedtaux,
   OC_INDEX node_start,OC_INDEX node_stop,
   int threadnumber) const;
  void CalcEnergy6NgbrBigAngMirror
  (const Oxs_MeshValue<ThreeVector>& spin,
   const Oxs_MeshValue<OC_REAL8m>& Ms_inverse,
   const Oxs_RectangularMesh* mesh,
   Oxs_ComputeEnergyDataThreaded& ocedt,
   Oxs_ComputeEnergyDataThreadedAux& ocedtaux,
   OC_INDEX node_start,OC_INDEX node_stop,
   int threadnumber) const;
  void CalcEnergy6NgbrZD2
  (const Oxs_MeshValue<ThreeVector>& spin,
   const Oxs_MeshValue<OC_REAL8m>& Ms_inverse,
   const Oxs_RectangularMesh* mesh,
   Oxs_ComputeEnergyDataThreaded& ocedt,
   Oxs_ComputeEnergyDataThreadedAux& ocedtaux,
   OC_INDEX node_start,OC_INDEX node_stop,
   int threadnumber) const;
  //KL(m)
  void CalcEnergy6NgbrMirror_zperiodic
  (const Oxs_MeshValue<ThreeVector>& spin,
   const Oxs_MeshValue<OC_REAL8m>& Ms_inverse,
   const Oxs_RectangularMesh* mesh,
   Oxs_ComputeEnergyDataThreaded& ocedt,
   Oxs_ComputeEnergyDataThreadedAux& ocedtaux,
   OC_INDEX node_start,OC_INDEX node_stop,
   int threadnumber) const;
  void CalcEnergy6NgbrMirror_lex_zperiodic
  (const Oxs_MeshValue<ThreeVector>& spin,
   const Oxs_MeshValue<OC_REAL8m>& Ms,
   const Oxs_RectangularMesh* mesh,
   Oxs_ComputeEnergyDataThreaded& ocedt,
   Oxs_ComputeEnergyDataThreadedAux& ocedtaux,
   OC_INDEX node_start,OC_INDEX node_stop,
   int threadnumber) const;
  void CalcEnergy6NgbrMirror_LLB
  (const Oxs_MeshValue<ThreeVector>& spin,
   const Oxs_MeshValue<OC_REAL8m>& Ms_inverse,
   const Oxs_MeshValue<OC_REAL8m>& Ms, // KL(m) we need it as well
   const Oxs_MeshValue<OC_REAL8m>& Me_T, // KL(m)
   const Oxs_RectangularMesh* mesh,
   Oxs_ComputeEnergyDataThreaded& ocedt,
   Oxs_ComputeEnergyDataThreadedAux& ocedtaux,
   OC_INDEX node_start,OC_INDEX node_stop,
   int threadnumber) const;
  void CalcEnergy6NgbrMirror_LLB_zp
  (const Oxs_MeshValue<ThreeVector>& spin,
   const Oxs_MeshValue<OC_REAL8m>& Ms_inverse,
   const Oxs_MeshValue<OC_REAL8m>& Ms, // KL(m) we need it as well
   const Oxs_MeshValue<OC_REAL8m>& Me_T, // KL(m)
   const Oxs_RectangularMesh* mesh,
   Oxs_ComputeEnergyDataThreaded& ocedt,
   Oxs_ComputeEnergyDataThreadedAux& ocedtaux,
   OC_INDEX node_start,OC_INDEX node_stop,
   int threadnumber) const;
  //
  void CalcEnergy12NgbrFree
  (const Oxs_MeshValue<ThreeVector>& spin,
   const Oxs_MeshValue<OC_REAL8m>& Ms_inverse,
   const Oxs_RectangularMesh* mesh,
   Oxs_ComputeEnergyDataThreaded& ocedt,
   Oxs_ComputeEnergyDataThreadedAux& ocedtaux,
   OC_INDEX node_start,OC_INDEX node_stop,
   int threadnumber) const;
  void CalcEnergy12NgbrZD1
  (const Oxs_MeshValue<ThreeVector>& spin,
   const Oxs_MeshValue<OC_REAL8m>& Ms_inverse,
   const Oxs_RectangularMesh* mesh,
   Oxs_ComputeEnergyDataThreaded& ocedt,
   Oxs_ComputeEnergyDataThreadedAux& ocedtaux,
   OC_INDEX node_start,OC_INDEX node_stop,
   int threadnumber) const;
  void CalcEnergy12NgbrZD1B
  (const Oxs_MeshValue<ThreeVector>& spin,
   const Oxs_MeshValue<OC_REAL8m>& Ms_inverse,
   const Oxs_RectangularMesh* mesh,
   Oxs_ComputeEnergyDataThreaded& ocedt,
   Oxs_ComputeEnergyDataThreadedAux& ocedtaux,
   OC_INDEX node_start,OC_INDEX node_stop,
   int threadnumber) const;
  void CalcEnergy12NgbrMirror
  (const Oxs_MeshValue<ThreeVector>& spin,
   const Oxs_MeshValue<OC_REAL8m>& Ms_inverse,
   const Oxs_RectangularMesh* mesh,
   Oxs_ComputeEnergyDataThreaded& ocedt,
   Oxs_ComputeEnergyDataThreadedAux& ocedtaux,
   OC_INDEX node_start,OC_INDEX node_stop,
   int threadnumber) const;
  void CalcEnergy26Ngbr
  (const Oxs_MeshValue<ThreeVector>& spin,
   const Oxs_MeshValue<OC_REAL8m>& Ms_inverse,
   const Oxs_RectangularMesh* mesh,
   Oxs_ComputeEnergyDataThreaded& ocedt,
   Oxs_ComputeEnergyDataThreadedAux& ocedtaux,
   OC_INDEX node_start,OC_INDEX node_stop,
   int threadnumber) const;

  // Supplied outputs, in addition to those provided by Oxs_Energy.
  void UpdateDerivedOutputs(const Oxs_SimState& state);
  // KL?? Is thi still used? Only in non-LLB mode? 
  // You mean: periodic-usage of kl_uniformexchange, yes?
  Oxs_ScalarOutput<Klm_UniformExchange> maxspinangle_output;
  Oxs_ScalarOutput<Klm_UniformExchange> stage_maxspinangle_output;
  Oxs_ScalarOutput<Klm_UniformExchange> run_maxspinangle_output;
  String MaxSpinAngleStateName() const {
    String dummy_name = InstanceName();
    dummy_name += ":Max Spin Angle";
    return dummy_name;
  }
  String StageMaxSpinAngleStateName() const {
    String dummy_name = InstanceName();
    dummy_name += ":Stage Max Spin Angle";
    return dummy_name;
  }
  String RunMaxSpinAngleStateName() const {
    String dummy_name = InstanceName();
    dummy_name += ":Run Max Spin Angle";
    return dummy_name;
  }
  // Do not control degree anymore (Max Spin Angle).
  // Instead we observe abolute magnetization change (units: A/m).
  // This outputs will be filled only in LLB mode KL??
  Oxs_ScalarOutput<Klm_UniformExchange> maxMdiff_output;
  Oxs_ScalarOutput<Klm_UniformExchange> stage_maxMdiff_output;
  Oxs_ScalarOutput<Klm_UniformExchange> run_maxMdiff_output;
  String MaxMdiffStateName() const {
    String dummy_name = InstanceName();
    dummy_name += ":Max m Diff";
    return dummy_name;
  }
  String StageMaxMdiffStateName() const {
    String dummy_name = InstanceName();
    dummy_name += ":Stage Max m Diff";
    return dummy_name;
  }
  String RunMaxMdiffStateName() const {
    String dummy_name = InstanceName();
    dummy_name += ":Run Max m Diff";
    return dummy_name;
  }
  mutable vector<OC_REAL8m> maxMdiffSq;
  // Maximum |M|/|Me|
  Oxs_ScalarOutput<Klm_UniformExchange> min_m_output;
  Oxs_ScalarOutput<Klm_UniformExchange> stage_min_m_output;
  Oxs_ScalarOutput<Klm_UniformExchange> run_min_m_output;
  Oxs_ScalarOutput<Klm_UniformExchange> max_m_output;  
  Oxs_ScalarOutput<Klm_UniformExchange> stage_max_m_output;  
  Oxs_ScalarOutput<Klm_UniformExchange> run_max_m_output;  
  String MinMStateName() const {
    String dummy_name = InstanceName();
    dummy_name += ":Min m";
    return dummy_name;
  }
  String StageMinMStateName() const {
    String dummy_name = InstanceName();
    dummy_name += ":Stage Min m";
    return dummy_name;
  }
  String RunMinMStateName() const {
    String dummy_name = InstanceName();
    dummy_name += ":Run Min m";
    return dummy_name;
  }
  String MaxMStateName() const {
    String dummy_name = InstanceName();
    dummy_name += ":Max m";
    return dummy_name;
  }
  String StageMaxMStateName() const {
    String dummy_name = InstanceName();
    dummy_name += ":Stage Max m";
    return dummy_name;
  }
  String RunMaxMStateName() const {
    String dummy_name = InstanceName();
    dummy_name += ":Run Max m";
    return dummy_name;
  }
  mutable vector<OC_REAL8m> min_m;
  mutable vector<OC_REAL8m> max_m;

protected:
  virtual void GetEnergy(const Oxs_SimState& state,
			 Oxs_EnergyData& oed) const {
    GetEnergyAlt(state,oed);
  }

  virtual void ComputeEnergy(const Oxs_SimState& state,
                             Oxs_ComputeEnergyData& oced) const {
    ComputeEnergyAlt(state,oced);
  }

  virtual void ComputeEnergyChunkInitialize
  (const Oxs_SimState& state,
   Oxs_ComputeEnergyDataThreaded& ocedt,
   Oc_AlignedVector<Oxs_ComputeEnergyDataThreadedAux>& thread_ocedtaux,
   int number_of_threads) const;

  virtual void ComputeEnergyChunkFinalize
  (const Oxs_SimState& state,
   const Oxs_ComputeEnergyDataThreaded& ocedt,
   Oc_AlignedVector<Oxs_ComputeEnergyDataThreadedAux>& thread_ocedtaux,
   int number_of_threads) const;

  virtual void ComputeEnergyChunk(const Oxs_SimState& state,
                                  Oxs_ComputeEnergyDataThreaded& ocedt,
                                  Oxs_ComputeEnergyDataThreadedAux& ocedtaux,
                                  OC_INDEX node_start,OC_INDEX node_stop,
                                  int threadnumber) const;

public:
  virtual const char* ClassName() const; // ClassName() is
  /// automatically generated by the OXS_EXT_REGISTER macro.
  Klm_UniformExchange(const char* name,     // Child instance id
		    Oxs_Director* newdtr, // App director
		    const char* argstr);  // MIF input block parameters
  virtual ~Klm_UniformExchange();
  virtual OC_BOOL Init();
};


#endif // _KLM_UNIFORMEXCHANGE
