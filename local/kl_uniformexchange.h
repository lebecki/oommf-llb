/* FILE: kl_uniformexchange.h            -*-Mode: c++-*-
 *
 * Uniform 6 neighbor exchange energy on rectangular mesh,
 * derived from Oxs_Energy class.
 * Based on file: uniformexchange.h
 * Extensions done by Kristof Lebecki, Universität Konstanz (Feb 2010):
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

class Klm_UniformExchange
  : public Oxs_ChunkEnergy, public Oxs_EnergyPreconditionerSupport {   
private:
  // KL(m) Actually, we support only "A_TYPE"
  // But I leave all these lines to keep the difference in code as small as possible
  enum ExchangeCoefType {
    A_UNKNOWN, A_TYPE, LEX_TYPE
  }  excoeftype;
  OC_REAL8m A;
  OC_REAL8m lex;

  enum ExchangeKernel { NGBR_UNKNOWN,
// KL(m) Based on NGBR_6_MIRROR, i.e. CalcEnergy6NgbrMirror
// (eventually on CalcEnergy6NgbrMirror_lex KL?)
			NGBR_6_LLB      // KL(m) LLB: influence of temperature
                      } kernel;
  /// Exchange formulation to use.  "unknown" is invalid; it
  /// is defined for error detection.
  /// NOTE: "kernel" is set inside the constructor, and should
  ///  be fixed thereafter.
  
  // Periodic boundaries
  mutable int xperiodic;
  mutable int yperiodic;
  mutable int zperiodic;

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

  // Support for threaded calculations
  mutable vector<OC_REAL8m> maxdot; // KL(m): here it is used
  // in a form of maxMdiffSq
  //mutable vector<OC_REAL8m> maxMdiffSq; <- my old naming
  mutable vector<OC_REAL8m> min_m;      // KL(m)
  mutable vector<OC_REAL8m> max_m;      // KL(m)

  // Calculation routines for each of the
  // aforementioned energy formulations.
  //KL(m)
  void CalcEnergy6NgbrMirror_LLB
  (const Oxs_MeshValue<ThreeVector>& spin,
   const Oxs_MeshValue<OC_REAL8m>& Ms_inverse,
   const Oxs_MeshValue<OC_REAL8m>& Ms, // KL(m) we need it as well
   const Oxs_MeshValue<OC_REAL8m>& Me_T, // KL(m)
   const Oxs_CommonRectangularMesh* mesh,
   Oxs_ComputeEnergyDataThreaded& ocedt,
   Oxs_ComputeEnergyDataThreadedAux& ocedtaux,
   OC_INDEX node_start,OC_INDEX node_stop,
   int threadnumber) const;

  // Supplied outputs, in addition to those provided by Oxs_Energy.
  void UpdateDerivedOutputs(const Oxs_SimState& state);
  // KL(m)
  // Do not control degree anymore (Max Spin Angle).
  // Instead we observe abolute magnetization change (units: A/m).
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
   Oxs_ComputeEnergyDataThreaded& ocedt,
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

  // Optional interface for conjugate-gradient evolver.
  virtual int IncrementPreconditioner(PreconditionerData& pcd);
};


#endif // _KLM_UNIFORMEXCHANGE
