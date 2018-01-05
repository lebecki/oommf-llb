/* FILE: kl_llbrungekuttaevolve.h                 -*-Mode: c++-*-
 *
 * This is an implementation of Landau-Lifshitz-Bloch dynamical
 * equation taking into account temperatures in the approximation
 * of the mean field - see Refs:
 *   [1] D. A. Garanin, 
 *       "Fokker-Planck and Landau-Lifshitz-Bloch equations for classical ferromagnets", 
 *       PHYSICAL REVIEW B 55(5), 3050 (1997).
 *   [2] N. Kazantseva, D. Hinzke, U. Nowak, R. W. Chantrell, U. Atxitia, 
 *       and O. Chubykalo-Fesenko, 
 *       "Towards multiscale modeling of magnetic materials: Simulations of FePt", 
 *       PHYSICAL REVIEW B 77, 184428 (2008).
 *
 * Important usage note.
 * The class Klm_LLB_RKEvolve must be used together with thermal correction
 * to the interactions, via usage of class:
 *   Klm_LLBextension
 * Additionally, if you use include interactions, you should compute them without
 * the assumption of |m|=1. Thus, taking the class Oxs_UniformExchange as an example,
 * you should specify kernel parameter value:
 *   kernel "26ngbr"
 *
 * Programming notes
 * - The |m|!=1 normalization error, norm_error, is not filled anymore.
 *   As we forget about the unity of |m|, this variable is constantly set to zero.
 *   Any calls to method ThreeVector::MakeUnit() are also removed.
 *
 * This file is a slightly modified version of rungekuttaevolve.h/.cc file.
 * I simply took it and made few modifications, to implement LLB instead of LL/LLG.
 * So, great thanks to the author of this original file (MD, I suppose).
 * Kristof Lebecki, 2010, Konstanz
 * KL(m) <- so are marked my important changes (changes of class names are not marked)
 *
 * Concrete evolver class, using Runge-Kutta steps
 *
 */

#ifndef _KLM_LLBRUNGEKUTTAEVOLVE
#define _KLM_LLBRUNGEKUTTAEVOLVE

#include <vector>

#include "kl_timeevolvervarms.h"
#include "key.h"
#include "mesh.h"
#include "meshvalue.h"
#include "scalarfield.h"
#include "output.h"

/* End includes */

#if REPORT_TIME
# ifndef REPORT_TIME_RKDEVEL
#  define REPORT_TIME_RKDEVEL 1
# endif
#endif

#if OOMMF_THREADS
  class _Klm_LLB_RKEvolve_RKFBase54_ThreadA; // Forward references
  class _Klm_LLB_RKEvolve_RKFBase54_ThreadB;
  class _Klm_LLB_RKEvolve_RKFBase54_ThreadC;
  class _Klm_LLB_RKEvolve_RKFBase54_ThreadD;
#endif

class Klm_LLB_RKEvolve:public Klm_TimeEvolverVarMs {
private:
#if REPORT_TIME_RKDEVEL
  mutable vector<Nb_StopWatch> timer;
  struct TimerCounts {
  public:
    OC_INT4m pass_count;
    unsigned long bytes;
    String name;
    TimerCounts() : pass_count(0), bytes(0) {}
    void Reset() { pass_count = 0; bytes = 0; }
  };
  mutable vector<TimerCounts> timer_counts;
#endif

  mutable OC_UINT4m mesh_id;     // Used by gamma and alpha meshvalues to
  void UpdateMeshArrays(const Oxs_Mesh*);   /// track changes in mesh.

  // KL(m). Additional input parameter, relative temparature, t=T/Tc
  OC_REAL8m relative_temperature;

  // Base step size control parameters
  OC_REAL8m min_timestep;           // Seconds
  OC_REAL8m max_timestep;           // Seconds

  const OC_REAL8m max_step_decrease;        // Safety size adjusment
  const OC_REAL8m max_step_increase_limit;  // bounds.
  const OC_REAL8m max_step_increase_adj_ratio;
  OC_REAL8m max_step_increase;
  /// NOTE: These bounds do not include step_headroom, which
  /// is applied at the end.

  // Error-based step size control parameters.  Each may be disabled
  // by setting to -1.  There is an additional step size control that
  // insures that energy is monotonically non-increasing (up to
  // estimated rounding error).
  OC_REAL8m allowed_error_rate;  // Step size is adjusted so
  /// that the estimated maximum error (across all spins) divided
  /// by the step size is smaller than this value.  The units
  /// internally are radians per second, converted from the value
  /// specified in the input MIF file, which is in deg/sec.
  /// KL(m): Units changed! Now MIF: (A/m*s)
  ///        internally and externally

  OC_REAL8m allowed_absolute_step_error; // Similar to allowed_error_rate,
  /// but without the step size adjustment.  Internal units are
  /// radians; MIF input units are degrees.
  /// KL(m): Units changed! Now MIF: (A/m)
  ///        internally and externally

  OC_REAL8m allowed_relative_step_error;
  // Step size is adjusted so that
  /// the estimated maximum error (across all spins) divided by
  /// [maximum dm/dt (across all spins) * step size] is smaller than
  /// this value.  This value is non-dimensional, representing the
  /// allowed relative (proportional) error, presumably in (0,1).

  OC_REAL8m expected_energy_precision; // Expected relative energy
  /// precision.

  // KL(m)
  OC_REAL8m expected_m_max_precision; // Cases |M|>|Ms_T0|
  OC_REAL8m expected_m_min; // Cases |M| close to zero

  OC_REAL8m reject_goal,reject_ratio;
  OC_REAL8m min_step_headroom,max_step_headroom;
  OC_REAL8m step_headroom; // Safety margin used in step size adjustment

  // Spatially variable Landau-Lifschitz-Gilbert gyromagnetic ratio
  // and damping coefficients.
  OC_BOOL do_precess;  // If false, then do pure damping
  OC_BOOL allow_signed_gamma; // If false, then force gamma negative
  enum GammaStyle { GS_INVALID, GS_LL, GS_G }; // Landau-Lifshitz or Gilbert
  GammaStyle gamma_style;
  Oxs_OwnedPointer<Oxs_ScalarField> gamma_init;
  mutable Oxs_MeshValue<OC_REAL8m> gamma;

  Oxs_OwnedPointer<Oxs_ScalarField> alpha_init;
  mutable Oxs_MeshValue<OC_REAL8m> alpha;

  // The next timestep is based on the error from the last step.  If
  // there is no last step (either because this is the first step,
  // or because the last state handled by this routine is different
  // from the incoming current_state), then timestep is calculated
  // so that max_dM_dt * timestep = start_dM,  ( KL(m) change )
  // or timestep = start_dt,
  // whichever is smaller.  Either can be disabled by setting <0.
  OC_REAL8m start_dM; // KL(m) See KLnote1
  OC_REAL8m start_dt;

  // Stepsize control for first step of each stage after the first.
  // Choices are to use start conditions (start_dM and/or start_dt),
  // use continuation from end of previous stage, or to automatically
  // select between the two methods depending on whether or not the
  // energy appears to be continuous across the stage boundary.
  enum StageInitStepControl { SISC_INVALID, SISC_START_COND,
			      SISC_CONTINUOUS, SISC_AUTO };
  StageInitStepControl stage_init_step_control;

  // Data cached from last state
  OC_UINT4m energy_state_id;
  Oxs_MeshValue<OC_REAL8m> energy;
  OC_REAL8m next_timestep;

  // Outputs
  void UpdateDerivedOutputs(const Oxs_SimState& state,
                            const Oxs_SimState* prevstate);
  void UpdateDerivedOutputs(const Oxs_SimState& state) {
    UpdateDerivedOutputs(state,NULL);
  }
  Oxs_ScalarOutput<Klm_LLB_RKEvolve> max_dM_dt_output;  // See KLnote1
  Oxs_ScalarOutput<Klm_LLB_RKEvolve> dE_dt_output;
  Oxs_ScalarOutput<Klm_LLB_RKEvolve> delta_E_output;
  // These outputs are moved to exchange module.
//  Oxs_ScalarOutput<Klm_LLB_RKEvolve> min_m_output;      // KL(m)
//  Oxs_ScalarOutput<Klm_LLB_RKEvolve> max_m_output;      // KL(m)
  Oxs_VectorFieldOutput<Klm_LLB_RKEvolve> dM_dt_output; // See KLnote1
  Oxs_VectorFieldOutput<Klm_LLB_RKEvolve> mxH_output;

  // Scratch space
  Oxs_MeshValue<OC_REAL8m> temp_energy;
  Oxs_MeshValue<ThreeVector> vtmpA;
  Oxs_MeshValue<ThreeVector> vtmpB;
  Oxs_MeshValue<ThreeVector> vtmpC;
  Oxs_MeshValue<ThreeVector> vtmpD;
  Oxs_MeshValue<ThreeVector> vtmpE; /**/
  // KL(m) space to hold field H
  // Above listed vtmp-tables are used simultaneously for two purposes:
  // - as Calculate_dM_dt-input parameter holding mxH values,
  // - as Calculate_dM_dt-output parameter giving dM_dt value.
  //
  // KL(m)?
  // kltmp-tables are used as Calculate_dm_dt-input parameter 
  // holding H values. This is most probably waste of space 
  // - only one table would be enough for that.
  // This input parameter is calculated by (earlier) 
  // GetEnergyDensity function.
  // KL(m)?
  Oxs_MeshValue<ThreeVector> kltmpA;
  Oxs_MeshValue<ThreeVector> kltmpB;
  Oxs_MeshValue<ThreeVector> kltmpC;
  Oxs_MeshValue<ThreeVector> kltmpD;
  Oxs_MeshValue<ThreeVector> kltmpX;

  // Utility functions
  void CheckCache(const Oxs_SimState& cstate);

  void AdjustState(OC_REAL8m hstep,
		   OC_REAL8m mstep,
		   const Oxs_SimState& old_state,
		   const Oxs_MeshValue<ThreeVector>& dM_dt,
		   Oxs_SimState& new_state,
		   OC_REAL8m& norm_error) const;
  // Export new state has time index from old_state + h,
  // and spins from old state + mstep*dM_dt.

#if OOMMF_THREADS
  void AdjustState
  (OC_REAL8m hstep,
   OC_REAL8m mstep,
   const Oxs_SimState& old_state,
   vector<OC_REAL8m>& b,
   vector<const Oxs_MeshValue<ThreeVector>*>& A,
   Oxs_MeshValue<ThreeVector>& kn, // target = b*source
   Oxs_SimState& new_state,
   OC_REAL8m& norm_error) const;
!!! m->M
  // Same as preceding AdjustState, except that dm_dt is an export,
  // which is filled via kn[i] = \sum_j b[j]*A[j][i] 

  void AdjustState
  (OC_REAL8m hstep,
   OC_REAL8m mstep,
   const Oxs_SimState& old_state,
   const Oxs_MeshValue<ThreeVector>& partial_sum,
   vector<OC_REAL8m>& b,
   vector<const Oxs_MeshValue<ThreeVector>*>& A,
   Oxs_SimState& new_state,
   OC_REAL8m& norm_error) const;
  // Same as preceding AdjustState, except that kn is not
  // exported, but is computed for internal use and computed
  // as
  //     kn[i] = partial_sum[i] + \sum_j b[j]*A[j][i] 

#endif // OOMMF_THREADS


  void UpdateTimeFields(const Oxs_SimState& cstate,
			Oxs_SimState& nstate,
			OC_REAL8m stepsize) const;

  void NegotiateTimeStep(const Klm_TimeDriver* driver,
			 const Oxs_SimState&  cstate,
			 Oxs_SimState& nstate,
			 OC_REAL8m stepsize,
			 OC_BOOL use_start_cond,
			 OC_BOOL& forcestep,
			 OC_BOOL& driver_set_step) const;

  OC_BOOL CheckError(OC_REAL8m global_error_order,OC_REAL8m error,
		  OC_REAL8m stepsize,OC_REAL8m reference_stepsize,
		  OC_REAL8m max_dM_dt,OC_REAL8m& new_stepsize);
  /// Returns 1 if step is good, 0 if error is too large.
  /// Export new_stepsize is set to suggested stepsize
  /// for next step.

  OC_REAL8m MaxDiff_mul(const Oxs_MeshValue<ThreeVector>& vecA,
                     const Oxs_MeshValue<OC_REAL8m>&      vecA_mul,
		     const Oxs_MeshValue<ThreeVector>& vecB,
                     const Oxs_MeshValue<OC_REAL8m>&      vecB_mul);
  /// Returns maximum difference between vectors multiplied by factors,
  /// in corresponding positions in two vector fields.

  void AdjustStepHeadroom(OC_INT4m step_reject);
  /// step_reject should be 0 or 1, reflecting whether the current
  /// step was rejected or not.  This routine updates reject_ratio
  /// and adjusts step_headroom appropriately.

  void ComputeEnergyChange(const Oxs_Mesh* mesh,
                           const Oxs_MeshValue<OC_REAL8m>& current_energy,
                           const Oxs_MeshValue<OC_REAL8m>& candidate_energy,
                           OC_REAL8m& dE,OC_REAL8m& var_dE,OC_REAL8m& total_E);
  /// Computes cellwise difference between energies, and variance.
  /// Export total_E is "current" energy (used for stepsize control).


  // Stepper routines:  If routine needs to compute the energy
  // at the new (final) state, then it should store the final
  // energy results in temp_energy, mxH in mxH_output.cache,
  // and dM_dt into the vtmpA scratch array, fill
  // the "Timestep lower bound", "Max dM/dt", "dE/dt", and
  // "pE/pt" derived data fields in nstate, and set the export
  // value new_energy_and_dmdt_computed true.  Otherwise the export
  // value should be set false, and the client routine is responsible
  // for obtaining these values as necessary.  (If possible, it is
  // better to let the client compute these values, because the
  // client may be able to defer computation until it has decided
  // whether or not to keep the step.)
  // KL(m)! should I read/understand it?

  // One would like to declare the step functions and pointer
  // to same via typedef's, but the MS VC++ 6.0 (& others?)
  // compiler doesn't handle member function typedef's properly---
  // it produces __cdecl linkage rather than instance member
  // linkage.  Typedef's on pointers to member functions work
  // okay, just not typedef's on member functions themselves.
  // So, instead we use a #define, which is ugly but portable.
#define RKStepFuncSig(NAME) \
  void NAME (                                            \
     OC_REAL8m stepsize,                                    \
     Oxs_ConstKey<Oxs_SimState> current_state,           \
     const Oxs_MeshValue<ThreeVector>& current_dm_dt,    \
     Oxs_Key<Oxs_SimState>& next_state,                  \
     OC_REAL8m& error_estimate,                             \
     OC_REAL8m& global_error_order,                         \
     OC_REAL8m& norm_error,                                 \
     OC_BOOL& new_energy_and_dmdt_computed)

  // Functions that calculate a single RK step
  RKStepFuncSig(TakeRungeKuttaFehlbergStep54);
  RKStepFuncSig(TakeRungeKuttaFehlbergStep54M);
  RKStepFuncSig(TakeRungeKuttaFehlbergStep54S);

  // Pointer set at runtime during instance initialization
  // to one of the above functions single RK step functions.
  RKStepFuncSig((Klm_LLB_RKEvolve::* rkstep_ptr));

  // Utility code used by the TakeRungeKuttaFehlbergStep54* routines.
  enum RKF_SubType { RKF_INVALID, RK547FC, RK547FM, RK547FS };
  void RungeKuttaFehlbergBase54(RKF_SubType method,
			   OC_REAL8m stepsize,
			   Oxs_ConstKey<Oxs_SimState> current_state,
			   const Oxs_MeshValue<ThreeVector>& current_dM_dt,
			   Oxs_Key<Oxs_SimState>& next_state,
			   OC_REAL8m& error_estimate,
			   OC_REAL8m& global_error_order,
			   OC_REAL8m& norm_error,
			   OC_BOOL& new_energy_and_dmdt_computed);

//  void Calculate_dm_dt
  void Calculate_dM_dt // See cc-file KLnote1
  (const Oxs_SimState& state_,
   const Oxs_MeshValue<ThreeVector>& mxH_,
   const Oxs_MeshValue<ThreeVector>& H_, // KL(m). Additionally needed by LLB
   OC_REAL8m pE_pt_,
   Oxs_MeshValue<ThreeVector>& dM_dt_,
   OC_REAL8m& max_dM_dt_,OC_REAL8m& dE_dt_,OC_REAL8m& min_timestep_);
  /// Imports: state_, mxH_, pE_pt
  /// Exports: dm_dt_, max_dm_dt_, dE_dt_, min_timestep_

#if OOMMF_THREADS
  // Thread-friendly version of Calculate_dm_dt breaks into
  // three parts: the Initialize and Finalize members here
  // and the Cmd member of _Klm_LLB_RKEvolve_RKFBase54_ThreadA
  // (which is defined in rungekuttaevolve.cc).
  void Initialize_Threaded_Calculate_dm_dt
  (const Oxs_SimState& state, // Import
   const Oxs_MeshValue<ThreeVector>& mxH, // Import
!!! m->M
   Oxs_MeshValue<ThreeVector>& dm_dt,     // Import ptr to export data
   vector<_Klm_LLB_RKEvolve_RKFBase54_ThreadA>& thread_data); // Export

  void Finalize_Threaded_Calculate_dm_dt
  (const vector<_Klm_LLB_RKEvolve_RKFBase54_ThreadA>& thread_data, // Import
   OC_REAL8m pE_pt,         // Import
   OC_REAL8m& max_dm_dt,    // Export
   OC_REAL8m& dE_dt,        // Export
   OC_REAL8m& min_timestep_export); // Export
#endif // OOMMF_THREADS

#if OOMMF_THREADS
  // Thread-friendly version of Calculate_dm_dt breaks into
  // three parts: the Initialize and Finalize members here
  // and the Cmd member of _Klm_LLB_RKEvolve_RKFBase54_ThreadB
  // (which is defined in rungekuttaevolve.cc).
  void Initialize_Threaded_Calculate_dm_dt
  (const Oxs_SimState& base_state, // Import
   Oxs_SimState& work_state,       // Import and export
   const Oxs_MeshValue<ThreeVector>& mxH, // Import
!!! m->M
   Oxs_MeshValue<ThreeVector>& dm_dt,     // Export; may be same as mxH
   vector<_Klm_LLB_RKEvolve_RKFBase54_ThreadB>& thread_data, // Export
   OC_REAL8m mstep,
   OC_REAL8m b_dm_dt,
   const vector<OC_REAL8m>& b,
   vector<const Oxs_MeshValue<ThreeVector>*>& A,
   Oxs_MeshValue<ThreeVector>& kn); // Export kn = b*A
  /// Also sets spins in work_state to base_state.spin + mstep*(b*A)

  void Finalize_Threaded_Calculate_dm_dt
  (const Oxs_SimState& base_state, // Import
   const vector<_Klm_LLB_RKEvolve_RKFBase54_ThreadB>& thread_data, // Import
   OC_REAL8m pE_pt,             // Import
   OC_REAL8m hstep,             // Import
   Oxs_SimState& work_state, // Export
   OC_REAL8m& max_dm_dt,        // Export
   OC_REAL8m& dE_dt,            // Export
   OC_REAL8m& min_timestep_export,  // Export
   OC_REAL8m& norm_error);      // Export
#endif // OOMMF_THREADS

#if OOMMF_THREADS
  // Thread-friendly version of Calculate_dm_dt breaks into
  // three parts: the Initialize and Finalize members here
  // and the Cmd member of _Klm_LLB_RKEvolve_RKFBase54_ThreadC
  // (which is defined in rungekuttaevolve.cc).
  void Initialize_Threaded_Calculate_dm_dt
  (const Oxs_SimState& base_state, // Import
   Oxs_SimState& work_state,       // Import and export
   const Oxs_MeshValue<ThreeVector>& mxH, // Import
!!! m->M
   Oxs_MeshValue<ThreeVector>& dm_dt,     // Export; may be same as mxH
   vector<_Klm_LLB_RKEvolve_RKFBase54_ThreadC>& thread_data, // Export
   OC_REAL8m mstep,
   OC_REAL8m b1_dm_dt,
   OC_REAL8m b2_dm_dt,
   const vector<OC_REAL8m>& b1,
   const vector<OC_REAL8m>& b2,
   vector<const Oxs_MeshValue<ThreeVector>*>& A,
   Oxs_MeshValue<ThreeVector>& kn); // Export kn = b2*A
  /// Also sets spins in work_state to base_state.spin + mstep*(b1*A)

  void Finalize_Threaded_Calculate_dm_dt
  (const Oxs_SimState& base_state, // Import
   const vector<_Klm_LLB_RKEvolve_RKFBase54_ThreadC>& thread_data, // Import
   OC_REAL8m pE_pt,             // Import
   OC_REAL8m hstep,             // Import
   Oxs_SimState& work_state, // Export
   OC_REAL8m& max_dm_dt,        // Export
   OC_REAL8m& dE_dt,            // Export
   OC_REAL8m& min_timestep_export,  // Export
   OC_REAL8m& norm_error);      // Export
#endif // OOMMF_THREADS

#if OOMMF_THREADS
  // Thread-friendly version of Calculate_dm_dt breaks into
  // three parts: the Initialize and Finalize members here
  // and the Cmd member of _Klm_LLB_RKEvolve_RKFBase54_ThreadD
  // (which is defined in rungekuttaevolve.cc).
  void Initialize_Threaded_Calculate_dm_dt
  (const Oxs_SimState& work_state,        // Import
   const Oxs_MeshValue<ThreeVector>& mxH, // Import
!!! m->M
   Oxs_MeshValue<ThreeVector>& dm_dt,     // Export
   vector<_Klm_LLB_RKEvolve_RKFBase54_ThreadD>& thread_data, // Import and export
   OC_REAL8m dc7, // Import
   const Oxs_MeshValue<ThreeVector>& dD13456); // Import

  void Finalize_Threaded_Calculate_dm_dt
  (const vector<_Klm_LLB_RKEvolve_RKFBase54_ThreadD>& thread_data, // Import
   OC_REAL8m pE_pt,          // Import
   OC_REAL8m& max_dm_dt,      // Export
   OC_REAL8m& dE_dt,          // Export
   OC_REAL8m& min_timestep_export,  // Export
   OC_REAL8m& max_dD_magsq);  // Export
#endif // OOMMF_THREADS

  // Declare but leave undefined copy constructor and assignment operator
  Klm_LLB_RKEvolve(const Klm_LLB_RKEvolve&);
  Klm_LLB_RKEvolve& operator=(const Klm_LLB_RKEvolve&);

public:
  virtual const char* ClassName() const; // ClassName() is
  /// automatically generated by the OXS_EXT_REGISTER macro.
  virtual OC_BOOL Init();
  Klm_LLB_RKEvolve(const char* name,     // Child instance id
		       Oxs_Director* newdtr, // App director
		       const char* argstr);  // MIF input block parameters
  virtual ~Klm_LLB_RKEvolve();

  virtual OC_BOOL
  InitNewStage(const Klm_TimeDriver* driver,
               Oxs_ConstKey<Oxs_SimState> state,
               Oxs_ConstKey<Oxs_SimState> prevstate);

  virtual  OC_BOOL
  Step(const Klm_TimeDriver* driver,
       Oxs_ConstKey<Oxs_SimState> current_state,
       const Oxs_DriverStepInfo& step_info,
       Oxs_Key<Oxs_SimState>& next_state);
  // Returns true if step was successful, false if
  // unable to step as requested.
};

#endif // _KLM_LLBRUNGEKUTTAEVOLVE
