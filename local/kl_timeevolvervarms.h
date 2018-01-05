/* FILE: kl_timeevolvervarms.h                 -*-Mode: c++-*-
 *
 * This is an evolver class complimentary to Klm_TimeDriver.
 * It contains no additional functionality - it only calls/refers
 * to Klm_TimeDriver, instead of Oxs_TimeDriver.
 * Actually a more sphisticated design change is probably needed
 * - see fragment of Mike donahue mail below. As for now, however,
 * I stick to this simple solution.
 * 
 *
 * Maybe the evolver should be broken up into two classes.  One class holds
 * the code, which btw doesn't know anything about Oxs_TimeDriver.  The
 * other class is just an interface class with the two member functions
 * InitNewStage and Step.  One can argue on the best way to order the
 * hierarchy of these classes, but in practice, if we want to not break
 * existing extensions (e.g., Anv_SpinTEvolve), then the interface class
 * would need to be a child of the code class, and concrete children
 * inherit off of the interface class.  In this scenario, you would just
 * need to write a replacement interface class, and have your evolver
 * child(ren) inherit off of your interface class.
 * 
 * I think this then frees up your driver class to do whatever you want,
 * because the time evolver code class doesn't #include or define
 * Oxs_TimeDriver.
 * 
 * So, if there is code in the Oxs_TimeDriver class that you want to use,
 * and if existing children of Oxs_TimeEvolver (e.g., Oxs_RungeKuttaEvolve)
 * would be happy using your driver through the Oxs_TimeDriver interface,
 * then you can inherit off of Oxs_TimeDriver.  In this case I think we
 * should re-factor Oxs_TimeDriver slighly to make it a proper base class;
 * in particular, pull VerifyAllInitArgsUsed() out of the Oxs_TimeDriver
 * constructor and place it into a concrete child class, but that is a
 * detail we can look at later.
 * 
 * If existing children of Oxs_TimeEvolver can't use your driver, or if
 * there is no code in Oxs_TimeDriver that you want to use, then inherit
 * off of Oxs_Driver and write a class parallel to Oxs_TimeDriver.
 *
 */

#ifndef _KLM_TIMEEVOLVERVARMS
#define _KLM_TIMEEVOLVERVARMS

#include "evolver.h"
#include "output.h"

/* End includes */

class Klm_TimeDriver; // Forward references
struct Oxs_DriverStepInfo;

class Klm_TimeEvolverVarMs:public Oxs_Evolver {
private:

  OC_UINT4m energy_calc_count; // Number of times GetEnergyDensity
  /// has been called in current problem run.

  Oxs_MeshValue<OC_REAL8m> temp_energy;     // Scratch space used by
  Oxs_MeshValue<ThreeVector> temp_field; // GetEnergyDensity().

  // Outputs maintained by this interface layer.  These are conceptually
  // public, but are specified private to force clients to use the
  // output_map interface.
  void UpdateEnergyOutputs(const Oxs_SimState&);
  void FillEnergyCalcCountOutput(const Oxs_SimState&);
  Oxs_ScalarOutput<Klm_TimeEvolverVarMs> total_energy_output;
  Oxs_ScalarFieldOutput<Klm_TimeEvolverVarMs> total_energy_density_output;
  Oxs_VectorFieldOutput<Klm_TimeEvolverVarMs> total_field_output;
  Oxs_ScalarOutput<Klm_TimeEvolverVarMs> energy_calc_count_output;

  // Disable copy constructor and assignment operator by declaring
  // them without defining them.
  Klm_TimeEvolverVarMs(const Klm_TimeEvolverVarMs&);
  Klm_TimeEvolverVarMs& operator=(const Klm_TimeEvolverVarMs&);

protected:

#if REPORT_TIME
  mutable Nb_StopWatch steponlytime;
#endif

  Klm_TimeEvolverVarMs(const char* name,      // Child instance id
                 Oxs_Director* newdtr);  // App director
  Klm_TimeEvolverVarMs(const char* name,
                 Oxs_Director* newdtr,
                 const char* argstr);      // MIF block argument string

  virtual OC_BOOL Init();  // All children of Klm_TimeEvolverVarMs *must*
  /// call this function in their Init() routines.  The main purpose
  /// of this function is to initialize output variables.

  void GetEnergyDensity(const Oxs_SimState& state,
                        Oxs_MeshValue<OC_REAL8m>& energy,
                        Oxs_MeshValue<ThreeVector>* mxH_req,
                        Oxs_MeshValue<ThreeVector>* H_req,
                        OC_REAL8m& pE_pt,
                        OC_REAL8m& total_E);

  void GetEnergyDensity(const Oxs_SimState& state,
			Oxs_MeshValue<OC_REAL8m>& energy,
			Oxs_MeshValue<ThreeVector>* mxH_req,
			Oxs_MeshValue<ThreeVector>* H_req,
			OC_REAL8m& pE_pt) {
    // This interface for backwards compatibility
    OC_REAL8m dummy_E;
    GetEnergyDensity(state,energy,mxH_req,H_req,pE_pt,dummy_E);
  }

public:
  virtual ~Klm_TimeEvolverVarMs();

  virtual OC_BOOL
  InitNewStage(const Klm_TimeDriver* /* driver */,
               Oxs_ConstKey<Oxs_SimState> /* state */,
               Oxs_ConstKey<Oxs_SimState> /* prevstate */) { return 1; }
  /// Default implementation is a NOP.  Children may override.
  /// NOTE: prevstate may be "INVALID".  Children should check
  ///       before use.

  // There are two versions of the Step routine, one (Step) where
  // next_state is an Oxs_Key<Oxs_SimState>& import/export value, and
  // one (TryStep) where next_state is an Oxs_ConstKey<Oxs_SimState>&
  // export-only value.  In the former, the caller initializes
  // next_state with a fresh state by calling
  // director->GetNewSimulationState(), which leaves new_state with a
  // write lock.  The latter interface is more recent; with it the
  // Step code should make its own calls to
  // director->GetNewSimulationState() to get working states as
  // needed, and pass the "next" state back to the caller through the
  // next_state reference.  This is more flexible, and should be
  // preferred by new code.  In particular, it allows the
  // evolver::Step() routine to return a reference to a pre-existing
  // state held in an Oxs_ConstKey<Oxs_SimState>, without having to
  // invoke a const_cast<> operation.  Be aware, however, that this
  // newer interface involves a transfer of "ownership" of the state.
  // In the original Oxs_Key interface the caller created next_state,
  // passed it to Step(), received it back, and ultimately dispose of
  // it.  In the newer Oxs_ConstKey interface, the Step() routine
  // creates next_state and passes back to the caller.  But this is
  // what the Oxs_Key class is designed for, so it should work fine.
  //
  // In both cases, the return value is true if step was successful,
  // false if unable to step as requested.  Also, the evolver object
  // is responsible for calling driver->FillState() to fill next_state
  // as needed.
  //
  // As a migration aid, a default implementations are provided.  The
  // Oxs_Driver code calls the newer interface; child evolver classes
  // should  override exactly one of these.
  //
  // Note: We can't give the same name to both routines because of implicit
  // conversion of Oxs_Key<T> to Oxs_ConstKey<T>.
  //
  virtual OC_BOOL
  Step(const Klm_TimeDriver* /* driver */,
       Oxs_ConstKey<Oxs_SimState> /* current_state */,
       const Oxs_DriverStepInfo& /* step_info */,
       Oxs_Key<Oxs_SimState>& /* next_state */)
 {
    throw Oxs_ExtError(this,
          "Programming error: Implementation of"
          " Oxs_TimeEvolver::Step(const Klm_TimeDriver*,Oxs_ConstKey<Oxs_SimState>,Oxs_Key<Oxs_SimState>&)"
          " not provided.");
    return 0;
  }

  virtual OC_BOOL
  TryStep(const Klm_TimeDriver* driver,
       Oxs_ConstKey<Oxs_SimState> current_state,
       const Oxs_DriverStepInfo& step_info,
       Oxs_ConstKey<Oxs_SimState>& next_state) {
    // Default implementation that wraps older Step() call.  New code
    // should override this implementation and ignore old Step()
    // interface.
    Oxs_Key<Oxs_SimState> temp_state;
    director->GetNewSimulationState(temp_state);
    OC_BOOL result = Step(driver,current_state,step_info,temp_state);
    temp_state.GetReadReference();
    next_state = temp_state;
    next_state.GetReadReference();
    return result;
  }
};

#endif // _KLM_TIMEEVOLVERVARMS
