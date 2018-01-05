/* FILE: kl_timedriver.cc            -*-Mode: c++-*-
 *
 * See appropriate h-file.
 *
 */

#include <string>

#include "nb.h"
#include "kl_timedriver.h"      // KL(m)
#include "scalarfield.h"        // KL(m)
#include "director.h"
#include "simstate.h"
#include "kl_timeevolvervarms.h"// KL(m)
#include "key.h"
#include "energy.h"		// Needed to make MSVC++ 5 happy

// KL(m)
#define KL_DEBUG 1

OC_USE_STRING;

// Oxs_Ext registration support
OXS_EXT_REGISTER(Klm_TimeDriver);

/* End includes */

// Constructor
Klm_TimeDriver::Klm_TimeDriver(
  const char* name,     // Child instance id
  Oxs_Director* newdtr, // App director
  const char* argstr)   // MIF input block parameters
  : Oxs_Driver(name,newdtr,argstr), max_dM_dt_obj_ptr(NULL)
{
  // Process arguments
  Oxs_OwnedPointer<Oxs_ScalarField> Msinit;

  OXS_GET_INIT_EXT_OBJECT("Ms_T0",Oxs_ScalarField,Msinit);
  // Fill Ms_T0 array, and verify that Ms_T0 is non-negative.
  Msinit->FillMeshValue(mesh_obj.GetPtr(), Ms_T0);
  for(OC_INDEX icell=0;icell<mesh_obj->Size();icell++)
    if(Ms_T0[icell]<0.0) {
      char buf[1024];
      Oc_Snprintf(buf,sizeof(buf),
		  "Negative Ms_T0 value (%g) detected at mesh index %u.",
		  static_cast<double>(Ms_T0[icell]),icell);
      throw Oxs_ExtError(this,String(buf));
    }

  OXS_GET_INIT_EXT_OBJECT("Ms_initial",Oxs_ScalarField,Msinit);
  // Fill Ms_inital array, and verify that 0 <= Ms_initial <= Ms_T0
  Msinit->FillMeshValue(mesh_obj.GetPtr(), Ms_initial);
  for(OC_INDEX icell=0;icell<mesh_obj->Size();icell++)
    if(Ms_initial[icell]<0.0 || Ms_initial[icell]>Ms_T0[icell]) {
      char buf[1024];
      Oc_Snprintf(buf,sizeof(buf),
		  "Value Ms_inital (%g) detected at mesh index %u"
		  " is either negative, or greater than Ms_T0 (%g)",
		  static_cast<double>(Ms_initial[icell]),icell,
		  static_cast<double>(Ms_T0[icell]));
      throw Oxs_ExtError(this,String(buf));
    }

  OXS_GET_INIT_EXT_OBJECT("evolver",Klm_TimeEvolverVarMs,evolver_obj);
  evolver_key.Set(evolver_obj.GetPtr());
  // Dependency lock on Klm_TimeEvolverVarMs object is
  // held until *this is destroyed.

  if(!HasInitValue("stopping_dM_dt")) {
    stopping_dM_dt.push_back(0.0); // Default is no control
  } else {
    GetGroupedRealListInitValue("stopping_dM_dt",stopping_dM_dt);
  }

  if(!HasInitValue("stopping_time")) {
    stopping_time.push_back(0.0); // Default is no control
  } else {
    GetGroupedRealListInitValue("stopping_time",stopping_time);
  }
  
  VerifyAllInitArgsUsed();

  last_timestep_output.Setup(
           this,InstanceName(),"Last time step","s",0,
	   &Klm_TimeDriver::Fill__last_timestep_output);
  simulation_time_output.Setup(
	   this,InstanceName(),"Simulation time","s",0,
	   &Klm_TimeDriver::Fill__simulation_time_output);

  last_timestep_output.Register(director,0);
  simulation_time_output.Register(director,0);

  // Reserve space for initial state (see GetInitialState() below)
  director->ReserveSimulationStateRequest(1);
}

Oxs_ConstKey<Oxs_SimState>
Klm_TimeDriver::GetInitialState() const
{
  Oxs_Key<Oxs_SimState> initial_state;
  director->GetNewSimulationState(initial_state);
  Oxs_SimState& istate = initial_state.GetWriteReference();
  SetStartValues(istate);

  // KL(m)
  // I assume NO RESTART situation.
  /// Restarting Ms-Var case has anyhow to be programmed separately.
  // Below I will overwrite one small thing coded in
  /// Oxs_Driver::SetStartValues copied below:
  //     istate.Ms.SetAsNonOwner(&Ms);
  //     istate.Ms_inverse.SetAsNonOwner(&Ms_inverse);
  // istate has no owned pointers. We have first to allocate them.
#if KL_DEBUG
  if( istate.Ms.IsOwner() ) {
    String msg="PROGRAMMING ERROR:"
      " istate.Ms.IsOwner() in Klm_TimeDriver::GetInitialState\n";
    throw Oxs_ExtError(this,msg);
  }
#endif // KL_DEBUG

  const Oxs_MeshValue<OC_REAL8m>* Ms_new_ptr;
  const Oxs_MeshValue<OC_REAL8m>* Ms_new_inv_ptr;
  Ms_new_ptr     = new const Oxs_MeshValue<OC_REAL8m>; // Allocate
  Ms_new_inv_ptr = new const Oxs_MeshValue<OC_REAL8m>;
  // We need here non-const pointer
// Old version: setting it by copying from Ms tables
//  *(const_cast<Oxs_MeshValue<OC_REAL8m>*>(Ms_new_ptr)) 
//                                    = *(Oxs_Driver::GetPtr_Ms());
//  *(const_cast<Oxs_MeshValue<OC_REAL8m>*>(Ms_new_inv_ptr))
//                                    = *(Oxs_Driver::GetPtr_Ms_inverse());
  // New version: setting it by copying from Ms_inital tables
  *(const_cast<Oxs_MeshValue<OC_REAL8m>*>(Ms_new_ptr))
                                    = Ms_initial;
  *(const_cast<Oxs_MeshValue<OC_REAL8m>*>(Ms_new_inv_ptr))
                                    = Ms_initial;
  // The Ms_new_inv_ptr table was copied mainly to allocate its space.
  /// Now, it has to be filled on a cell-by-cell basis.
  for(OC_INDEX icell=0;icell<mesh_obj->Size();icell++) {
    if(Ms_initial[icell]==0.0) {
      // Special case handling
      (*(const_cast<Oxs_MeshValue<OC_REAL8m>*>(Ms_new_inv_ptr)))[icell]=0.0;
    } else {
      (*(const_cast<Oxs_MeshValue<OC_REAL8m>*>(Ms_new_inv_ptr)))[icell]=
                                         1.0/Ms_initial[icell];
    }
  }

  // Set ownership for _these_ pointers
  istate.Ms.SetAsOwner(Ms_new_ptr);
  istate.Ms_inverse.SetAsOwner(Ms_new_inv_ptr);
  // We do similar thing with Ms_T0, Me_T vel Ms
  istate.Set_reference_Ms_Ptrs(&Ms_T0, this->GetPtr_Ms());
  // KL(m) end

  initial_state.GetReadReference();  // Release write lock.
  /// The read lock will be automatically released when the
  /// key "initial_state" is destroyed.
  return initial_state;
}

OC_BOOL Klm_TimeDriver::Init()
{ 
  Oxs_Driver::Init();  // Run init routine in parent.
  /// This will call Klm_TimeDriver::GetInitialState().

  // Get pointer to output object providing max dm/dt data
  const Klm_TimeEvolverVarMs* evolver = evolver_key.GetPtr();
  if(evolver==NULL) {
    throw Oxs_ExtError(this,"PROGRAMMING ERROR: No evolver found?");
  }
  
// KL(m)!
 
  String output_name = String(evolver->InstanceName());
  output_name += String(":Max dM/dt");
  max_dM_dt_obj_ptr
    =  director->FindOutputObjectExact(output_name.c_str());
  if(max_dM_dt_obj_ptr==NULL) {
    throw Oxs_ExtError(this,"Unable to identify unique"
                         " Max dM/dt output object");
  }

  return 1;
}

Klm_TimeDriver::~Klm_TimeDriver()
{}

void Klm_TimeDriver::StageRequestCount
(unsigned int& min,
 unsigned int& max) const
{ // Number of stages wanted by driver

  Oxs_Driver::StageRequestCount(min,max);

  unsigned int count = static_cast<OC_UINT4m>(stopping_dM_dt.size());
  if(count>min) min=count;
  if(count>1 && count<max) max=count;
  // Treat length 1 lists as imposing no upper constraint.

  count =  static_cast<OC_UINT4m>(stopping_time.size());
  if(count>min) min=count;
  if(count>1 && count<max) max=count;
  // Treat length 1 lists as imposing no upper constraint.
}

OC_BOOL
Klm_TimeDriver::ChildIsStageDone(const Oxs_SimState& state) const
{
  OC_UINT4m stage_index = state.stage_number;

  // Stage time check
  OC_REAL8m stop_time=0.;
  if(stage_index >= stopping_time.size()) {
    stop_time = stopping_time[stopping_time.size()-1];
  } else {
    stop_time = stopping_time[stage_index];
  }
  if(stop_time>0.0
     && stop_time-state.stage_elapsed_time<=stop_time*OC_REAL8_EPSILON*2) {
    return 1; // Stage done
  }

  // dM_dt check
  Tcl_Interp* mif_interp = director->GetMifInterp();
  if(max_dM_dt_obj_ptr==NULL ||
     max_dM_dt_obj_ptr->Output(&state,mif_interp,0,NULL) != TCL_OK) {
    String msg=String("Unable to obtain Max dM/dt output: ");
    if(max_dM_dt_obj_ptr==NULL) {
      msg += String("PROGRAMMING ERROR: max_dM_dt_obj_ptr not set."
		    " Driver Init() probably not called.");
    } else {
      msg += String(Tcl_GetStringResult(mif_interp));
    }
    throw Oxs_ExtError(this,msg.c_str());
  }
  OC_BOOL err;
  OC_REAL8m max_dM_dt = Nb_Atof(Tcl_GetStringResult(mif_interp),err);
  if(err) {
    String msg=String("Error detected in StageDone method"
		      " --- Invalid Max dM/dt output: ");
    msg += String(Tcl_GetStringResult(mif_interp));
    throw Oxs_ExtError(this,msg.c_str());
  }
  OC_REAL8m stop_dM_dt=0.;
  if(stage_index >= stopping_dM_dt.size()) {
    stop_dM_dt = stopping_dM_dt[stopping_dM_dt.size()-1];
  } else {
    stop_dM_dt = stopping_dM_dt[stage_index];
  }
  if(stop_dM_dt>0.0 && max_dM_dt <= stop_dM_dt) {
    return 1; // Stage done
  }

  // If control gets here, then stage not done
  return 0;
}

OC_BOOL
Klm_TimeDriver::ChildIsRunDone(const Oxs_SimState& /* state */) const
{
  // No child-specific checks at this time...
  return 0; // Run not done
}

void Klm_TimeDriver::FillStateSupplemental(Oxs_SimState& work_state) const
{
  OC_REAL8m work_step = work_state.last_timestep;
  OC_REAL8m base_time = work_state.stage_elapsed_time - work_step;

  // Insure that step does not go past stage stopping time
  OC_UINT4m stop_index = work_state.stage_number;
  OC_REAL8m stop_value=0.0;
  if(stop_index >= stopping_time.size()) {
    stop_value = stopping_time[stopping_time.size()-1];
  } else {
    stop_value = stopping_time[stop_index];
  }
  if(stop_value>0.0) {
    OC_REAL8m timediff = stop_value-work_state.stage_elapsed_time;
    if(timediff<=0) { // Over step
      // In the degenerate case where dm_dt=0, work_step will be
      // large (==1) and work_state.stage_elapsed_time will also
      // be large.  In that case, timediff will be numerically
      // poor because stop_value << work_state.stage_elapsed_time.
      // Check for this, and adjust sums accordingly.
      if(work_step>stop_value) { // Degenerate case
        work_step -= work_state.stage_elapsed_time;
        work_step += stop_value;
      } else {                   // Normal case
        work_step += timediff;
      }
      if(work_step<=0.0) work_step = stop_value*OC_REAL8_EPSILON; // Safety
      work_state.last_timestep = work_step;
      work_state.stage_elapsed_time = stop_value;
    } else if(timediff < 2*stop_value*OC_REAL8_EPSILON) {
      // Under step, but close enough for government work
      work_state.last_timestep += timediff;
      work_state.stage_elapsed_time = stop_value;
    } else if(0.25*work_step>timediff) {
      // Getting close to stage boundary.  Foreshorten.
      OC_REAL8m tempstep = (3*work_step+timediff)*0.25;
      work_state.last_timestep = tempstep;
      work_state.stage_elapsed_time = base_time+tempstep;
    }
  }
}

OC_BOOL
Klm_TimeDriver::Step
(Oxs_ConstKey<Oxs_SimState> base_state,
 const Oxs_DriverStepInfo& stepinfo,
 Oxs_ConstKey<Oxs_SimState>& next_state)
{ // Returns true if step was successful, false if
  // unable to step as requested.

  // Put write lock on evolver in order to get a non-const
  // pointer.  Use a temporary variable, temp_key, so
  // write lock is automatically removed when temp_key
  // is destroyed.
  Oxs_Key<Klm_TimeEvolverVarMs> temp_key = evolver_key;
  Klm_TimeEvolverVarMs& evolver = temp_key.GetWriteReference();
  return evolver.TryStep(this,base_state,stepinfo,next_state);
}

OC_BOOL
Klm_TimeDriver::InitNewStage
(Oxs_ConstKey<Oxs_SimState> state,
 Oxs_ConstKey<Oxs_SimState> prevstate)
{
  // Put write lock on evolver in order to get a non-const
  // pointer.  Use a temporary variable, temp_key, so
  // write lock is automatically removed when temp_key
  // is destroyed.
  Oxs_Key<Klm_TimeEvolverVarMs> temp_key = evolver_key;
  Klm_TimeEvolverVarMs& evolver = temp_key.GetWriteReference();
  return evolver.InitNewStage(this,state,prevstate);
}


////////////////////////////////////////////////////////////////////////
// State-based outputs, maintained by the driver.  These are
// conceptually public, but are specified private to force
// clients to use the output_map interface in Oxs_Director.

#define OSO_FUNC(NAME) \
void Klm_TimeDriver::Fill__##NAME##_output(const Oxs_SimState& state) \
{ NAME##_output.cache.state_id=state.Id(); \
  NAME##_output.cache.value=state.NAME; }

OSO_FUNC(last_timestep)

void
Klm_TimeDriver::Fill__simulation_time_output(const Oxs_SimState& state)
{
  simulation_time_output.cache.state_id = state.Id();
  simulation_time_output.cache.value =
    state.stage_start_time + state.stage_elapsed_time;
}
