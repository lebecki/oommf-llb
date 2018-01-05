/* FILE: kl_anisotropy.cc            -*-Mode: c++-*-
 *
 * See appropriate h-file.
 *
 */

#include <limits>
#include <string>

#include "oc.h"
#include "nb.h"
#include "threevector.h"
#include "director.h"
#include "simstate.h"
#include "ext.h"
#include "key.h"
#include "mesh.h"
#include "meshvalue.h"
#include "uniformscalarfield.h"
#include "uniformvectorfield.h"
#include "rectangularmesh.h"  // For QUAD-style integration
#include "kl_anisotropy.h" // KL(m)
#include "kl_llb_util.h" // KL(m) 
#include "energy.h"		// Needed to make MSVC++ 5 happy

OC_USE_STRING;

// Oxs_Ext registration support
OXS_EXT_REGISTER(Klm_Anisotropy);

/* End includes */

// Constructor
Klm_Anisotropy::Klm_Anisotropy(
  const char* name,     // Child instance id
  Oxs_Director* newdtr, // App director
  const char* argstr)   // MIF input block parameters
  : Oxs_ChunkEnergy(name,newdtr,argstr), mesh_id(0),
    K1_is_uniform(0),axis_is_uniform(0),uniform_K1_value(0.0),
    uniform_chi_x_value(0.0), // KL(m)
    uniform_chi_y_value(0.0), // KL(m)
    uniform_chi_z_value(0.0), // KL(m)
    integration_method(UNKNOWN_INTEG)
{
  // Process arguments
  // KL(m)
  if(HasInitValue("K1")) // K1 is not a required parameter anymore
  { // Standard OOMMF behavior
    // Here K1 and axis are required
    OXS_GET_INIT_EXT_OBJECT("K1",Oxs_ScalarField,K1_init);
    Oxs_UniformScalarField* tmpK1ptr
      = dynamic_cast<Oxs_UniformScalarField*>(K1_init.GetPtr());
    if(tmpK1ptr) {
      // Initialization is via a uniform field; set up uniform
      // K1 variables.
      K1_is_uniform = 1;
      uniform_K1_value = tmpK1ptr->SoleValue();
    }
    OXS_GET_INIT_EXT_OBJECT("axis",Oxs_VectorField,axis_init);
    Oxs_UniformVectorField* tmpaxisptr
      = dynamic_cast<Oxs_UniformVectorField*>(axis_init.GetPtr());
    if(tmpaxisptr) {
      // Initialization is via a uniform field.  For convenience,
      // modify the size of the field components to norm 1, as
      // required for the axis specification.  This allows the
      // user to specify the axis direction as, for example, {1,1,1},
      // as opposed to {0.57735027,0.57735027,0.57735027}, or
      //
      //      Specify Oxs_UniformVectorField {
      //        norm 1 
      //        vector { 1 1 1 } 
      //    }
      // Also setup uniform axis variables
      tmpaxisptr->SetMag(1.0);
      axis_is_uniform = 1;
      uniform_axis_value = tmpaxisptr->SoleValue();
    }
    if(HasInitValue("chi_x") || HasInitValue("chi_y") || HasInitValue("chi_z")) {
      // Either K1 of susceptibility, not both
      String msg=String("\nYou must either specify susceptibility"
          " (via \"chi_x\", \"chi_y\", or \"chi_z\"),"
          " or anisotropy constant K1. Not both.");
      throw Oxs_ExtError(this,msg.c_str());
    }
  }
  else
  { // LLB-behavior
    K1_init.SetAsNonOwner(NULL);
    axis_init.SetAsNonOwner(NULL);
    // In such a case at least one chi_[xyz] must be supplied
    // Present implementation supports only uniform chi_perpendicular
    uniform_chi_x_value = GetRealInitValue("chi_x",0.0);
    uniform_chi_y_value = GetRealInitValue("chi_y",0.0);
    uniform_chi_z_value = GetRealInitValue("chi_z",0.0);
    if(uniform_chi_x_value<0. || uniform_chi_y_value<0. || uniform_chi_z_value<0.) {
      char buf[4096];
      Oc_Snprintf(buf,sizeof(buf),
          "\n\"chi_x\", \"chi_y\", and \"chi_z\""
          " (%g,%g,%g) cannot be smaller than zero.",
          uniform_chi_x_value, uniform_chi_y_value, uniform_chi_z_value);
      throw Oxs_ExtError(this,buf);
    }
    if(uniform_chi_x_value==0. && uniform_chi_y_value==0. && uniform_chi_z_value==0.) {
      char buf[4096];
      Oc_Snprintf(buf,sizeof(buf),
          "\nYou have to specify at least one susceptibility, "
          " \"chi_x\", \"chi_y\", or \"chi_z\" (%g,%g,%g).",
          uniform_chi_x_value, uniform_chi_y_value, uniform_chi_z_value);
      throw Oxs_ExtError(this,buf);
    }
    // As for now, we restrict ourselves to two equal susceptibilities,
    // according to Kazansteva, Phys. Rev. B 77, 184428 (2008),
    // somehow contrary to Garanin, Phys. Rev. B 55, 3050 (1997).
    if(!((uniform_chi_x_value==0. && uniform_chi_y_value==uniform_chi_z_value) ||
         (uniform_chi_y_value==0. && uniform_chi_x_value==uniform_chi_z_value) ||
         (uniform_chi_z_value==0. && uniform_chi_x_value==uniform_chi_y_value)
      )) {
      String msg=String("\nAccording to Kazansteva, Phys. Rev. B 77, 184428 (2008),"
          " currently the only suppported case is, when two out of three:"
          " {\"chi_x\", \"chi_y\", \"chi_z\"} are specified."
          " Additionally, these two values must be equal.\n");
      throw Oxs_ExtError(this,msg.c_str());
    }
  }
  
  String integration_request = GetStringInitValue("integration","rect");
  if(integration_request.compare("rect")==0) {
    integration_method = RECT_INTEG;
  } else if(integration_request.compare("quad")==0) {
    integration_method = QUAD_INTEG;
  } else if(integration_request.compare("LLB")==0) {
    integration_method = LLB_INTEG;
    // For LLB it is necessary to supply chi_perp
    if(0.0==uniform_chi_x_value && 0.0==uniform_chi_y_value && 0.0==uniform_chi_z_value) {
      String msg=String("\nYou must specify susceptibility"
          " (via \"chi_x\", \"chi_y\", or \"chi_z\"),"
          "  for LLB integration mode.");
      throw Oxs_ExtError(this,msg.c_str());
    }
  } else {
    String msg=String("\nInvalid integration request: ")
      + integration_request
      + String("\n Should be either \"rect\", \"quad\" or \"LLB\".");
    throw Oxs_ExtError(this,msg.c_str());
  }
  VerifyAllInitArgsUsed();
}

OC_BOOL Klm_Anisotropy::Init()
{
  // LLB consistency
  if(VarMs_driver(this->director) && integration_method!=LLB_INTEG) {
    String msg=String("\nYou use LLB driver+evolver."
      " In such a case you have to set \"integration LLB\".");
    throw Oxs_ExtError(this,msg.c_str());
  }
  if(!VarMs_driver(this->director) && integration_method==LLB_INTEG) {
    String msg=String("\nYou must use LLB driver+evolver"
      " if you specify \"integration LLB\".");
    throw Oxs_ExtError(this,msg.c_str());
  }

  mesh_id = 0;
  K1.Release();
  axis.Release();
  return Oxs_ChunkEnergy::Init();
}

void Klm_Anisotropy::RectIntegEnergy
(const Oxs_SimState& state,
 Oxs_ComputeEnergyDataThreaded& ocedt,
 Oxs_ComputeEnergyDataThreadedAux& ocedtaux,
 OC_INDEX node_start,OC_INDEX node_stop
 ) const
{
  const Oxs_Mesh* mesh = state.mesh;
  const Oxs_MeshValue<OC_REAL8m>& Ms_inverse = *(state.Ms_inverse);
  const Oxs_MeshValue<ThreeVector>& spin = state.spin;

  Nb_Xpfloat energy_sum = 0.0;

  OC_REAL8m k = 0.0;
  if(K1_is_uniform) k = uniform_K1_value;

  ThreeVector unifaxis;
  if(axis_is_uniform) unifaxis = uniform_axis_value;

  for(OC_INDEX i=node_start;i<node_stop;++i) {
    if(!K1_is_uniform) k = K1[i];
    OC_REAL8m field_mult = (2.0/MU0)*k*Ms_inverse[i];
    if(field_mult==0.0) {
      if(ocedt.energy) (*ocedt.energy)[i] = 0.0;
      if(ocedt.H)      (*ocedt.H)[i].Set(0.,0.,0.);
      if(ocedt.mxH)    (*ocedt.mxH)[i].Set(0.,0.,0.);
      continue;
    }
    const ThreeVector& axisi = (axis_is_uniform ? unifaxis : axis[i]);
    if(k<=0) {
      // Easy plane (hard axis)
      OC_REAL8m dot = spin[i].x*axisi.x
        + spin[i].y*axisi.y + spin[i].z*axisi.z;
      ThreeVector H = (field_mult*dot)*axisi;
      OC_REAL8m ei = -k*dot*dot;

      OC_REAL8m tx = spin[i].y*H.z - spin[i].z*H.y; // mxH
      OC_REAL8m ty = spin[i].z*H.x - spin[i].x*H.z;
      OC_REAL8m tz = spin[i].x*H.y - spin[i].y*H.x;

      energy_sum += ei * mesh->Volume(i);
      if(ocedt.energy)       (*ocedt.energy)[i] = ei;
      if(ocedt.energy_accum) (*ocedt.energy_accum)[i] += ei;
      if(ocedt.H)       (*ocedt.H)[i] = H;
      if(ocedt.H_accum) (*ocedt.H_accum)[i] += H;
      if(ocedt.mxH)       (*ocedt.mxH)[i] = ThreeVector(tx,ty,tz);
      if(ocedt.mxH_accum) (*ocedt.mxH_accum)[i] += ThreeVector(tx,ty,tz);

    } else {
      // Easy axis case.  For improved accuracy, we want to report
      // energy as -k*(dot*dot-1), where dot = axis * spin.  But
      // dot*dot-1 suffers from bad loss of precision if spin is
      // nearly parallel to axis.  The are a couple of ways around
      // this.  Recall that both spin and axis are unit vectors.
      // Then from the cross product:
      //            (axis x spin)^2 = 1 - dot*dot
      // The cross product requires 6 mults and 3 adds, and
      // the norm squared takes 3 mult and 2 adds
      //            => 9 mults + 5 adds.
      // Another option is to use
      //            (axis - spin)^2 = 2*(1-dot) 
      //     so  1 - dot*dot = t*(2-t)
      //                where t = 0.5*(axis-spin)^2.
      // The op count here is 
      //            => 5 mults + 6 adds.
      // Another advantage to the second approach is you get 'dot', as
      // opposed to dot*dot, which saves a sqrt if dot is needed.  The
      // downside is that if axis and spin are anti-parallel, then you
      // want to use (axis+spin)^2 rather than (axis-spin)^2.  I did
      // some single-spin test runs and the performance of the two
      // methods was about the same.  Below we use the cross-product
      // formulation. -mjd, 28-Jan-2001
      OC_REAL8m dot = spin[i].x*axisi.x
        + spin[i].y*axisi.y + spin[i].z*axisi.z;
      OC_REAL8m scale = field_mult*dot;

      OC_REAL8m tx = spin[i].y*axisi.z - spin[i].z*axisi.y;
      OC_REAL8m ty = spin[i].z*axisi.x - spin[i].x*axisi.z;
      OC_REAL8m tz = spin[i].x*axisi.y - spin[i].y*axisi.x;

      OC_REAL8m ei = k*(tx*tx+ty*ty+tz*tz);
      energy_sum += ei * mesh->Volume(i);

      if(ocedt.energy)       (*ocedt.energy)[i] = ei;
      if(ocedt.energy_accum) (*ocedt.energy_accum)[i] += ei;
      if(ocedt.H)       (*ocedt.H)[i] = scale*axisi;
      if(ocedt.H_accum) (*ocedt.H_accum)[i] += scale*axisi;
      if(ocedt.mxH)       (*ocedt.mxH)[i]
                            = ThreeVector(scale*tx,scale*ty,scale*tz);
      if(ocedt.mxH_accum) (*ocedt.mxH_accum)[i]
                           += ThreeVector(scale*tx,scale*ty,scale*tz);
    }
  }
  ocedtaux.energy_total_accum += energy_sum.GetValue();
}

void Klm_Anisotropy::LLBIntegEnergy
(const Oxs_SimState& state,
 Oxs_ComputeEnergyDataThreaded& ocedt,
 Oxs_ComputeEnergyDataThreadedAux& ocedtaux,
 OC_INDEX node_start,OC_INDEX node_stop
 ) const
{
  const Oxs_Mesh* mesh = state.mesh;
  const Oxs_MeshValue<OC_REAL8m>& Ms = *(state.Ms);
  const Oxs_MeshValue<ThreeVector>& spin = state.spin;

  Nb_Xpfloat energy_sum = 0.0;

  const OC_REAL8m field_x_mult = (uniform_chi_x_value==0.0 ? 0.0 : 
                          -1.0/uniform_chi_x_value);
  const OC_REAL8m field_y_mult = (uniform_chi_y_value==0.0 ? 0.0 : 
                          -1.0/uniform_chi_y_value);
  const OC_REAL8m field_z_mult = (uniform_chi_z_value==0.0 ? 0.0 : 
                          -1.0/uniform_chi_z_value);

  const OC_REAL8m energy_x_mult = (uniform_chi_x_value==0.0 ? 0.0 : 
                       MU0/(2.0*uniform_chi_x_value));
  const OC_REAL8m energy_y_mult = (uniform_chi_y_value==0.0 ? 0.0 : 
                       MU0/(2.0*uniform_chi_y_value));
  const OC_REAL8m energy_z_mult = (uniform_chi_z_value==0.0 ? 0.0 : 
                       MU0/(2.0*uniform_chi_z_value));

  for(OC_INDEX i=node_start;i<node_stop;++i) {
    if(Ms[i]==0.0) {
      if(ocedt.energy) (*ocedt.energy)[i] = 0.0;
      if(ocedt.H)      (*ocedt.H)[i].Set(0.,0.,0.);
      if(ocedt.mxH)    (*ocedt.mxH)[i].Set(0.,0.,0.);
      continue;
    }
    {
      // According to Garanin, Phys. Rev. B 55, 3050 (1997):
      //   E = mu0/2 * (Mx^2/chi_x + My^2/chi_y + Mz^2/chi_z)
      //   H = (-Mx/chi_x, -My/chi_y, -Mz/chi_z)
      
      const OC_REAL8m Mx = spin[i].x*Ms[i];
      const OC_REAL8m My = spin[i].y*Ms[i];
      const OC_REAL8m Mz = spin[i].z*Ms[i];
      const ThreeVector H = ThreeVector(field_x_mult*Mx, field_y_mult*My, field_z_mult*Mz);
      
      const OC_REAL8m ei = energy_x_mult*Mx*Mx + energy_y_mult*My*My + energy_z_mult*Mz*Mz;
      energy_sum += ei * mesh->Volume(i);

      if(ocedt.energy)        (*ocedt.energy)[i]        = ei;
      if(ocedt.energy_accum)  (*ocedt.energy_accum)[i] += ei;
      if(ocedt.H)             (*ocedt.H)[i]             = H;
      if(ocedt.H_accum)       (*ocedt.H_accum)[i]      += H;
      if(ocedt.mxH)           (*ocedt.mxH)[i]           = spin[i]^H;
      if(ocedt.mxH_accum)     (*ocedt.mxH_accum)[i]    += spin[i]^H;
    }
  }
  ocedtaux.energy_total_accum += energy_sum.GetValue();
}

void Klm_Anisotropy::ComputeEnergyChunk
(const Oxs_SimState& state,
 Oxs_ComputeEnergyDataThreaded& ocedt,
 Oxs_ComputeEnergyDataThreadedAux& ocedtaux,
 OC_INDEX node_start,OC_INDEX node_stop,
 int threadnumber
 ) const
{
  if(node_stop>state.mesh->Size() || node_start>node_stop) {
    throw Oxs_ExtError(this,"Programming error:"
                       " Invalid node_start/node_stop values");
  }

  if(mesh_id !=  state.mesh->Id()) {
    // This is either the first pass through, or else mesh
    // has changed.  Initialize/update data fields.
    // NB: At a lower level, this may potentially involve calls back
    // into the Tcl interpreter.  Per Tcl spec, only the thread
    // originating the interpreter is allowed to make calls into it, so
    // only threadnumber == 0 can do this processing.  Any other thread
    // must block until that processing is complete.
    thread_control.Lock();
    if(Oxs_ThreadError::IsError()) {
      if(thread_control.count>0) {
        // Release a blocked thread
        thread_control.Notify();
      }
      thread_control.Unlock();
      return; // What else?
    }
    if(threadnumber != 0) {
      if(mesh_id != state.mesh->Id()) {
        // If above condition is false, then the main thread came
        // though and initialized everything between the time of
        // the previous check and this thread's acquiring of the
        // thread_control mutex; in which case, "never mind".
        // Otherwise:
        ++thread_control.count; // Multiple threads may progress to this
        /// point before the main thread (threadnumber == 0) grabs the
        /// thread_control mutex.  Keep track of how many, so that
        /// afterward they may be released, one by one.  (The main
        /// thread will Notify control_wait.cond once; after that
        /// as each waiting thread is released, the newly released
        /// thread sends a Notify to wake up the next one.
        thread_control.Wait(0);
        --thread_control.count;
        int condcheckerror=0;
        if(mesh_id !=  state.mesh->Id()) {
          // Error?
          condcheckerror=1;
          Oxs_ThreadPrintf(stderr,"Invalid condition in"
                           " Klm_Anisotropy::ComputeEnergyChunk(),"
                           " thread number %d\n",threadnumber);
        }
        if(thread_control.count>0) {
          // Free a waiting thread.
          thread_control.Notify();
        }
        thread_control.Unlock();
        if(condcheckerror || Oxs_ThreadError::IsError()) {
          return; // What else?
        }
      } else {
        if(thread_control.count>0) {
          // Free a waiting thread.  (Actually, it can occur that the
          // thread_control will be grabbed by another thread that is
          // blocked at the first thread_control mutex Lock() call above
          // rather than on the ConditionWait, in which case this
          // ConditionNotify will be effectively lost.  But that is
          // okay, because then *that* thread will Notify when it
          // releases the mutex.)
          thread_control.Notify();
        }
        thread_control.Unlock();
      }
    } else {
      // Main thread (threadnumber == 0)
      try {
        // KL(m) Do it only for non-LLB mode
        if(integration_method != LLB_INTEG) {
          K1_init->FillMeshValue(state.mesh,K1);
          axis_init->FillMeshValue(state.mesh,axis);
          const OC_UINT4m size = state.mesh->Size();
          for(OC_UINT4m i=0;i<size;i++) {
            // Check that axis is a unit vector:
            const OC_REAL8m eps = 1e-14;
            if(axis[i].MagSq()<eps*eps) {
              throw Oxs_ExtError(this,"Invalid initialization detected:"
                                 " Zero length anisotropy axis");
            } else {
              axis[i].MakeUnit();
            }
          }
        }
        mesh_id = state.mesh->Id();
      } catch(Oxs_ExtError& err) {
        // Leave unmatched mesh_id as a flag to check
        // Oxs_ThreadError for an error.
        Oxs_ThreadError::SetError(String(err));
        if(thread_control.count>0) {
          thread_control.Notify();
        }
        thread_control.Unlock();
        throw;
      } catch(String& serr) {
        // Leave unmatched mesh_id as a flag to check
        // Oxs_ThreadError for an error.
        Oxs_ThreadError::SetError(serr);
        if(thread_control.count>0) {
          thread_control.Notify();
        }
        thread_control.Unlock();
        throw;
      } catch(const char* cerr) {
        // Leave unmatched mesh_id as a flag to check
        // Oxs_ThreadError for an error.
        Oxs_ThreadError::SetError(String(cerr));
        if(thread_control.count>0) {
          thread_control.Notify();
        }
        thread_control.Unlock();
        throw;
      } catch(...) {
        // Leave unmatched mesh_id as a flag to check
        // Oxs_ThreadError for an error.
        Oxs_ThreadError::SetError(String("Error in "
            "Klm_Anisotropy::ComputeEnergyChunk"));
        if(thread_control.count>0) {
          thread_control.Notify();
        }
        thread_control.Unlock();
        throw;
      }
      if(thread_control.count>0) {
        // Free a waiting thread.  (Actually, it can occur that the
        // thread_control will be grabbed by another thread that is
        // blocked at the first thread_control mutex Lock() call above
        // rather than on the ConditionWait, in which case this
        // ConditionNotify will be effectively lost.  But that is
        // okay, because then *that* thread will Notify when it
        // releases the mutex.)
        thread_control.Notify();
      }
      thread_control.Unlock();
    }
  }

  if(integration_method == RECT_INTEG) { // KL(m)
    RectIntegEnergy(state,ocedt,ocedtaux,node_start,node_stop);
    return;
  }
  else if(integration_method == LLB_INTEG) { // KL(m)
    LLBIntegEnergy(state,ocedt,ocedtaux,node_start,node_stop);
    return;
  }

  // Otherwise, use QUAD_INTEG.  This implementation is not especially
  // efficient, and could presumably be sped up by a good bit if
  // needed.

  // The QUAD_INTEG code requires a rectangular mesh
  const Oxs_RectangularMesh* mesh
    = dynamic_cast<const Oxs_RectangularMesh*>(state.mesh);
  if(mesh==NULL) {
    throw Oxs_ExtError(this,
          "Import mesh to Klm_Anisotropy::GetEnergy()"
          " is not an Oxs_RectangularMesh object.  Quad integration"
          " method requires a rectangular mesh.");
  }

  // Do base evaluations, writing results into scratch space
  Oxs_ComputeEnergyDataThreaded work_ocedt(state);
  Oxs_ComputeEnergyDataThreadedAux work_ocedtaux;
  if(ocedt.energy) {
    work_ocedt.energy = work_ocedt.scratch_energy = ocedt.energy;
  } else {
    ocedt.scratch_energy->AdjustSize(state.mesh); // Thread-safe
    work_ocedt.energy = work_ocedt.scratch_energy = ocedt.scratch_energy;
  }
  if(ocedt.H) {
    work_ocedt.H      = work_ocedt.scratch_H      = ocedt.H;
  } else {
    ocedt.scratch_H->AdjustSize(state.mesh); // Thread-safe
    work_ocedt.H      = work_ocedt.scratch_H      = ocedt.scratch_H;
  }
  RectIntegEnergy(state,work_ocedt,work_ocedtaux,node_start,node_stop);

  // Edge correction if higher-order integration method requested.
  // See mjd's NOTES II, pp 178-181, Aug-2002.
  // NOTE: For short dimension lengths, all cells are modified.
  Oxs_MeshValue<OC_REAL8m>& energy = *work_ocedt.energy;
  Oxs_MeshValue<ThreeVector>& field = *work_ocedt.H;

  const OC_INDEX xdim = mesh->DimX();
  const OC_INDEX ydim = mesh->DimY();
  const OC_INDEX zdim = mesh->DimZ();
  const OC_INDEX xydim = xdim*ydim;
  OC_INDEX x,y,z;

  OC_INDEX xstart,ystart,zstart;
  mesh->GetCoords(node_start,xstart,ystart,zstart);
  OC_INDEX xstop,ystop,zstop;
  mesh->GetCoords(node_stop-1,xstop,ystop,zstop);


  // x-axis.  Note special case handling for short lengths.  Also note
  // very lazy handling of node start/stop conditions.  With luck, the
  // compiler can figure it out.  Otherwise, if we ever care to have
  // this code run fast, the condition handling should be re-written.
  if(xdim>=6) {
    for(z=zstart;z<=zstop;++z) {
      OC_INDEX ya = (z==zstart ? ystart : 0);
      OC_INDEX yb = (z==zstop  ? ystop : ydim-1);
      for(y=ya;y<=yb;++y) {
        OC_INDEX i = mesh->Index(0,y,z); // Get base linear address
        if(node_start<=i && i<node_stop) {
          energy[i]   *= 26./24.; // Left face
          field[i]    *= 26./24.;
        }
        if(node_start<=i+1 && i+1<node_stop) {
          energy[i+1] *= 21./24.;
          field[i+1]  *= 21./24.;
        }
        if(node_start<=i+2 && i+2<node_stop) {
          energy[i+2] *= 25./24.;
          field[i+2]  *= 25./24.;
        }

        i += xdim-3;
        if(node_start<=i && i<node_stop) {
          energy[i]   *= 25./24.;
          field[i]    *= 25./24.;
        }
        if(node_start<=i+1 && i+1<node_stop) {
          energy[i+1] *= 21./24.;
          field[i+1]  *= 21./24.;
        }
        if(node_start<=i+2 && i+2<node_stop) {
          energy[i+2] *= 26./24.; // Right face
          field[i+2]  *= 26./24.;
        }
      }
    }
  } else if(xdim==5) {
    for(z=zstart;z<=zstop;++z) {
      OC_INDEX ya = (z==zstart ? ystart : 0);
      OC_INDEX yb = (z==zstop  ? ystop : ydim-1);
      OC_INDEX i = mesh->Index(0,ya,z); // Get base linear address
      for(y=ya;y<=yb;++y) {
        if(node_start<=i && i<node_stop) {
          energy[i]   *= 26./24.;
          field[i]    *= 26./24.;
        }
        if(node_start<=i+1 && i+1<node_stop) {
          energy[i+1] *= 21./24.;
          field[i+1]  *= 21./24.;
        }
        if(node_start<=i+2 && i+2<node_stop) {
          energy[i+2] *= 26./24.;
          field[i+2]  *= 26./24.;
        }
        if(node_start<=i+3 && i+3<node_stop) {
          energy[i+3] *= 21./24.;
          field[i+3]  *= 21./24.;
        }
        if(node_start<=i+4 && i+4<node_stop) {
          energy[i+4] *= 26./24.;
          field[i+4]  *= 26./24.;
        }
        i += 5;
      }
    }
  } else if(xdim==4) {
    for(z=zstart;z<=zstop;++z) {
      OC_INDEX ya = (z==zstart ? ystart : 0);
      OC_INDEX yb = (z==zstop  ? ystop : ydim-1);
      OC_INDEX i = mesh->Index(0,ya,z); // Get base linear address
      for(y=ya;y<=yb;++y) {
        if(node_start<=i && i<node_stop) {
          energy[i]   *= 26./24.;
          field[i]    *= 26./24.;
        }
        if(node_start<=i+1 && i+1<node_stop) {
          energy[i+1] *= 22./24.;
          field[i+1]  *= 22./24.;
        }
        if(node_start<=i+2 && i+2<node_stop) {
          energy[i+2] *= 22./24.;
          field[i+2]  *= 22./24.;
        }
        if(node_start<=i+3 && i+3<node_stop) {
          energy[i+3] *= 26./24.;
          field[i+3]  *= 26./24.;
        }
        i += 4;
      }
    }
  } else if(xdim==3) {
    for(z=zstart;z<=zstop;++z) {
      OC_INDEX ya = (z==zstart ? ystart : 0);
      OC_INDEX yb = (z==zstop  ? ystop : ydim-1);
      OC_INDEX i = mesh->Index(0,ya,z); // Get base linear address
      for(y=ya;y<=yb;++y) {
        if(node_start<=i && i<node_stop) {
          energy[i]   *= 27./24.;
          field[i]    *= 27./24.;
        }
        if(node_start<=i+1 && i+1<node_stop) {
          energy[i+1] *= 18./24.;
          field[i+1]  *= 18./24.;
        }
        if(node_start<=i+2 && i+2<node_stop) {
          energy[i+2] *= 27./24.;
          field[i+2]  *= 27./24.;
        }
        i += 3;
      }
    }
  }
  // Quadratic fit requires 3 points, so no higher order method
  // available if xdim<3.
    
  // y-axis.  Note special case handling for short lengths.
  if(ydim>=6) {
    for(z=zstart;z<=zstop;++z) {
      // Front face
      OC_INDEX i = mesh->Index(0,0,z); // Get base linear address
      for(x=0;x<xdim;x++) { // y==0
        if(node_start<=i && i<node_stop) {
          energy[i] *= 26./24.;
          field[i]  *= 26./24.;
        }
        ++i;
      } // NB: At end of loop, i wraps around to next x-row.
      for(x=0;x<xdim;x++) { // y==1
        if(node_start<=i && i<node_stop) {
          energy[i] *= 21./24.;
          field[i]  *= 21./24.;
        }
        ++i;
      }
      for(x=0;x<xdim;x++) { // y==2
        if(node_start<=i && i<node_stop) {
          energy[i] *= 25./24.;
          field[i]  *= 25./24.;
        }
        ++i;
      }
      // Back face
      i = mesh->Index(0,ydim-3,z);
      for(x=0;x<xdim;x++) { // y==ydim-3
        if(node_start<=i && i<node_stop) {
          energy[i] *= 25./24.;
          field[i]  *= 25./24.;
        }
        ++i;
      }
      for(x=0;x<xdim;x++) { // y==ydim-2
        if(node_start<=i && i<node_stop) {
          energy[i] *= 21./24.;
          field[i]  *= 21./24.;
        }
        ++i;
      }
      for(x=0;x<xdim;x++) { // y==ydim-1
        if(node_start<=i && i<node_stop) {
          energy[i] *= 26./24.;
          field[i]  *= 26./24.;
        }
        ++i;
      }
    }
  } else if(ydim==5) {
    for(z=zstart;z<=zstop;++z) {
      OC_INDEX i = mesh->Index(0,0,z); // Get base linear address
      for(x=0;x<xdim;x++) { // y==0
        if(node_start<=i && i<node_stop) {
          energy[i] *= 26./24.;
          field[i]  *= 26./24.;
        }
        ++i;
      }
      for(x=0;x<xdim;x++) { // y==1
        if(node_start<=i && i<node_stop) {
          energy[i] *= 21./24.;
          field[i]  *= 21./24.;
        }
        ++i;
      }
      for(x=0;x<xdim;x++) { // y==2
        if(node_start<=i && i<node_stop) {
          energy[i] *= 26./24.;
          field[i]  *= 26./24.;
        }
        ++i;
      }
      for(x=0;x<xdim;x++) { // y==3
        if(node_start<=i && i<node_stop) {
          energy[i] *= 21./24.;
          field[i]  *= 21./24.;
        }
        ++i;
      }
      for(x=0;x<xdim;x++) { // y==4
        if(node_start<=i && i<node_stop) {
          energy[i] *= 26./24.;
          field[i]  *= 26./24.;
        }
        ++i;
      }
    }
  } else if(ydim==4) {
    for(z=zstart;z<=zstop;++z) {
      OC_INDEX i = mesh->Index(0,0,z); // Get base linear address
      for(x=0;x<xdim;x++) { // y==0
        if(node_start<=i && i<node_stop) {
          energy[i] *= 26./24.;
          field[i]  *= 26./24.;
        }
        ++i;
      }
      for(x=0;x<xdim;x++) { // y==1
        if(node_start<=i && i<node_stop) {
          energy[i] *= 22./24.;
          field[i]  *= 22./24.;
        }
        ++i;
      }
      for(x=0;x<xdim;x++) { // y==2
        if(node_start<=i && i<node_stop) {
          energy[i] *= 22./24.;
          field[i]  *= 22./24.;
        }
        ++i;
      }
      for(x=0;x<xdim;x++) { // y==3
        if(node_start<=i && i<node_stop) {
          energy[i] *= 26./24.;
          field[i]  *= 26./24.;
        }
        ++i;
      }
    }
  } else if(ydim==3) {
    for(z=zstart;z<=zstop;++z) {
      OC_INDEX i = mesh->Index(0,0,z); // Get base linear address
      for(x=0;x<xdim;x++) { // y==0
        if(node_start<=i && i<node_stop) {
          energy[i] *= 27./24.;
          field[i]  *= 27./24.;
        }
        ++i;
      }
      for(x=0;x<xdim;x++) { // y==1
        if(node_start<=i && i<node_stop) {
          energy[i] *= 18./24.;
          field[i]  *= 18./24.;
        }
        ++i;
      }
      for(x=0;x<xdim;x++) { // y==2
        if(node_start<=i && i<node_stop) {
          energy[i] *= 27./24.;
          field[i]  *= 27./24.;
        }
        ++i;
      }
    }
  }
  // Quadratic fit requires 3 points, so no higher order method
  // available if ydim<3.

    
  // z-axis.  Note special case handling for short lengths.
  if(zdim>=6) {
    // Bottom face, z==0
    OC_INDEX i = node_start;
    while(i<xydim && i<node_stop) { // z==0
      energy[i] *= 26./24.;
      field[i]  *= 26./24.;
      ++i;
    } // NB: At end of loop, i wraps around to next xy-plane.
    while(i<2*xydim && i<node_stop) { // z==1
      energy[i] *= 21./24.;
      field[i]  *= 21./24.;
      ++i;
    }
    while(i<3*xydim && i<node_stop) { // z==2
      energy[i] *= 25./24.;
      field[i]  *= 25./24.;
      ++i;
    }
    // Top face, z==zdim-1
    OC_INDEX itop = mesh->Index(0,0,zdim-3);
    i = (node_start<itop ? itop : node_start);
    while(i<itop+xydim && i<node_stop) { // z==zdim-3
      energy[i] *= 25./24.;
      field[i]  *= 25./24.;
      ++i;
    }
    while(i<itop+2*xydim && i<node_stop) { // z==zdim-2
      energy[i] *= 21./24.;
      field[i]  *= 21./24.;
      ++i;
    }
    while(i<node_stop) { // z==zdim-1
      energy[i] *= 26./24.;
      field[i]  *= 26./24.;
      ++i;
    }
  } else if(zdim==5) {
    OC_INDEX i = node_start;
    while(i<xydim && i<node_stop) { // z==0; bottom face
      energy[i] *= 26./24.;
      field[i]  *= 26./24.;
      ++i;
    }
    while(i<2*xydim && i<node_stop) { // z==1
      energy[i] *= 21./24.;
      field[i]  *= 21./24.;
      ++i;
    }
    while(i<3*xydim && i<node_stop) { // z==2
      energy[i] *= 26./24.;
      field[i]  *= 26./24.;
      ++i;
    }
    while(i<4*xydim && i<node_stop) { // z==3
      energy[i] *= 21./24.;
      field[i]  *= 21./24.;
      ++i;
    }
    while(i<node_stop) { // z==4; top face
      energy[i] *= 26./24.;
      field[i]  *= 26./24.;
      ++i;
    }
  } else if(zdim==4) {
    OC_INDEX i = node_start;
    while(i<xydim && i<node_stop) { // z==0; bottom face
      energy[i] *= 26./24.;
      field[i]  *= 26./24.;
      ++i;
    }
    while(i<2*xydim && i<node_stop) { // z==1
      energy[i] *= 22./24.;
      field[i]  *= 22./24.;
      ++i;
    }
    while(i<3*xydim && i<node_stop) { // z==2
      energy[i] *= 22./24.;
      field[i]  *= 22./24.;
      ++i;
    }
    while(i<node_stop) { // z==3; top face
      energy[i] *= 26./24.;
      field[i]  *= 26./24.;
      ++i;
    }
  } else if(zdim==3) {
    OC_INDEX i = node_start;
    while(i<xydim && i<node_stop) { // z==0; bottom face
      energy[i] *= 27./24.;
      field[i]  *= 27./24.;
      ++i;
    }
    while(i<2*xydim && i<node_stop) { // z==1
      energy[i] *= 18./24.;
      field[i]  *= 18./24.;
      ++i;
    }
    while(i<node_stop) { // z==2; top face
      energy[i] *= 27./24.;
      field[i]  *= 27./24.;
      ++i;
    }
  }
  // Quadratic fit requires 3 points, so no higher order method
  // available if zdim<3.

  // Copy from scratch to real buffers.
  Nb_Xpfloat energy_sum = 0.0;
  if(ocedt.energy_accum) {
    for(OC_INDEX i=node_start;i<node_stop;i++) {
      (*ocedt.energy_accum)[i] += energy[i];
      energy_sum += energy[i];
    }
  } else {
    for(OC_INDEX i=node_start;i<node_stop;i++) {
      energy_sum += energy[i];
    }
  }
  ocedtaux.energy_total_accum += energy_sum.GetValue() * mesh->Volume(0); // All cells
  /// have same volume in an Oxs_RectangularMesh.

  const Oxs_MeshValue<ThreeVector>& spin = state.spin;
  if(ocedt.mxH_accum && ocedt.H_accum) {
    for(OC_INDEX i=node_start; i<node_stop; ++i) {
      (*ocedt.H_accum)[i] += field[i];
      ThreeVector temp = spin[i];
      temp ^= field[i];
      (*ocedt.mxH_accum)[i] += temp;
    }
  } else if(ocedt.mxH_accum) {
    for(OC_INDEX i=node_start; i<node_stop; ++i) {
      ThreeVector temp = spin[i];
      temp ^= field[i];
      (*ocedt.mxH_accum)[i] += temp;
    }
  } else if(ocedt.H_accum) {
    for(OC_INDEX i=node_start; i<node_stop; ++i) {
      (*ocedt.H_accum)[i] += field[i];
    }
  }
  if(ocedt.mxH) {
    for(OC_INDEX i=node_start; i<node_stop; ++i) {
      ThreeVector temp = spin[i];
      temp ^= field[i];
      (*ocedt.mxH)[i] = temp;
    }
  }
}
