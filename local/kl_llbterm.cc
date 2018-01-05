/* FILE: uuanisotropy.cc            -*-Mode: c++-*- 
 *
 * Uniform Uniaxial Anisotropy, derived from Oxs_Energy class.
 *
 */

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
#include "rectangularmesh.h"
#include "kl_llbterm.h"
#include "kl_llb_util.h"
#include "energy.h"		// Needed to make MSVC++ 5 happy

OC_USE_STRING;

// Oxs_Ext registration support
OXS_EXT_REGISTER(Klm_LLB_Term);

/* End includes */


// Constructor
Klm_LLB_Term::Klm_LLB_Term(
  const char* name,     // Child instance id
  Oxs_Director* newdtr, // App director
  const char* argstr)   // MIF input block parameters
  : Oxs_ChunkEnergy(name,newdtr,argstr),
    chi_parallel(0.0) // KL(m)
{
  // Process arguments
  chi_parallel = GetRealInitValue("chi_parallel");
  if(chi_parallel<=0.0) {
    char buf[4096];
    Oc_Snprintf(buf,sizeof(buf),
		  "\nInvalid parameter value:"
		  " Specified \"chi_parallel\" is %g (should be >0.0)",
		  chi_parallel);
    throw Oxs_Ext::Error(this,buf);
  }

  VerifyAllInitArgsUsed();
  
  // LLB_Term consistency
  Set_Use_LLB_Term();
}

OC_BOOL Klm_LLB_Term::Init()
{
  // LLB consistency
  if(!VarMs_driver(this->director)) {
    String msg=String("\nYou must use LLB driver+evolver"
      " if you use Klm_LLB_Term.");
    throw Oxs_ExtError(this,msg.c_str());
  }

  return Oxs_ChunkEnergy::Init();
}

void Klm_LLB_Term::ComputeEnergyChunk
(const Oxs_SimState& state,
 Oxs_ComputeEnergyDataThreaded& ocedt,
 Oxs_ComputeEnergyDataThreadedAux& ocedtaux,
 OC_INDEX node_start,OC_INDEX node_stop,
 int /* threadnumber */
 ) const
{
  const OC_INDEX size = state.mesh->Size();
  if(size<1) {
    return;
  }

  const Oxs_MeshValue<ThreeVector>& spin = state.spin;
  const Oxs_MeshValue<OC_REAL8m>& Ms = *(state.Ms);
  const Oxs_MeshValue<OC_REAL8m>& Me_T = *(state.GetPtr_Me_T());

  Nb_Xpfloat energy_sum = 0.0;
  for(OC_INDEX i=node_start;i<node_stop;++i) {
    if(0.0 == Ms[i]) {
      if(ocedt.energy) (*ocedt.energy)[i] = 0.0;
      if(ocedt.H)      (*ocedt.H)[i].Set(0.,0.,0.);
      if(ocedt.mxH)    (*ocedt.mxH)[i].Set(0.,0.,0.);
    } else {
    
      const OC_REAL8m chi_2_MeT2 = 1.0/(2.0*chi_parallel*Me_T[i]*Me_T[i]);
      const OC_REAL8m me_mi = Me_T[i]*Me_T[i] - Ms[i]*Ms[i];
      const OC_REAL8m ei = (chi_2_MeT2/4.0)*MU0*me_mi*me_mi;
      const ThreeVector hi = chi_2_MeT2*me_mi*Ms[i]*spin[i];

      energy_sum += ei * state.mesh->Volume(i);

      if(ocedt.energy)       (*ocedt.energy)[i] = ei;
      if(ocedt.energy_accum) (*ocedt.energy_accum)[i] += ei;
      if(ocedt.H)       (*ocedt.H)[i] = hi;
      if(ocedt.H_accum) (*ocedt.H_accum)[i] += hi;
      // m || H, thus mxH= 0
      if(ocedt.mxH)       (*ocedt.mxH)[i] = ThreeVector(0.,0.,0.); // probably not needed here
      //if(ocedt.mxH_accum) (*ocedt.mxH_accum)[i] += ThreeVector(0.0,0.0,0.0);
    }
  }

  ocedtaux.energy_total_accum += energy_sum.GetValue();
}
