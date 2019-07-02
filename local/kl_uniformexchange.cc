/* FILE: kluniformexchange.cc            -*-Mode: c++-*- 
 *
 * see: kl_uniformexchange.h
 *
 */

#include <ctype.h>

#include <string>

#include "nb.h"
#include "director.h"
#include "mesh.h"
#include "meshvalue.h"
#include "oxswarn.h"
#include "simstate.h"
#include "threevector.h"
#include "rectangularmesh.h"
#include "kl_uniformexchange.h"
#include "kl_llb_util.h"
#include "energy.h"             // Needed to make MSVC++ 5 happy

OC_USE_STRING;

// Oxs_Ext registration support
OXS_EXT_REGISTER(Klm_UniformExchange);

/* End includes */

// Revision information, set via CVS keyword substitution
static const Oxs_WarningMessageRevisionInfo revision_info
  (__FILE__,
   "$Revision: 1.75 $",
   "$Date: 2016/03/04 22:08:58 $",
   "$Author: lebecki $",
   "Kristof M. Lebecki (lebecki@fuw.edu.pl)");

#undef SKIP_PROBLEM_CODE
/// The 8.0 release of the Intel icc/ia64/Lintel compiler (and perhaps
/// others) gets hung up trying to compile the code wrapped up below in
/// the "SKIP_PROBLEM_CODE" #ifdef sections, at least when "-O3"
/// optimization is enabled.  These blocks provide features that are
/// undocumented at this time (April-2004), so for now we just disable
/// compiling of these blocks.  Remove the above #define to enable
/// compilation.

// Constructor
Klm_UniformExchange::Klm_UniformExchange(
  const char* name,     // Child instance id
  Oxs_Director* newdtr, // App director
  const char* argstr)   // MIF input block parameters
  : Oxs_ChunkEnergy(name,newdtr,argstr),
    excoeftype(A_UNKNOWN), A(-1.), lex(-1.),
    kernel(NGBR_UNKNOWN), 
    xperiodic(0),yperiodic(0),zperiodic(0),
    mesh_id(0), energy_density_error_estimate(-1)
{
  // Process arguments
  OC_BOOL has_A = HasInitValue("A");
  OC_BOOL has_lex = HasInitValue("lex");
  if(has_A && has_lex) {
    throw Oxs_ExtError(this,"Invalid exchange coefficient request:"
			 " both A and lex specified; only one should"
			 " be given.");
  } else if(has_lex) {
    lex = GetRealInitValue("lex");
    excoeftype = LEX_TYPE;
  } else {
    A = GetRealInitValue("A");
    excoeftype = A_TYPE;
  }

  String kernel_request = GetStringInitValue("kernel","6ngbr");
  if(kernel_request.compare("6ngbrLLB")==0) {
    kernel = NGBR_6_LLB;
  } else {
    String msg=String("Invalid kernel request: ")
      + kernel_request
      + String("\n Should be 6ngbrLLB."); // KL(m)
    throw Oxs_ExtError(this,msg.c_str());
  }

//  if(A_TYPE != excoeftype && NGBR_6_MIRROR != kernel) {
//    throw Oxs_ExtError(this,"Invalid exchange coefficient+kernel"
//			 " combination; lex specification currently"
//			 " only supported using 6ngbr kernel.");
//  }

// KL(m)
  if(A_TYPE != excoeftype) {
    throw Oxs_Ext::Error(this,"Invalid exchange coefficient+kernel"
			 " combination; lex specification currently"
			 " not supported.");
  }

  VerifyAllInitArgsUsed();

  // Setup outputs, KL(m)
  min_m_output.Setup(this,InstanceName(),"Min m","",1,
			    &Klm_UniformExchange::UpdateDerivedOutputs);
  min_m_output.Register(director,0);  
  stage_min_m_output.Setup(this,InstanceName(),"Stage Min m","",1,
			    &Klm_UniformExchange::UpdateDerivedOutputs);
  stage_min_m_output.Register(director,0);  
  run_min_m_output.Setup(this,InstanceName(),"Run Min m","",1,
			    &Klm_UniformExchange::UpdateDerivedOutputs);
  run_min_m_output.Register(director,0);  
  //
  max_m_output.Setup(this,InstanceName(),"Max m","",1,
			    &Klm_UniformExchange::UpdateDerivedOutputs);
  max_m_output.Register(director,0);  
  stage_max_m_output.Setup(this,InstanceName(),"Stage Max m","",1,
			    &Klm_UniformExchange::UpdateDerivedOutputs);
  stage_max_m_output.Register(director,0);  
  run_max_m_output.Setup(this,InstanceName(),"Run Max m","",1,
			    &Klm_UniformExchange::UpdateDerivedOutputs);
  run_max_m_output.Register(director,0);  
  //
  maxMdiff_output.Setup(this,InstanceName(),"Max m Diff","",1,
			    &Klm_UniformExchange::UpdateDerivedOutputs);
  maxMdiff_output.Register(director,0);
  stage_maxMdiff_output.Setup(this,InstanceName(),"Stage Max m Diff","",1,
			    &Klm_UniformExchange::UpdateDerivedOutputs);
  stage_maxMdiff_output.Register(director,0);
  run_maxMdiff_output.Setup(this,InstanceName(),"Run Max m Diff","",1,
			    &Klm_UniformExchange::UpdateDerivedOutputs);
  run_maxMdiff_output.Register(director,0);
}

Klm_UniformExchange::~Klm_UniformExchange()
{}

OC_BOOL Klm_UniformExchange::Init()
{
  // LLB consistency
  if(VarMs_driver(this->director) && 
                  kernel!=NGBR_6_LLB) {
    String msg=String("\nYou use LLB driver+evolver."
      " In such a case you have to set kernel to"
      " \"6ngbrLLB\".");
    throw Oxs_ExtError(this,msg.c_str());
  }
  if(!VarMs_driver(this->director) && 
                   (kernel==NGBR_6_LLB)) {
    String msg=String("\nYou must use LLB driver+evolver"
      " if you specify kernel \"6ngbrLLB\".");
    throw Oxs_ExtError(this,msg.c_str());
  }
  
  mesh_id = 0;
  xcoef.Free();
  ycoef.Free();
  zcoef.Free();
  energy_density_error_estimate = -1;
  return Oxs_ChunkEnergy::Init();
}

// KL(m)
void
Klm_UniformExchange::CalcEnergy6NgbrMirror_LLB
// Usually*, see J. Miltat and M. J. Donahue, in Handbook of Magnetism and Advanced 
//   Magnetic Materials, edited by H. Kronmüller and S. Parkin (Wiley-Interscience, 
//   Chichester, 2007), Vol. 2, p. 742. Eq.(64,65):
//   E_i = -A.|V_i| \sum_j m_i.(m_j-m_i)/delta_ij^2,
//   H_i = 2A/(mu_0.M_s) \sum_j (m_j-m_i)/delta_ij^2,
//   where delta_ij = EdgeLengthX(), EdgeLengthY(), or EdgeLengthZ(),
//         |m| = 1.
// For LLB we transform it to:
//   E_i = -A.|V_i| \sum_j M_i.(M_j-M_i)/(delta_ij^2.Me(T)^2),
//   H_i = 2A/(mu_0.Me(T)^2) \sum_j (M_j-M_i)/delta_ij^2,
// Implementation changes in:     NEW    >instead-of>   OLD
// - hmult :                      -2/mu0                -2/(mu0*Ms[i])
// - base, sum (xlast, temp*):    Ms[.]*spin[.]         spin[.]
// - dot:                         spin[.]*base/Ms[i]    spin[.]*base
// - tx, ty, tz:                  base/Ms[i] x sum      base x sum
// * This is a field-based approach, see Miltat+Donahue. For energy-based
//   approach so-called 26-neighbors schema is necessary, see M. J. Donahue and D. G. 
//   Porter, Physica B 343, 177 (2004).

(const Oxs_MeshValue<ThreeVector>& spin,
 const Oxs_MeshValue<OC_REAL8m>& Ms_inverse,
 const Oxs_MeshValue<OC_REAL8m>& Ms, // KL(m)
 const Oxs_MeshValue<OC_REAL8m>& Me_T, // KL(m)
 const Oxs_CommonRectangularMesh* mesh,
 Oxs_ComputeEnergyDataThreaded& ocedt,
 Oxs_ComputeEnergyDataThreadedAux& ocedtaux,
 OC_INDEX node_start,OC_INDEX node_stop,
 int threadnumber) const
{

#ifndef NDEBUG
  if(node_stop>mesh->Size() || node_start>node_stop) {
    throw Oxs_ExtError(this,"Programming error:"
                       " Invalid node_start/node_stop values");
  }
#endif

  OC_INDEX xdim = mesh->DimX();
  OC_INDEX ydim = mesh->DimY();
  OC_INDEX zdim = mesh->DimZ();
  OC_INDEX xydim = xdim*ydim;
  OC_INDEX xyzdim = xdim*ydim*zdim;

  OC_REAL8m wgtx = -A/(mesh->EdgeLengthX()*mesh->EdgeLengthX());
  OC_REAL8m wgty = -A/(mesh->EdgeLengthY()*mesh->EdgeLengthY());
  OC_REAL8m wgtz = -A/(mesh->EdgeLengthZ()*mesh->EdgeLengthZ());
  
  Nb_Xpfloat energy_sum = 0.0;
  // KL(m) maxMdiffSq is now called maxdot
  // KL!?? I really do not see advantage of having three thread_maxdot_?
  // variables instead of one thread_maxdot_x
  // Maybe ask Mike?
  OC_REAL8m thread_maxdot = maxdot[threadnumber];
  OC_REAL8m thread_maxdot_x = thread_maxdot;
  OC_REAL8m thread_maxdot_y = thread_maxdot;
  OC_REAL8m thread_maxdot_z = thread_maxdot;
  // KL(m)
  OC_REAL8m thread_min_m = min_m[threadnumber];
  OC_REAL8m thread_max_m = max_m[threadnumber];
  // Note: For maxdot (and others) calculation, it suffices
  // to check spin[j]*spin[i] for j>i, or j<i, or various mixes of the
  // two.

  OC_INDEX x,y,z;
  mesh->GetCoords(node_start,x,y,z);

  OC_INDEX i = node_start;
  while(i<node_stop) {
    OC_INDEX xstop = xdim;
    if(xdim-x>node_stop-i) xstop = x + (node_stop-i);
    ThreeVector xlast(0.,0.,0.);
    if(x>0 || (x==0 && xperiodic)) {
      OC_INDEX j = i-1;
      if(x==0) j += xdim;
      if(Ms_inverse[j]!=0.0) xlast = (Ms[j]*spin[j]-Ms[i]*spin[i]); // KL(m)
    }
    while(x<xstop) {
      if(0 == Ms_inverse[i]) {
        if(ocedt.energy) (*ocedt.energy)[i] = 0.0;
        if(ocedt.H)      (*ocedt.H)[i].Set(0.,0.,0.);
        if(ocedt.mxH)    (*ocedt.mxH)[i].Set(0.,0.,0.);
        xlast.Set(0.,0.,0.);
      } else {
#if defined(__GNUC__) && __GNUC__ == 4 \
  && defined(__GNUC_MINOR__) && __GNUC_MINOR__ <= 1
        Oc_Nop(Ms_inverse[i]); // i686-apple-darwin9-gcc-4.0.1 (OS X Leopard)
        /// (others?) has an optimization bug that segfaults in the
        /// following code block.  Problem appears to be with the
        /// -fstrict-aliasing option (see demag-threaded.cc, line 1664),
        /// which is thrown in by -fast on Mac OS X/x86.
#endif
        // KL(m)
        OC_REAL8m m = Ms[i]/Me_T[i];
        if(m<thread_min_m) thread_min_m = m;
        if(m>thread_max_m) thread_max_m = m;
//      const OC_REAL8m hmult = (-2/MU0) * Ms_inverse[i];
        const OC_REAL8m hmult = (-2/MU0); // KL(m)

        const ThreeVector base = Ms[i]*spin[i]; // KL(m); also added const
        const OC_REAL8m   MeT2 = Me_T[i]*Me_T[i]; // KL(m)
        ThreeVector sum = xlast;
        if(x<xdim-1 || (x==xdim-1 && xperiodic)) {
          OC_INDEX j = i + 1;
          if(x==xdim-1) j -= xdim;
          if(Ms_inverse[j]!=0) {
            xlast = (base - Ms[j]*spin[j]); //KL(m)
            sum -= xlast;
            OC_REAL8m dot = xlast.MagSq()/MeT2 ; // KL(m), Mdiff normalization "on".
            if(dot>thread_maxdot_x) thread_maxdot_x = dot;
          }
        }
        sum *= wgtx/MeT2; //KL(m)
 
        ThreeVector tempy(0.,0.,0.);
        if(y>0 || (y==0 && yperiodic)) {
          OC_INDEX j = i-xdim;
          if(y==0) j += xydim;
          if(Ms_inverse[j]!=0.0) {
            tempy = (Ms[j]*spin[j] - base); // KL(m)
            OC_REAL8m dot = tempy.MagSq()/MeT2 ; // KL(m), Mdiff normalization "on".
            if(dot>thread_maxdot_y) thread_maxdot_y = dot;
          }
        }
        if(y<ydim-1 || (y==ydim-1 && yperiodic)) {
          OC_INDEX j = i+xdim;
          if(y==ydim-1) j -= xydim;
          if(Ms_inverse[j]!=0.0) {
            tempy += (Ms[j]*spin[j] - base); // KL(m)
          }
        }
        sum += tempy*(wgty/MeT2); //KL(m)

        ThreeVector tempz(0.,0.,0.);
        if(z>0 || (z==0 && zperiodic)) {
          OC_INDEX j = i-xydim;
          if(z==0) j+= xyzdim;
          if(Ms_inverse[j]!=0.0) {
            tempz = (Ms[j]*spin[j] - base); // KL(m)
            OC_REAL8m dot = tempz.MagSq()/MeT2 ; // KL(m), Mdiff normalization "on".
            if(dot>thread_maxdot_z) thread_maxdot_z = dot;
          }
        }
        if(z<zdim-1 || (z==zdim-1 && zperiodic)) {
          OC_INDEX j = i+xydim;
          if(z==zdim-1) j -= xyzdim;
          if(Ms_inverse[j]!=0.0) {
            tempz += (Ms[j]*spin[j] - base); // KL(m)
          }
        }
        sum += tempz*(wgtz/MeT2); // KL(m)

        OC_REAL8m ei = base.x*sum.x + base.y*sum.y + base.z*sum.z;
        sum.x *= hmult;  sum.y *= hmult;   sum.z *= hmult;

        // The following computation of mxH is wasted if ocedt.mxH and
        // ocedt.mxH_accum are NULL (as happens if one is using the
        // old-style GetEnergy interface), but timing tests imply that
        // the loss in that case is small, but the win is significant
        // for the case where ocedt.mxH_accum is non-NULL (as happens
        // if one is using the new-style
        // ComputeEnergy/AccumEnergyAndTorque interface).
        // KL(m)
        // sum= H, base= |M|.m
        // Thus, to get m^H you need bas^sum/|M|
        OC_REAL8m tx = (base.y*sum.z - base.z*sum.y)*Ms_inverse[i]; // KL(m)
        OC_REAL8m ty = (base.z*sum.x - base.x*sum.z)*Ms_inverse[i]; // KL(m)
        OC_REAL8m tz = (base.x*sum.y - base.y*sum.x)*Ms_inverse[i]; // KL(m)
        
        energy_sum += ei;
        if(ocedt.energy)       (*ocedt.energy)[i] = ei;
        if(ocedt.energy_accum) (*ocedt.energy_accum)[i] += ei;
        if(ocedt.H)            (*ocedt.H)[i] = sum;
        if(ocedt.H_accum)      (*ocedt.H_accum)[i] += sum;
        if(ocedt.mxH)          (*ocedt.mxH)[i] = ThreeVector(tx,ty,tz);
        if(ocedt.mxH_accum)    (*ocedt.mxH_accum)[i] += ThreeVector(tx,ty,tz);
      }
      ++i;
      ++x;
    }
    x=0;
    if((++y)>=ydim) {
      y=0;
      ++z;
    }
  }
  ocedtaux.energy_total_accum += energy_sum * mesh->Volume(0); 
  /// All cells have same volume in an Oxs_CommonRectangularMesh.

  if(thread_maxdot_x>thread_maxdot) thread_maxdot = thread_maxdot_x;
  if(thread_maxdot_y>thread_maxdot) thread_maxdot = thread_maxdot_y;
  if(thread_maxdot_z>thread_maxdot) thread_maxdot = thread_maxdot_z;
  maxdot[threadnumber] = thread_maxdot;
  
  // KL(m)
  min_m[threadnumber] = thread_min_m;
  max_m[threadnumber] = thread_max_m;
}

void Klm_UniformExchange::ComputeEnergyChunkInitialize
(const Oxs_SimState& state,
 Oxs_ComputeEnergyDataThreaded& ocedt,
 Oc_AlignedVector<Oxs_ComputeEnergyDataThreadedAux>& /* thread_ocedtaux */,
 int number_of_threads) const
{
  if(maxdot.size() != (vector<OC_REAL8m>::size_type)number_of_threads) {
    maxdot.resize(number_of_threads);
  }
  for(int i=0;i<number_of_threads;++i) {
    maxdot[i] = 0.0; // Minimum possible value for (m_i-m_j).MagSq()
  }
  
  // KL(m)
  if(min_m.size() != (vector<OC_REAL8m>::size_type)number_of_threads) {
    min_m.resize(number_of_threads);
  }
  if(max_m.size() != (vector<OC_REAL8m>::size_type)number_of_threads) {
    max_m.resize(number_of_threads);
  }
  for(int i=0;i<number_of_threads;++i) {
    min_m[i] = OC_REAL8m_MAX; // Maximum possible value for |m|
    max_m[i] = 0.0; // Minimum possible value for |m|
  }
  
  // (Re)-initialize coefficients if mesh has changed.  Note: Unlike
  // initialization code in some of the other Oxs_ChunkEnergy classes,
  // this initialization does not make calls back into the Tcl
  // interpreter.  Therefore, this initialization can be done by any one
  // thread; it doesn't have to be threadnumber == -1 (main thread).  So
  // we can manage with a simple mutex lock (as opposed to a more
  // complicated scheme involving ConditionWait calls.)
  const Oxs_CommonRectangularMesh* mesh
    = dynamic_cast<const Oxs_CommonRectangularMesh*>(state.mesh);
  if(mesh==NULL) {
    String msg=String("Object ")
      + String(state.mesh->InstanceName())
      + String(" is not a rectangular mesh.");
    throw Oxs_ExtError(this,msg);
  }
  if(mesh_id != mesh->Id()) {
    mesh_id = mesh->Id();
    //if(kernel == NGBR_12_ZD1B) {
    //  const OC_INDEX xdim = mesh->DimX();
    //  const OC_INDEX ydim = mesh->DimY();
    //  const OC_INDEX zdim = mesh->DimZ();
    //  if((1<xdim && xdim<5) || (1<ydim && ydim<5)
    //     || (1<zdim && zdim<5)) {
    //    char buf[1024];
    //    Oc_Snprintf(buf,sizeof(buf),
    //                "Each dimension must be ==1 or >=5 for 12ngbrZD1b kernel."
    //                " (Actual dimensions: xdim=%u, ydim=%u, zdim=%u.)",
    //                xdim,ydim,zdim);
    //    throw Oxs_ExtError(this,buf);
    //  }
    //  InitCoef_12NgbrZD1(xdim,xinteg,xcoef);
    //  InitCoef_12NgbrZD1(ydim,yinteg,ycoef);
    //  InitCoef_12NgbrZD1(zdim,zinteg,zcoef);
    //}
    OC_REAL8m minedge = OC_REAL8m_MAX;
    if(mesh->DimX() > 1) {
      minedge = mesh->EdgeLengthX();
    }
    if(mesh->DimY() > 1 && mesh->EdgeLengthY()<minedge) {
      minedge = mesh->EdgeLengthY();
    }
    if(mesh->DimZ() > 1 && mesh->EdgeLengthZ()<minedge) {
      minedge = mesh->EdgeLengthZ();
    }
    OC_REAL8m working_A = 0.0;
    if(excoeftype == A_TYPE) {
      working_A = A;
    } else if(excoeftype == LEX_TYPE) {
      OC_REAL8m lexMs = lex*state.max_absMs;
      if(lexMs>0) {
        working_A = 0.5*MU0*lexMs*lexMs;
      }
    } else {
      throw Oxs_ExtError(this,"Unsupported ExchangeCoefType.");
    }
    energy_density_error_estimate
      = 16*OC_REAL8m_EPSILON*working_A/minedge/minedge;
    // Worse case prefactor should be larger than 16, but in practice
    // error is probably smaller.
  }
  ocedt.energy_density_error_estimate = energy_density_error_estimate;  
}

void Klm_UniformExchange::ComputeEnergyChunkFinalize
(const Oxs_SimState& state,
 Oxs_ComputeEnergyDataThreaded& /* ocedt */,
 Oc_AlignedVector<Oxs_ComputeEnergyDataThreadedAux>&
    /* thread_ocedtaux */,
 int number_of_threads) const
{
  // KL(m) Min |M|/|Me|  and  Max |M|/|Me|   **************************
  OC_REAL8m total_min_m = OC_REAL8m_MAX;
  OC_REAL8m total_max_m = 0.0;
  for(int i=0;i<number_of_threads;++i) {
    if(min_m[i]<total_min_m) total_min_m = min_m[i];
    if(max_m[i]>total_max_m) total_max_m = max_m[i];
  }

  OC_REAL8m dummy_value_min_m;
  OC_REAL8m dummy_value_max_m;
  String min_m_name = MinMStateName();
  String max_m_name = MaxMStateName();
  if(state.GetDerivedData(min_m_name,dummy_value_min_m)) {
    // Here comes long original comment of Mike... ;)
#ifndef NDEBUG
    static Oxs_WarningMessage maxanglesetZ(3);
    maxanglesetZ.Send(revision_info,OC_STRINGIFY(__LINE__),
                     "Programming error?"
                     " Klm_UniformExchange min-m set twice.");
#endif
    // KL(m) clarifying it...
    // Appeared during -restart 1 - maybe only in such a case?
    // Past: switched off. I have no time now to investigate it :(
    //    if(1==0 && fabs(dummy_value_min_m-total_min_m)>OC_REAL8_EPSILON) {
    // Now: switched on. Maybe we will clarify it?
    if(fabs(dummy_value_min_m-total_min_m)>OC_REAL8_EPSILON) {
      char errbuf[1024];
      Oc_Snprintf(errbuf,sizeof(errbuf),
                  "Programming error (?):"
                  " Klm_UniformExchange min-m set to"
                  " two different values;"
                  " orig val=%#.17g, new val=%#.17g",
                  dummy_value_min_m,total_min_m);
      throw Oxs_ExtError(this,errbuf);
    }
  } else {
    state.AddDerivedData(min_m_name,total_min_m);
  }
  if(state.GetDerivedData(max_m_name,dummy_value_max_m)) {
#ifndef NDEBUG
    static Oxs_WarningMessage maxanglesetZ(3);
    maxanglesetZ.Send(revision_info,OC_STRINGIFY(__LINE__),
                     "Programming error?"
                     " Klm_UniformExchange max-m set twice.");
#endif
    // KL(m) clarifying it...
    // Appeared during -restart 1 - maybe only in such a case?
    // Past: switched off. I have no time now to investigate it :(
    //   if(1==0 && fabs(dummy_value_max_m-total_max_m)>OC_REAL8_EPSILON) {
    // Now: switched on. Maybe we will clarify it?
    if(fabs(dummy_value_max_m-total_max_m)>OC_REAL8_EPSILON) {
      char errbuf[1024];
      Oc_Snprintf(errbuf,sizeof(errbuf),
                  "Programming error (?):"
                  " Klm_UniformExchange max-m set to"
                  " two different values;"
                  " orig val=%#.17g, new val=%#.17g",
                  dummy_value_max_m,total_max_m);
      throw Oxs_ExtError(this,errbuf);
    }
  } else {
    state.AddDerivedData(max_m_name,total_max_m);
  }

  // Run and stage min/max data depend on data from the previous state.
  // In the case that the energy (and hence run/stage data)
  // for the current state was computed previously, then the previous
  // state may have been dropped.  So, compute and save run and stage
  // data iff they are not already computed.

  // Check stage and run min/max data from previous state
  const Oxs_SimState* oldstate = NULL;
  OC_REAL8m stage_min_m = OC_REAL8m_MAX;
  OC_REAL8m run_min_m = OC_REAL8m_MAX;
  OC_REAL8m stage_max_m = -1;
  OC_REAL8m run_max_m = -1;
  String s_min_m_name = StageMinMStateName();
  String r_min_m_name = RunMinMStateName();
  String s_max_m_name = StageMaxMStateName();
  String r_max_m_name = RunMaxMStateName();
  if(state.previous_state_id &&
     0 != (oldstate
      = director->FindExistingSimulationState(state.previous_state_id)) ) {
    if(oldstate->stage_number != state.stage_number) {
      stage_min_m = OC_REAL8m_MAX;
      stage_max_m = 0.0;
    } else {
      if(oldstate->GetDerivedData(s_min_m_name,dummy_value_min_m)) {
        stage_min_m = dummy_value_min_m;
      }
      if(oldstate->GetDerivedData(s_max_m_name,dummy_value_max_m)) {
        stage_max_m = dummy_value_max_m;
      }
    }
    if(oldstate->GetDerivedData(r_min_m_name,dummy_value_min_m)) {
      run_min_m = dummy_value_min_m;
    }
    if(oldstate->GetDerivedData(r_max_m_name,dummy_value_max_m)) {
      run_max_m = dummy_value_max_m;
    }
  }
  if(stage_min_m>total_min_m) stage_min_m = total_min_m;
  if(run_min_m>total_min_m)   run_min_m = total_min_m;
  if(stage_max_m<total_max_m) stage_max_m = total_max_m;
  if(run_max_m<total_max_m)   run_max_m = total_max_m;

  // Stage data
  if(!state.GetDerivedData(s_min_m_name,dummy_value_min_m)) {
    state.AddDerivedData(s_min_m_name,stage_min_m);
  }
  if(!state.GetDerivedData(s_max_m_name,dummy_value_max_m)) {
    state.AddDerivedData(s_max_m_name,stage_max_m);
  }

  // Run data
  if(!state.GetDerivedData(r_min_m_name,dummy_value_min_m)) {
    state.AddDerivedData(r_min_m_name,run_min_m);
  }
  if(!state.GetDerivedData(r_max_m_name,dummy_value_max_m)) {
    state.AddDerivedData(r_max_m_name,run_max_m);
  }
  
  // KL(m) Max Diff M, stored in maxdot ********************************************
  OC_REAL8m total_maxdot = 0.0;
  // Remember: we store in maxdot the squares ...
  for(int i=0;i<number_of_threads;++i) {
    if(maxdot[i]>total_maxdot) total_maxdot = maxdot[i];
  }
  // ... so, we have to correct it here 
  const OC_REAL8m total_maxMdiff = sqrt(total_maxdot);

  OC_REAL8m dummy_value_maxMdiff;
  String maxMdiff_name =  MaxMdiffStateName();
  if(state.GetDerivedData(maxMdiff_name,dummy_value_maxMdiff)) {
    // Ideally, energy values would never be computed more than once
    // for any one state, but in practice it seems inevitable that
    // such will occur on occasion.  For example, suppose a user
    // requests output on a state obtained by a stage crossing (as
    // opposed to a state obtained through a normal intrastage step);
    // a subsequent ::Step operation will re-compute the energies
    // because not all the information needed by the step transition
    // machinery is cached from an energy computation.  Even user
    // output requests on post ::Step states is problematic if some of
    // the requested output is not cached as part of the step
    // proceedings.  A warning is put into place below for debugging
    // purposes, but in general an error is raised only if results
    // from the recomputation are different than originally.
#ifndef NDEBUG
    static Oxs_WarningMessage maxanglesetZ(3);
    maxanglesetZ.Send(revision_info,OC_STRINGIFY(__LINE__),
                     "Programming error?"
                     " Klm_UniformExchange maxMdiff set twice.");
#endif
    // KL(m) clarifying it...
    // Appeared during -restart 1 - maybe only in such a case?
    // Past: switched off. I have no time now to investigate it :(
    //   if(1==0 && fabs(dummy_value_maxMdiff-total_maxMdiff)>OC_REAL8_EPSILON) {
    // Now: switched on. Maybe we will clarify it?
    if(fabs(dummy_value_maxMdiff-total_maxMdiff)>OC_REAL8_EPSILON) {
      char errbuf[1024];
      Oc_Snprintf(errbuf,sizeof(errbuf),
                  "Programming error (?):"
                  " Klm_UniformExchange maxMdiff set to"
                  " two different values;"
                  " orig val=%#.17g, new val=%#.17g",
                  dummy_value_maxMdiff,total_maxMdiff);
      throw Oxs_ExtError(this,errbuf);
    }
  } else {
    state.AddDerivedData(maxMdiff_name,total_maxMdiff);
  }

  // Run and stage maxMdiff data depend on data from the previous state.
  // In the case that the energy (and hence run/stage data)
  // for the current state was computed previously, then the previous
  // state may have been dropped.  So, compute and save run and stage
  // data iff they are not already computed.

  // Check stage and run maxMdiff data from previous state
  // const Oxs_SimState* oldstate = NULL; KL(m): defined already
  OC_REAL8m stage_maxMdiff = -1;
  OC_REAL8m run_maxMdiff = -1;
  String s_maxMdiff_name  = StageMaxMdiffStateName();
  String r_maxMdiff_name  = RunMaxMdiffStateName();
  if(state.previous_state_id &&
     0 != (oldstate
      = director->FindExistingSimulationState(state.previous_state_id)) ) {
    if(oldstate->stage_number != state.stage_number) {
      stage_maxMdiff = 0.0;
    } else {
      if(oldstate->GetDerivedData(s_maxMdiff_name,dummy_value_maxMdiff)) {
        stage_maxMdiff = dummy_value_maxMdiff;
      }
    }
    if(oldstate->GetDerivedData(r_maxMdiff_name,dummy_value_maxMdiff)) {
      run_maxMdiff = dummy_value_maxMdiff;
    }
  }
  if(stage_maxMdiff<total_maxMdiff) stage_maxMdiff = total_maxMdiff;
  if(run_maxMdiff<total_maxMdiff)   run_maxMdiff = total_maxMdiff;

  // Stage data
  if(!state.GetDerivedData(s_maxMdiff_name,dummy_value_maxMdiff)) {
    state.AddDerivedData(s_maxMdiff_name,stage_maxMdiff);
  }

  // Run data
  if(!state.GetDerivedData(r_maxMdiff_name,dummy_value_maxMdiff)) {
    state.AddDerivedData(r_maxMdiff_name,run_maxMdiff);
  }  
}

void Klm_UniformExchange::ComputeEnergyChunk
(const Oxs_SimState& state,
 Oxs_ComputeEnergyDataThreaded& ocedt,
 Oxs_ComputeEnergyDataThreadedAux& ocedtaux,
 OC_INDEX node_start,
 OC_INDEX node_stop,
 int threadnumber
 ) const
{
  const OC_INDEX size = state.mesh->Size();
  if(size<1) {
    return;
  }

  const Oxs_MeshValue<ThreeVector>& spin = state.spin;
  const Oxs_MeshValue<OC_REAL8m>& Ms = *(state.Ms);
  const Oxs_MeshValue<OC_REAL8m>& Ms_inverse = *(state.Ms_inverse);
  
  // KL(m)
  // We need sometimes access to equilibrium magnetisation length, Me(T)
  const Oxs_MeshValue<OC_REAL8m>& Me_T = *(state.GetPtr_Me_T());

  const Oxs_CommonRectangularMesh* mesh
    = dynamic_cast<const Oxs_CommonRectangularMesh*>(state.mesh);
  if(mesh==NULL) {
    String msg=String("Object ")
      + String(state.mesh->InstanceName())
      + String(" is not a rectangular mesh.");
    throw Oxs_ExtError(this,msg);
  }

  // Check periodicity.  Note that the following kernels have not been
  // upgraded to supported periodic meshes:
  //   NGBR_12_FREE, NGBR_12_ZD1, NGBR_12_ZD1B, NGBR_26
  // This is checked for and reported in the individual arms of the
  // kernel if-test below.
  const Oxs_RectangularMesh* rmesh 
    = dynamic_cast<const Oxs_RectangularMesh*>(mesh);
  const Oxs_PeriodicRectangularMesh* pmesh
    = dynamic_cast<const Oxs_PeriodicRectangularMesh*>(mesh);
  if(pmesh!=NULL) {
    // Rectangular, periodic mesh
    xperiodic = pmesh->IsPeriodicX();
    yperiodic = pmesh->IsPeriodicY();
    zperiodic = pmesh->IsPeriodicZ();
  } else if (rmesh!=NULL) {
    xperiodic=0; yperiodic=0; zperiodic=0;
  } else {
    String msg=String("Unknown mesh type: \"")
      + String(ClassName())
      + String("\".");
    throw Oxs_ExtError(this,msg.c_str());
  }
  
  // Note: Might want to consider subclassing exchange energies,
  //  to replace this if-block with virtual function pointers.
  if(kernel == NGBR_6_LLB) {
    CalcEnergy6NgbrMirror_LLB(spin,Ms_inverse,Ms,Me_T,mesh,ocedt,ocedtaux,
                       node_start,node_stop,threadnumber);
  } else {
    throw Oxs_ExtError(this,"Invalid kernel type detected."
                         "  (Programming error)");
  }

}

void Klm_UniformExchange::UpdateDerivedOutputs(const Oxs_SimState& state)
{
  // KL(m) ***************************************
  min_m_output.cache.state_id
    = stage_min_m_output.cache.state_id
    = run_min_m_output.cache.state_id
    //
    = max_m_output.cache.state_id
    = stage_max_m_output.cache.state_id
    = run_max_m_output.cache.state_id
    //
    = maxMdiff_output.cache.state_id
    = stage_maxMdiff_output.cache.state_id
    = run_maxMdiff_output.cache.state_id
    = 0;  // Mark change in progress

  String dummy_nameKL = MinMStateName();
  if(!state.GetDerivedData(dummy_nameKL.c_str(),
                           min_m_output.cache.value)) {
    // Error; This should always be set.  For now, just set the value to
    // -1, but in the future should consider throwing an exception.
    min_m_output.cache.value = -1.0;
  }

  dummy_nameKL = StageMinMStateName();
  if(!state.GetDerivedData(dummy_nameKL.c_str(),
                           stage_min_m_output.cache.value)) {
    // Error; This should always be set.  For now, just set the value to
    // -1, but in the future should consider throwing an exception.
    stage_min_m_output.cache.value = -1.0;
  }

  dummy_nameKL = RunMinMStateName();
  if(!state.GetDerivedData(dummy_nameKL.c_str(),
                           run_min_m_output.cache.value)) {
    // Error; This should always be set.  For now, just set the value to
    // -1, but in the future should consider throwing an exception.
    run_min_m_output.cache.value = -1.0;
  }
  //
  dummy_nameKL = MaxMStateName();
  if(!state.GetDerivedData(dummy_nameKL.c_str(),
                           max_m_output.cache.value)) {
    // Error; This should always be set.  For now, just set the value to
    // -1, but in the future should consider throwing an exception.
    max_m_output.cache.value = -1.0;
  }

  dummy_nameKL = StageMaxMStateName();
  if(!state.GetDerivedData(dummy_nameKL.c_str(),
                           stage_max_m_output.cache.value)) {
    // Error; This should always be set.  For now, just set the value to
    // -1, but in the future should consider throwing an exception.
    stage_max_m_output.cache.value = -1.0;
  }

  dummy_nameKL = RunMaxMStateName();
  if(!state.GetDerivedData(dummy_nameKL.c_str(),
                           run_max_m_output.cache.value)) {
    // Error; This should always be set.  For now, just set the value to
    // -1, but in the future should consider throwing an exception.
    run_max_m_output.cache.value = -1.0;
  }
  //
  dummy_nameKL = MaxMdiffStateName();
  if(!state.GetDerivedData(dummy_nameKL.c_str(),
                           maxMdiff_output.cache.value)) {
    // Error; This should always be set.  For now, just set the value to
    // -1, but in the future should consider throwing an exception.
    maxMdiff_output.cache.value = -1.0;
  }

  dummy_nameKL = StageMaxMdiffStateName();
  if(!state.GetDerivedData(dummy_nameKL.c_str(),
                           stage_maxMdiff_output.cache.value)) {
    // Error; This should always be set.  For now, just set the value to
    // -1, but in the future should consider throwing an exception.
    stage_maxMdiff_output.cache.value = -1.0;
  }

  dummy_nameKL = RunMaxMdiffStateName();
  if(!state.GetDerivedData(dummy_nameKL.c_str(),
                           run_maxMdiff_output.cache.value)) {
    // Error; This should always be set.  For now, just set the value to
    // -1, but in the future should consider throwing an exception.
    run_maxMdiff_output.cache.value = -1.0;
  }
    
  min_m_output.cache.state_id
    = stage_min_m_output.cache.state_id
    = run_min_m_output.cache.state_id
    //
    = max_m_output.cache.state_id
    = stage_max_m_output.cache.state_id
    = run_max_m_output.cache.state_id
    //
    = maxMdiff_output.cache.state_id
    = stage_maxMdiff_output.cache.state_id
    = run_maxMdiff_output.cache.state_id
    = state.Id(); 
  // KL(m) end *************************************  
}

// Optional interface for conjugate-gradient evolver.
// NOTE: At present, not all of the kernels support preconditioning.
// For details on this code, see NOTES VI, 21-July-2011, pp 10-11.
OC_INT4m
Klm_UniformExchange::IncrementPreconditioner(PreconditionerData& pcd)
{
  // Just using the parameter to avoid a compiler warning about an "unused variable"
  const Oxs_SimState& state = *(pcd.state);
        
  const String msg=String("IncrementPreconditioner not implemented!")
      + String(" (Unimportant string: ")
      + String(state.mesh->InstanceName())
      + String(".)");
  throw Oxs_ExtError(this,msg);	
	
  return 0; // No preconditioning support
}

