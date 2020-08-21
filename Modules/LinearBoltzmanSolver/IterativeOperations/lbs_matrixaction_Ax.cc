#include "lbs_matrixaction_Ax.h"
#include "../Tools/ksp_data_context.h"
#include "ChiMesh/SweepUtilities/SweepScheduler/sweepscheduler.h"

#include "../../DiffusionSolver/Solver/diffusion_solver.h"

typedef chi_mesh::sweep_management::SweepScheduler MainSweepScheduler;
//###################################################################
/**Computes the action of the transport matrix on a vector.*/
int NPTMatrixAction_Ax(Mat matrix, Vec krylov_vector, Vec Ax)
{
  KSPDataContext* context;
  MatShellGetContext(matrix,&context);

  LinearBoltzman::Solver* solver = context->solver;
  SweepChunk* sweep_chunk = context->sweep_chunk;
  LBSGroupset* groupset  = context->groupset;
  MainSweepScheduler* sweepScheduler = context->sweepScheduler;

  //============================================= Copy krylov vector into local
  solver->DisAssembleVector(groupset,
                            krylov_vector,
                            solver->phi_old_local.data());

  //============================================= Setting the source using
  //                                             updated phi_old
  solver->SetSource(context->group_set_num,
                    LinearBoltzman::SourceFlags::USE_DLINV_SOURCE);

  //============================================= Sweeping the new source
  groupset->angle_agg->ResetDelayedPsi();

  solver->phi_new_local.assign(solver->phi_new_local.size(),0.0);
  sweepScheduler->Sweep(sweep_chunk);

  //=================================================== Apply WGDSA
  if (groupset->apply_wgdsa)
  {
    solver->AssembleWGDSADeltaPhiVector(groupset,
                                        solver->phi_old_local.data(),
                                        solver->phi_new_local.data());
    ((chi_diffusion::Solver*)groupset->wgdsa_solver)->ExecuteS(true,false);
    solver->DisAssembleWGDSADeltaPhiVector(groupset,
                                           solver->phi_new_local.data());
  }
  if (groupset->apply_tgdsa)
  {
    solver->AssembleTGDSADeltaPhiVector(groupset,
                                        solver->phi_old_local.data(),
                                        solver->phi_new_local.data());
    ((chi_diffusion::Solver*)groupset->tgdsa_solver)->ExecuteS(true,false);
    solver->DisAssembleTGDSADeltaPhiVector(groupset,
                                           solver->phi_new_local.data());
  }


  solver->AssembleVector(groupset,
                         context->x_temp,
                         solver->phi_new_local.data());


  //============================================= Computing action
  VecWAXPY(Ax,-1.0,context->x_temp,krylov_vector);

  return 0;
}

//###################################################################
/**Computes the action of the transport matrix on a vector.*/
int NPTMatrixAction_Ax_POD_b(Mat matrix, Vec input_vector, Vec Ax)
{
  KSPDataContext* context;
  MatShellGetContext(matrix,&context);

  LinearBoltzman::Solver* solver = context->solver;
  SweepChunk* sweep_chunk = context->sweep_chunk;
  LBSGroupset* groupset  = context->groupset;
  MainSweepScheduler* sweepScheduler = context->sweepScheduler;

  //============================================= Copy krylov vector into local
  solver->DisAssembleVector(groupset,
                            input_vector,
                            solver->phi_old_local.data());

  //============================================= Setting the source using
  //                                             updated phi_old
  solver->SetSource(context->group_set_num,
                    false, // Apply material src???
                    true); // Suppress MSphi ???

  //============================================= Sweeping the new source
  groupset->angle_agg->ResetDelayedPsi();

  solver->phi_new_local.assign(solver->phi_new_local.size(),0.0);
  sweepScheduler->Sweep(sweep_chunk);

  solver->AssembleVector(groupset,
                         context->x_temp,
                         solver->phi_new_local.data());


  //============================================= Computing action
  VecWAXPY(Ax,-1.0,context->x_temp,input_vector);

  return 0;
}

//###################################################################
/**Computes the action of the transport matrix on a vector.*/
int NPTMatrixAction_Ax_POD(Mat matrix, Vec input_vector, Vec Ax)
{
  KSPDataContext* context;
  MatShellGetContext(matrix,&context);

  LinearBoltzman::Solver* solver = context->solver;
  SweepChunk* sweep_chunk = context->sweep_chunk;
  LBSGroupset* groupset  = context->groupset;
  MainSweepScheduler* sweepScheduler = context->sweepScheduler;

  //============================================= Copy krylov vector into local
  solver->DisAssembleVector(groupset,
                            input_vector,
                            solver->phi_old_local.data());

  //============================================= Setting the source using
  //                                             updated phi_old
  solver->SetSource(context->group_set_num,
                    false, // Apply material src? It is MSphi in this case
                    true); // Dont iterate to include scattering

  //============================================= Sweeping the new source
  groupset->angle_agg->ResetDelayedPsi();

  solver->phi_new_local.assign(solver->phi_new_local.size(),0.0);
  sweepScheduler->Sweep(sweep_chunk);

  solver->AssembleVector(groupset,
                         context->x_temp,
                         solver->phi_new_local.data());


  //============================================= Computing action
  VecWAXPY(Ax,-1.0,context->x_temp,input_vector);

  return 0;
}