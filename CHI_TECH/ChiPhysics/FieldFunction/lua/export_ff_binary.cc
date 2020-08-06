#include "ChiLua/chi_lua.h"

#include "ChiPhysics/chi_physics.h"
#include "ChiPhysics/FieldFunction/fieldfunction.h"

#include "ChiMath/UnknownManager/unknown_manager.h"
#include "ChiMath/SpatialDiscretization/PiecewiseLinear/pwl.h"

#include "ChiMesh/MeshContinuum/chi_meshcontinuum.h"

#include "chi_log.h"
extern ChiLog& chi_log;

#include "chi_mpi.h"
extern ChiMPI& chi_mpi;

#include <fstream>

//###################################################################
/**This is a lua test function.
\param argument1 Any Argument of any type.
\ingroup LuaGeneralUtilities
 */
int chiExportFieldFunctionToBinary(lua_State* L)
{
  auto& physics_handler = ChiPhysics::GetInstance();

  int num_args = lua_gettop(L);
  if (num_args != 2)
    LuaPostArgAmountError("chiExportFieldFunctionToBinary",2,num_args);

  int handle =  lua_tonumber(L,1);
  const char* file_name = lua_tostring(L,2);

  auto ff = physics_handler.fieldfunc_stack[handle];
  auto uk_man = ff->unknown_manager;
  auto pwl = (SpatialDiscretization_PWL*)ff->spatial_discretization;

  std::ofstream file;
  file.open(std::string(file_name)+"_"+
            std::to_string(chi_mpi.location_id)+
            std::string(".txt"));

  int M = uk_man->unknowns.size();
  int G = uk_man->unknowns.front()->num_components;

  file << M << "\n";
  file << G << "\n";

  auto grid = pwl->ref_grid;

  for (auto& cell : grid->local_cells)
  {
    file << cell.global_id << "\n";
    
    for (int v=0; v<cell.vertex_ids.size(); ++v)
    {
      for (int m=0; m<M; ++m)
      {
        for (int g=0; g<G; ++g)
        {
          int ig = pwl->MapDFEMDOF(&cell,v,uk_man,0,g);
          int il = pwl->MapDFEMDOFLocal(&cell,v,uk_man,0,g);
          file << ig << " " << (*ff->field_vector_local)[il] << "\n";
        }//for g
      }//for m
    }//for v
  }//for cell

  file.close();


  return 0;
}