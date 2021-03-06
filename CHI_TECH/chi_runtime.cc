/** @file Runtime file*/
#include "chi_runtime.h"
#include <iostream>

#include "ChiConsole/chi_console.h"
#include "ChiMath/chi_math.h"
#include "ChiPhysics/chi_physics.h"
#include "ChiMesh/MeshHandler/chi_meshhandler.h"

#include <chi_mpi.h>
#include <chi_log.h>
#include "ChiTimer/chi_timer.h"

ChiConsole  ChiConsole::instance;
ChiMath     ChiMath::instance;
ChiMPI      ChiMPI::instance;
ChiLog      ChiLog::instance;
ChiPhysics ChiPhysics::instance;

ChiConsole&  chi_console = ChiConsole::GetInstance();
ChiMath&     chi_math_handler = ChiMath::GetInstance();
ChiMPI&      chi_mpi = ChiMPI::GetInstance();
ChiLog&      chi_log = ChiLog::GetInstance();
ChiPhysics&  chi_physics_handler = ChiPhysics::GetInstance();

/** Global stack of handlers */
std::vector<chi_mesh::MeshHandler*>  chi_meshhandler_stack;
int chi_current_mesh_handler=-1;

ChiTimer    chi_program_timer;

//=============================================== Global variables
bool            chi_termination_posted = false;
std::string     chi_input_file_name;
bool            chi_sim_option_interactive = true;




//############################################### Argument parser
/**Parses input arguments.
\param argc int    Number of arguments supplied.
\param argv char** Array of strings representing each argument.
 */
void ChiTech::ParseArguments(int argc, char** argv)
{
  bool input_file_found = false;
  for (int i=1; i<argc; i++)
  {
    std::string argument(argv[i]);

    if ((argument.find('=') == std::string::npos) && (!input_file_found) )
    {
      chi_input_file_name = argument;
      input_file_found = true;
      chi_sim_option_interactive = false;
    }//no =
    else if (argument.find('=') != std::string::npos)
    {
      luaL_dostring(chi_console.consoleState, argument.c_str());
    }//=
    //================================================ No-graphics option
    if ((argument.find("-b")!=std::string::npos)  )
    {
      chi_sim_option_interactive = false;
    }//-b
    //================================================ Verbosity
    if (argument.find("-v") != std::string::npos)
    {
      if ((i+1) >= argc)
      {
        std::cerr << "Invalid option used with command line argument "
                     "-v. Options are 0,1 or 2." << std::endl;
        exit(EXIT_FAILURE);
      }
      else
      {
        std::string v_option(argv[i+1]);
        try {
          int level = std::stoi(v_option);
          chi_log.SetVerbosity(level);
        }
        catch (const std::invalid_argument& e)
        {
          std::cerr << "Invalid option used with command line argument "
                       "-v. Options are 0,1 or 2." << std::endl;
          exit(EXIT_FAILURE);
        }
      }

    }//-v
  }//for argument
}

//############################################### Initialize ChiTech
/**Initializes all necessary items for ChiTech.
\param argc int    Number of arguments supplied.
\param argv char** Array of strings representing each argument.
 */
int ChiTech::Initialize(int argc, char** argv)
{
  ParseArguments(argc, argv);
  
  int location_id, number_processes;

  MPI_Init (&argc, &argv);                           /* starts MPI */
  MPI_Comm_rank (MPI_COMM_WORLD, &location_id);      /* get current process id */
  MPI_Comm_size (MPI_COMM_WORLD, &number_processes); /* get number of processes */

  chi_mpi.location_id = location_id;
  chi_mpi.process_count = number_processes;

  chi_console.PostMPIInfo(location_id, number_processes);
  chi_mpi.Initialize();

  chi_physics_handler.InitPetSc(argc,argv);

  return 0;
}

//############################################### Finalize ChiTech
/**Finalizes ChiTech.
 * */
void ChiTech::Finalize()
{
  PetscFinalize();
  MPI_Finalize();
}

//############################################### Interactive interface
/**Runs the interactive chitech engine*/
int ChiTech::RunInteractive(int argc, char** argv)
{
  chi_log.Log(LOG_0)
    << chi_program_timer.GetLocalDateTimeString()
    << " Running ChiTech in interactive-mode with "
    << chi_mpi.process_count << " processes.";

  chi_log.Log(LOG_0)
    << "ChiTech number of arguments supplied: "
    << argc - 1;

  if ( not chi_input_file_name.empty() )
    chi_console.ExecuteFile(chi_input_file_name.c_str(), argc, argv);

  chi_console.RunConsoleLoop();

  chi_log.Log(LOG_0)
    << "Final program time " << chi_program_timer.GetTimeString();

  chi_log.Log(LOG_0)
    << chi_program_timer.GetLocalDateTimeString()
    << " ChiTech finished execution.";

  return 0;
}

//############################################### Batch interface
/**Runs ChiTech in pure batch mode. Start then finish.*/
int ChiTech::RunBatch(int argc, char** argv)
{
  chi_log.Log(LOG_0)
    << chi_program_timer.GetLocalDateTimeString()
    << " Running ChiTech in batch-mode with "
    << chi_mpi.process_count << " processes.";

  chi_log.Log(LOG_0)
    << "ChiTech number of arguments supplied: "
    << argc - 1;

  if (argc<=1)
    chi_log.Log(LOG_0)
      << "\nUsage: exe inputfile [options values]\n"
      << "\n"
      << "     -v    Level of verbosity. Default 0. Can be either 0, 1 or 2.\n"
      << "     a=b   Executes argument as a lua string.\n\n\n";

  int error_code = 0;
  if ( not chi_input_file_name.empty() )
    error_code = chi_console.ExecuteFile(chi_input_file_name.c_str(), argc, argv);

  chi_log.Log(LOG_0)
    << "Final program time " << chi_program_timer.GetTimeString();

  chi_log.Log(LOG_0)
    << chi_program_timer.GetLocalDateTimeString()
    << " ChiTech finished execution of " << chi_input_file_name;

  return error_code;
}


