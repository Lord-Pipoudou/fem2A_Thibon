#include <iostream>
#include <string>
#include <vector>

#include "src/fem.h"
#include "src/mesh.h"
#include "src/solver.h"
#include "src/tests.h"
#include "src/simu.h"
/*
/* Global variables */
std::vector< std::string > arguments;

/* To parse command line arguments */
bool flag_is_used(
    const std::string& flag,
    const std::vector< std::string >& arguments )
{
    for( int i = 0; i < arguments.size(); ++i ) {
        if( flag == arguments[i] ) {
            return true;
        }
    }
    return false;
}

using namespace FEM2A;

void run_tests()
{
    const bool t_opennl = false;
    const bool t_lmesh = false;
    const bool t_io = false;
    const bool t_quad = false;
    const bool t_edge = false;
    const bool t_triangle = false;
    const bool t_jac_edge = false;
    const bool t_jac_triangle = false;
    const bool t_assemble_matrix = false;
    const bool t_assemble_vector = false;
    const bool t_K = true;
    const bool t_F = true;

    if( t_opennl ) test_opennl();
    if( t_lmesh ) Tests::test_load_mesh();
    if( t_io ) Tests::test_load_save_mesh();
    if( t_quad ) Tests::test_quadrature();
    if( t_edge ) Tests::test_element_edge();
    if( t_triangle ) Tests::test_element_triangle();
    if( t_jac_edge ) Tests::test_jacob_edge();
    if( t_jac_triangle ) Tests::test_jacob_triangle();
    if( t_assemble_matrix) Tests::test_assemble_elementary_matrix();
    if( t_assemble_vector) Tests::test_assemble_elementary_vector();
    if( t_K) Tests::test_assemble_Ke_K();
    if( t_F ) Tests::test_F();
}


void run_simu()
{

    const bool simu_pure_dirichlet = false;
    const bool simu_source_dirichlet = false;
    const bool simu_sinus_bump_dirichlet = false;
    const bool simu_neumann = false;
    const bool simu_mug = true;

    const bool verbose = flag_is_used( "-v", arguments )
        || flag_is_used( "--verbose", arguments );

    if( simu_pure_dirichlet ) {
        Simu::pure_dirichlet_pb("data/square_fine.mesh", verbose);
    }
    if( simu_source_dirichlet ) {
        Simu::source_dirichlet_pb("data/square_fine.mesh", verbose);
    }
    if( simu_sinus_bump_dirichlet ) {
        Simu::sinus_bump_dirichlet_pb("data/square.mesh", verbose);
    }
    
    if( simu_neumann ) {
        Simu::neumann_pb("data/square_fine.mesh", verbose);
    }
    
    if( simu_mug ) {
        Simu::mug_pb("data/mug_1.mesh", verbose);
    }
    
}

int main( int argc, const char * argv[] )
{
    /* Command line parsing */
    for( int i = 1; i < argc; ++i ) {
        arguments.push_back( std::string(argv[i]) );
    }

    /* Show usage if asked or no arguments */
    if( arguments.size() == 0 || flag_is_used("-h", arguments)
        || flag_is_used("--help", arguments) ) {
        std::cout << "Usage: ./fem2a [options]" << std::endl
            << "Options: " << std::endl;
        std::cout << " -h, --help:        show usage" << std::endl;
        std::cout << " -t, --run-tests:   run the tests" << std::endl;
        std::cout << " -s, --run-simu:    run the simulations" << std::endl;
        std::cout << " -v, --verbose:     print lots of details" << std::endl;
        return 0;
    }

    /* Run the tests if asked */
    if( flag_is_used("-t", arguments)
        || flag_is_used("--run-tests", arguments) ) {
        run_tests();
    }

    /* Run the simulation if asked */
    if( flag_is_used("-s", arguments)
        || flag_is_used("--run-simu", arguments) ) {
        run_simu();
    }

    return 0;
}

