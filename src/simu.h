#pragma once

#include "mesh.h"
#include "fem.h"
#include <math.h>
#include <cmath>
#include <iostream>

namespace FEM2A {
    namespace Simu {

        //#################################
        //  Useful functions
        //#################################

        double unit_fct( vertex v )
        {
            return 1.;
        }

        double zero_fct( vertex v )
        {
            return 0.;
        }

        double xy_fct( vertex v )
        {
            return v.x + v.y;
        }

        //#################################
        //  Simulations
        //#################################

        void pure_dirichlet_pb( const std::string& mesh_filename, bool verbose )
        {
            std::cout << "Solving a pure Dirichlet problem \n" << std::endl;
            
            Mesh M;
            M.load(mesh_filename);
            M.set_attribute( unit_fct,  1, true);
            std::vector <double> F(M.nb_vertices());
            SparseMatrix K(M.nb_vertices());
			Quadrature quad;
			quad = quad.get_quadrature(0, false);
			ShapeFunctions shape = ShapeFunctions(2,1);
            
            for (int i = 0; i < M.nb_triangles(); ++i){
            	ElementMapping element( M, false, i);
            	DenseMatrix Ke;
            	Ke.set_size(3,3);
            	std::vector < double > Fe;
            	assemble_elementary_matrix( element, shape, quad, unit_fct, Ke);
            	local_to_global_matrix( M, i, Ke, K);
            	assemble_elementary_vector( element, shape, quad, unit_fct, Fe);
            	local_to_global_vector( M, false, i, Fe, F);
            }
            
            std::vector < bool > attr_dirich;
            std::vector < double > values;
            for (int i = 0; i < M.nb_edges(); ++i){
            	attr_dirich.push_back(true);
            	for (int j = 0; j < 2; ++j){
            		int index = M.get_edge_vertex_index(i, j); 
            		vertex vert = M.get_vertex( index ); 
            		values.push_back(vert.x + vert.y);
            	}
            }
        	apply_dirichlet_boundary_conditions( M, attr_dirich, values, K, F);
        	
        	std::vector <double> U;
        	solve(K, F, U);
        	//if ( verbose ) {
        		for (int i = 0; i < U.size(); ++i){ 
                	std::cout << U[i] << std::endl;
                }
            //}
            
        }

    }

}
