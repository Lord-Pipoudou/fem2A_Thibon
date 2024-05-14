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
        
        double sin_bump( vertex v )
        {
        	return 2*pow(2, M_PI)*sin(M_PI*v.x)*sin(M_PI*v.y);
        }
        
        double neum( vertex v )
        {
        	return sin(M_PI*v.y);
        }
        
        double droite( vertex v )
        {
        	if (abs(v.x - 1) < 0.001){
        		return 1;
        	}
        	else {
	        	return 0;
	        }
        }
        
        double gauche( vertex v )
        {
        	if (abs(v.x) < 0.001){
        		return 1;
        	}
        	else {
	        	return 0;
	        }
        }
        
        
        

        void pure_dirichlet_pb( const std::string& mesh_filename, bool verbose )
        {
            std::cout << "Solving a pure Dirichlet problem \n" << std::endl;
            
            Mesh M;
            M.load(mesh_filename);
            M.set_attribute( unit_fct,  0, true);
            std::vector <double> F(M.nb_vertices(), 0);
            SparseMatrix K(M.nb_vertices());
			Quadrature quad;
			quad = quad.get_quadrature(2, false);
			ShapeFunctions shape = ShapeFunctions(2,1);
            
            //sur les triangles
            for (int i = 0; i < M.nb_triangles(); ++i){
            	ElementMapping element( M, false, i);
            	DenseMatrix Ke;
            	Ke.set_size(3,3);
            	std::vector < double > Fe;
            	assemble_elementary_matrix( element, shape, quad, unit_fct, Ke);
            	local_to_global_matrix( M, i, Ke, K);
            	/*
            	assemble_elementary_vector( element, shape, quad, unit_fct, Fe);
            	local_to_global_vector( M, false, i, Fe, F);*/
            }
            
            std::vector < bool > attr_dirich;
            std::vector < double > values(M.nb_vertices());
            attr_dirich.push_back(true);
            for (int i = 0; i < M.nb_vertices(); ++i){
            	values[i] = xy_fct(M.get_vertex(i));
            }
        	apply_dirichlet_boundary_conditions( M, attr_dirich, values, K, F);
        	std::vector <double> U;
        	solve(K, F, U);
        	save_solution(U, "pure_dirichlet_square_fine.bb");
        	if ( verbose ) {
        		for (int i = 0; i < U.size(); ++i){ 
                	std::cout << U[i] << std::endl;
                }
            }
            
        }
        
        void source_dirichlet_pb( const std::string& mesh_filename, bool verbose )
        {
            std::cout << "Solving a source Dirichlet problem \n" << std::endl;
            
            Mesh M;
            M.load(mesh_filename);
            M.set_attribute( unit_fct,  0, true);
            std::vector <double> F(M.nb_vertices(), 0);
            SparseMatrix K(M.nb_vertices());
			Quadrature quad;
			quad = quad.get_quadrature(2, false);
			ShapeFunctions shape = ShapeFunctions(2,1);
            
            //sur les triangles
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
            std::vector < double > values(M.nb_vertices());
            attr_dirich.push_back(true);
            for (int i = 0; i < M.nb_vertices(); ++i){
            	values[i] = zero_fct(M.get_vertex(i));
            }
        	apply_dirichlet_boundary_conditions( M, attr_dirich, values, K, F);
        	std::vector <double> U;
        	solve(K, F, U);
        	save_solution(U, "source_dirichlet_square_fine.bb");
        	if ( verbose ) {
        		for (int i = 0; i < U.size(); ++i){ 
                	std::cout << U[i] << std::endl;
                }
            }
            
        }
        
        void sinus_bump_dirichlet_pb( const std::string& mesh_filename, bool verbose )
        {
            std::cout << "Solving a sinus bump Dirichlet problem \n" << std::endl;
            
            Mesh M;
            M.load(mesh_filename);
            M.set_attribute( unit_fct,  0, true);
            std::vector <double> F(M.nb_vertices(), 0);
            SparseMatrix K(M.nb_vertices());
			Quadrature quad;
			quad = quad.get_quadrature(2, false);
			ShapeFunctions shape = ShapeFunctions(2,1);
            
            //sur les triangles
            for (int i = 0; i < M.nb_triangles(); ++i){
            	ElementMapping element( M, false, i);
            	DenseMatrix Ke;
            	Ke.set_size(3,3);
            	std::vector < double > Fe;
            	assemble_elementary_matrix( element, shape, quad, unit_fct, Ke);
            	local_to_global_matrix( M, i, Ke, K);
            	
            	assemble_elementary_vector( element, shape, quad, sin_bump, Fe);
            	local_to_global_vector( M, false, i, Fe, F);
            }
            
            std::vector < bool > attr_dirich;
            std::vector < double > values(M.nb_vertices());
            attr_dirich.push_back(true);
            for (int i = 0; i < M.nb_vertices(); ++i){
            	values[i] = zero_fct(M.get_vertex(i));
            }
        	apply_dirichlet_boundary_conditions( M, attr_dirich, values, K, F);
        	std::vector <double> U;
        	solve(K, F, U);
        	save_solution(U, "sinus_bump_square.bb");
        	if ( verbose ) {
        		for (int i = 0; i < U.size(); ++i){ 
                	std::cout << U[i] << std::endl;
                }
            }
            
        }
        
        void neumann_pb( const std::string& mesh_filename, bool verbose )
        {
            std::cout << "Solving a neumann problem \n" << std::endl;
            
            Mesh M;
            M.load(mesh_filename);
            M.set_attribute( unit_fct,  2, true);
            M.set_attribute( droite,  0, true);
            M.set_attribute( gauche,  1, true);

            std::vector <double> F(M.nb_vertices(), 0);
            SparseMatrix K(M.nb_vertices());
			Quadrature quad;
			Quadrature quad_1D;
			quad = quad.get_quadrature(2, false);
			ShapeFunctions shape = ShapeFunctions(2,1);
			quad_1D = quad_1D.get_quadrature(2, true);
			ShapeFunctions shape_1D = ShapeFunctions(1,1);
            
            //sur les triangles
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
            std::vector < double > values(M.nb_vertices());
            attr_dirich.push_back(true);
            attr_dirich.push_back(false);
            attr_dirich.push_back(false);
            for (int i = 0; i < M.nb_vertices(); ++i){
            	values[i] = zero_fct(M.get_vertex(i));
            }
            for (int i = 0; i < M.nb_edges(); ++i){
				if( M.get_edge_attribute(i) == 1 ){
					ElementMapping element( M, true, i);
					std::vector < double > Fe;
					assemble_elementary_neumann_vector(element, shape_1D, quad_1D, neum, Fe);
					local_to_global_vector( M, true, i, Fe, F);
				}
			}
			
			apply_dirichlet_boundary_conditions( M, attr_dirich, values, K, F);
			std::vector <double> U;
			solve(K, F, U);
			save_solution(U, "neumann_square_fine.bb");
			if ( verbose ) {
				for (int i = 0; i < U.size(); ++i){ 
					std::cout << U[i] << std::endl;
				}
            }
            
        }
        

    }

}
