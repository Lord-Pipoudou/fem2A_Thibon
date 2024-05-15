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
        	return 2*pow(M_PI, 2)*sin(M_PI*v.x)*sin(M_PI*v.y);
        }
        
        double soluce_sin_bump( vertex v )
        {
        	return sin(M_PI*v.x)*sin(M_PI*v.y);
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
        
        double mug_eau_bouill( vertex v )
        {
        	double eps = 0.001; 
        	
        	// Pour le cas ou le vertex se situe sur le fond du mug
        	if (abs(v.y - 1) < eps){
        	
        		if (v.x > 1 - eps && v.x < 20 + eps){
        			return 1;
        		}
        		
        		else{ return 0; }
        	}
        	
        	// Pour le cas où le vertex se situe sur la paroi verticale intérieure en contact avec l'eau bouillante
        	else {
        		if (abs(v.x - 1) < eps || abs(v.x - 20) < eps){
        		
        			if (v.y > 1 - eps && v.y < 10 + eps){
        				return 1;
        			}
        			
        			else{ return 0; }
        		}
        		
        		else{ return 0; }
        	}
        }
        
        double mug_neum( vertex v )
        {
        	return -0.1;
        }
        
        
       // -------------------------------------------------------------------- Dirichlet Pur

        void pure_dirichlet_pb( const std::string& mesh_filename, bool verbose )
        {
            std::cout << "Solving a pure Dirichlet problem \n" << std::endl;
            
            // CRÉATION DU MAILLAGE
            Mesh M;
            M.load(mesh_filename);
            M.set_attribute( unit_fct,  0, true); // on applique à tout le bord les conditions de dirichlet
            
            // INITIALISATION DE LA MATRICE K ET DU VECTEUR F (inutile ici)
            std::vector <double> F(M.nb_vertices(), 0);
            SparseMatrix K(M.nb_vertices());
            
            // CRÉATION DE LA QUADRATURE D'ORDRE 2 ET DES FONCTIONS DE FORME
			Quadrature quad;
			quad = quad.get_quadrature(2, false);
			ShapeFunctions shape = ShapeFunctions(2,1);
            
            //sur chaque triangle
            for (int i = 0; i < M.nb_triangles(); ++i){
            	ElementMapping element( M, false, i);
            	DenseMatrix Ke;
            	Ke.set_size(3,3);
            	std::vector < double > Fe;
            	
            	// CRÉATION DE LA MATRICE Ke ET SON IMPLÉMENTATION DANS K
            	assemble_elementary_matrix( element, shape, quad, unit_fct, Ke);
            	local_to_global_matrix( M, i, Ke, K);
            }
            
            // CONFIGURATION DE attr_dirich (0 --> dirichlet)
            std::vector < bool > attr_dirich;
            std::vector < double > values(M.nb_vertices());
            attr_dirich.push_back(true);
            
            // APPLICATION DES CONDITIONS DE DIRICHLET
            for (int i = 0; i < M.nb_vertices(); ++i){
            	values[i] = xy_fct(M.get_vertex(i));
            }
            
        	apply_dirichlet_boundary_conditions( M, attr_dirich, values, K, F);
        	
        	// CALCUL DE LA SOLUTION
        	std::vector <double> U;
        	solve(K, F, U);
        	save_solution(U, "pure_dirichlet_square_fine.bb");
        	
        	// AFFICHE LA SOLUTION DANS LE TERMINAL SI -v
        	if ( verbose ) {
        		std::cout << "Affichage de la solution :" << std::endl;
        		for (int i = 0; i < U.size(); ++i){ 
                	std::cout << U[i] << std::endl;
                }
            }
            
        }
        
        // -------------------------------------------------------------------- Dirichlet Source
        
        void source_dirichlet_pb( const std::string& mesh_filename, bool verbose )
        {
            std::cout << "Solving a source Dirichlet problem \n" << std::endl;
            
            // CRÉATION DU MAILLAGE
            Mesh M;
            M.load(mesh_filename);
            M.set_attribute( unit_fct,  0, true); // on applique à tout le bord les conditions de dirichlet
            
            // INITIALISATION DE LA MATRICE K ET DU VECTEUR F
            std::vector <double> F(M.nb_vertices(), 0);
            SparseMatrix K(M.nb_vertices());
            
            // CRÉATION DE LA QUADRATURE D'ORDRE 2 ET DES FONCTIONS DE FORME
			Quadrature quad;
			quad = quad.get_quadrature(2, false);
			ShapeFunctions shape = ShapeFunctions(2,1);
            
            //sur chaque triangle
            for (int i = 0; i < M.nb_triangles(); ++i){
            	ElementMapping element( M, false, i);
            	DenseMatrix Ke;
            	Ke.set_size(3,3);
            	std::vector < double > Fe;
            	
            	// CRÉATION DE LA MATRICE Ke ET SON IMPLÉMENTATION DANS K
            	assemble_elementary_matrix( element, shape, quad, unit_fct, Ke);
            	local_to_global_matrix( M, i, Ke, K);
            	
            	// CRÉATION DU VECTEUR Fe ET SON IMPLÉMENTATION DANS F (terme source)
            	assemble_elementary_vector( element, shape, quad, unit_fct, Fe);
            	local_to_global_vector( M, false, i, Fe, F);
            }
            
            // CONFIGURATION DE attr_dirich (0 --> dirichlet)
            std::vector < bool > attr_dirich;
            std::vector < double > values(M.nb_vertices());
            attr_dirich.push_back(true);
            
            // APPLICATION DES CONDITIONS DE DIRICHLET
            for (int i = 0; i < M.nb_vertices(); ++i){
            	values[i] = zero_fct(M.get_vertex(i));
            }
            
        	apply_dirichlet_boundary_conditions( M, attr_dirich, values, K, F);
        	
        	// CALCUL DE LA SOLUTION
        	std::vector <double> U;
        	solve(K, F, U);
        	save_solution(U, "source_dirichlet_square_fine.bb");
        	
        	// AFFICHE LA SOLUTION DANS LE TERMINAL SI -v
        	if ( verbose ) {
        		std::cout << "Affichage de la solution :" << std::endl;
        		for (int i = 0; i < U.size(); ++i){ 
                	std::cout << U[i] << std::endl;
                }
            }
            
        }
        
        // -------------------------------------------------------------------- Sinus Bump
        
        void sinus_bump_dirichlet_pb( const std::string& mesh_filename, bool verbose )
        {
            std::cout << "Solving a sinus bump Dirichlet problem \n" << std::endl;
            bool diff = true; // Fait la différence entre la solution analytique et la solution obtenue par la méthode des éléments finis
            
            // CRÉATION DU MAILLAGE
            Mesh M;
            M.load(mesh_filename);
            M.set_attribute( unit_fct,  0, true); // on applique à tout le bord les conditions de dirichlet
            
            // INITIALISATION DE LA MATRICE K ET DU VECTEUR F
            std::vector <double> F(M.nb_vertices(), 0);
            SparseMatrix K(M.nb_vertices());
            
            // CRÉATION DE LA QUADRATURE D'ORDRE 2 ET DES FONCTIONS DE FORME
			Quadrature quad;
			quad = quad.get_quadrature(2, false);
			ShapeFunctions shape = ShapeFunctions(2,1);
            
            //sur chaque triangle
            for (int i = 0; i < M.nb_triangles(); ++i){
            	ElementMapping element( M, false, i);
            	DenseMatrix Ke;
            	Ke.set_size(3,3);
            	std::vector < double > Fe;
            	
            	// CRÉATION DE LA MATRICE Ke ET SON IMPLÉMENTATION DANS K
            	assemble_elementary_matrix( element, shape, quad, unit_fct, Ke);
            	local_to_global_matrix( M, i, Ke, K);
            	
            	// CRÉATION DU VECTEUR Fe ET SON IMPLÉMENTATION DANS F (terme source)
            	assemble_elementary_vector( element, shape, quad, sin_bump, Fe);
            	local_to_global_vector( M, false, i, Fe, F);
            }
            
            // CONFIGURATION DE attr_dirich (0 --> dirichlet)
            std::vector < bool > attr_dirich;
            std::vector < double > values(M.nb_vertices());
            attr_dirich.push_back(true);
            
            // APPLICATION DES CONDITIONS DE DIRICHLET
            for (int i = 0; i < M.nb_vertices(); ++i){
            	values[i] = zero_fct(M.get_vertex(i));
            }
            
        	apply_dirichlet_boundary_conditions( M, attr_dirich, values, K, F);
        	
        	// CALCUL DE LA SOLUTION
        	std::vector <double> U;
        	solve(K, F, U);
        	
        	// CALCUL DE LA DIFFERENCE ENTRE LA SOLUTION ANALYTIQUE ET CELLE TROUVÉE (si diff = true)
        	if( diff ){
        	  	for (int i = 0; i < M.nb_vertices(); ++i){
            		U[i] -= soluce_sin_bump( M.get_vertex(i) );
            	}
            }
            
        	save_solution(U, "sinus_bump_diff_square.bb");
        	
        	// AFFICHE LA SOLUTION DANS LE TERMINAL SI -v
        	if ( verbose ) {
        		std::cout << "Affichage de la solution :" << std::endl;
        		for (int i = 0; i < U.size(); ++i){ 
                	std::cout << U[i] << std::endl;
                }
            }
            
            
        }
        
        // -------------------------------------------------------------------- Neumann
        
        void neumann_pb( const std::string& mesh_filename, bool verbose )
        {
            std::cout << "Solving a neumann problem \n" << std::endl;
            
            // CRÉATION DU MAILLAGE
            Mesh M;
            M.load(mesh_filename);
            M.set_attribute( unit_fct,  2, true); 
            M.set_attribute( droite,  0, true); // on applique au bord de droite les conditions de dirichlet
            M.set_attribute( gauche,  1, true); // on applique au bord de gauche les conditions de neumann

			// INITIALISATION DE LA MATRICE K ET DU VECTEUR F
            std::vector <double> F(M.nb_vertices(), 0);
            SparseMatrix K(M.nb_vertices());
            
            // INITIALISATION DES QUADRATURES 1D ET 2D
			Quadrature quad;
			Quadrature quad_1D;
			
			// CRÉATION DE LA QUADRATURE D'ORDRE 2 ET DES FONCTIONS DE FORME
			quad = quad.get_quadrature(2, false);
			ShapeFunctions shape = ShapeFunctions(2,1);
			
			// CRÉATION DE LA QUADRATURE 1D D'ORDRE 2 ET DES FONCTIONS DE FORME
			quad_1D = quad_1D.get_quadrature(2, true);
			ShapeFunctions shape_1D = ShapeFunctions(1,1);
            
            //sur chaque triangle
            for (int i = 0; i < M.nb_triangles(); ++i){
            	ElementMapping element( M, false, i);
            	DenseMatrix Ke;
            	Ke.set_size(3,3);
            	std::vector < double > Fe;
            	
            	// CRÉATION DE LA MATRICE Ke ET SON IMPLÉMENTATION DANS K
            	assemble_elementary_matrix( element, shape, quad, unit_fct, Ke);
            	local_to_global_matrix( M, i, Ke, K);
            	
            	// CRÉATION DU VECTEUR Fe ET SON IMPLÉMENTATION DANS F (terme source)
            	assemble_elementary_vector( element, shape, quad, unit_fct, Fe);
            	local_to_global_vector( M, false, i, Fe, F);
            }
            
            // CONFIGURATION DE attr_dirich (0 --> dirichlet, 1 --> neumann, 2 --> neumann nul)
            std::vector < bool > attr_dirich;
            std::vector < double > values(M.nb_vertices());
            attr_dirich.push_back(true);
            attr_dirich.push_back(false);
            attr_dirich.push_back(false);
            
            for (int i = 0; i < M.nb_vertices(); ++i){
            	values[i] = zero_fct(M.get_vertex(i));
            }
            
            // APPLICATION DES CONDITIONS DE NEUMANN
            for (int i = 0; i < M.nb_edges(); ++i){
				if( M.get_edge_attribute(i) == 1 ){
					ElementMapping element( M, true, i);
					std::vector < double > Fe;
					assemble_elementary_neumann_vector(element, shape_1D, quad_1D, neum, Fe);
					local_to_global_vector( M, true, i, Fe, F);
				}
			}
			
			// APPLICATION DES CONDITIONS DE DIRICHLET
			apply_dirichlet_boundary_conditions( M, attr_dirich, values, K, F);
			
			// CALCUL DE LA SOLUTION
			std::vector <double> U;
			solve(K, F, U);
			save_solution(U, "neumann_square_fine.bb");
			
			// AFFICHE LA SOLUTION DANS LE TERMINAL SI -v
			if ( verbose ) {
				std::cout << "Affichage de la solution :" << std::endl;
				for (int i = 0; i < U.size(); ++i){ 
					std::cout << U[i] << std::endl;
				}
            }
            
        }
        
        // -------------------------------------------------------------------- Mug
        
        void mug_pb( const std::string& mesh_filename, bool verbose )
        {
            std::cout << "Solving a neumann problem \n" << std::endl;
            
            // CRÉATION DU MAILLAGE
            Mesh M;
            M.load(mesh_filename);
            M.set_attribute( unit_fct,  1, true); // on applique au bord qui ne sont pas au contact de l'eau bouillante des conditions de neumann
            M.set_attribute( mug_eau_bouill,  0, true); // on applique au bord qui sont au contact de l'eau des conditions de dirichlet

			// INITIALISATION DE LA MATRICE K ET DU VECTEUR F
            std::vector <double> F(M.nb_vertices(), 0);
            SparseMatrix K(M.nb_vertices());
            
            // INITIALISATION DES QUADRATURES 1D ET 2D
			Quadrature quad;
			Quadrature quad_1D;
			
			// CRÉATION DE LA QUADRATURE D'ORDRE 2 ET DES FONCTIONS DE FORME
			quad = quad.get_quadrature(2, false);
			ShapeFunctions shape = ShapeFunctions(2,1);
			
			// CRÉATION DE LA QUADRATURE 1D D'ORDRE 2 ET DES FONCTIONS DE FORME
			quad_1D = quad_1D.get_quadrature(2, true);
			ShapeFunctions shape_1D = ShapeFunctions(1,1);
            
            //sur chaque triangle
            for (int i = 0; i < M.nb_triangles(); ++i){
            	ElementMapping element( M, false, i);
            	DenseMatrix Ke;
            	Ke.set_size(3,3);
            	std::vector < double > Fe;
            	
            	// CRÉATION DE LA MATRICE Ke ET SON IMPLÉMENTATION DANS K
            	assemble_elementary_matrix( element, shape, quad, unit_fct, Ke);
            	local_to_global_matrix( M, i, Ke, K);
            	
            	// CRÉATION DU VECTEUR Fe ET SON IMPLÉMENTATION DANS F (terme source)
            	assemble_elementary_vector( element, shape, quad, zero_fct, Fe);
            	local_to_global_vector( M, false, i, Fe, F);
            }
            
            // CONFIGURATION DE attr_dirich (0 --> dirichlet, 1 --> neumann, 2 --> neumann nul)
            std::vector < bool > attr_dirich;
            std::vector < double > values(M.nb_vertices());
            attr_dirich.push_back(true);
            attr_dirich.push_back(false);
            
            for (int i = 0; i < M.nb_vertices(); ++i){
            	values[i] = 100*mug_eau_bouill(M.get_vertex(i));
            }
            
            // APPLICATION DES CONDITIONS DE NEUMANN
            for (int i = 0; i < M.nb_edges(); ++i){
				if( M.get_edge_attribute(i) == 1 ){
					ElementMapping element( M, true, i);
					std::vector < double > Fe;
					assemble_elementary_neumann_vector(element, shape_1D, quad_1D, mug_neum, Fe);
					local_to_global_vector( M, true, i, Fe, F);
				}
			}
			
			// APPLICATION DES CONDITIONS DE DIRICHLET
			apply_dirichlet_boundary_conditions( M, attr_dirich, values, K, F);
			
			// CALCUL DE LA SOLUTION
			std::vector <double> U;
			solve(K, F, U);
			save_solution(U, "mug_1.bb");
			
			// AFFICHE LA SOLUTION DANS LE TERMINAL SI -v
			if ( verbose ) {
				std::cout << "Affichage de la solution :" << std::endl;
				for (int i = 0; i < U.size(); ++i){ 
					std::cout << U[i] << std::endl;
				}
            }
            
        }

    }

}
