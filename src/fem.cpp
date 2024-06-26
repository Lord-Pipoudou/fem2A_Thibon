#include "fem.h"
#include "mesh.h"

#include <iomanip>
#include <iostream>
#include <cmath>
#include <limits>
#include <stdlib.h>
#include <assert.h>
#include <vector>


namespace FEM2A {

    void print( const std::vector<double>& x )
    {
        for ( int i = 0; i < x.size(); ++i ) {
            std::cout << x[i] << " ";
        }
        std::cout << std::endl;
    }

    /****************************************************************/
    /* Implementation of Quadrature */
    /****************************************************************/
    int Quadrature::nb_points() const
    {
        return wxy_.size() / 3 ;
    }

    vertex Quadrature::point( int i ) const
    {
        assert( i < nb_points() ) ;
        vertex v ;
        v.x = wxy_[3 * i + 1] ;
        v.y = wxy_[3 * i + 2] ;
        return v ;
    }

    double Quadrature::weight( int i ) const
    {
        assert( i < nb_points() ) ;
        return wxy_[3 * i + 0] ;
    }

    const double triangle_P0[3] = {
        0.5, 0.333333333333333, 0.333333333333333
    };

    const double triangle_P2[9] = {
        0.166666666666667, 0.166666666666667, 0.166666666666667,
        0.166666666666667, 0.166666666666667, 0.666666666666667,
        0.166666666666667, 0.666666666666667, 0.166666666666667
    };

    const double triangle_P4[18] = {
        0.0549758718276609, 0.0915762135097707, 0.0915762135097707,
        0.0549758718276609, 0.0915762135097707, 0.816847572980459,
        0.0549758718276609, 0.816847572980459, 0.0915762135097707,
        0.111690794839006, 0.445948490915965, 0.445948490915965,
        0.111690794839006, 0.445948490915965, 0.10810301816807,
        0.111690794839006, 0.10810301816807, 0.445948490915965
    };

    const double triangle_P6[36] = {
        0.0254224531851034, 0.0630890144915022, 0.0630890144915022,
        0.0254224531851034, 0.0630890144915022, 0.873821971016996,
        0.0254224531851034, 0.873821971016996, 0.0630890144915022,
        0.0583931378631897, 0.24928674517091, 0.24928674517091,
        0.0583931378631897, 0.24928674517091, 0.501426509658179,
        0.0583931378631897, 0.501426509658179, 0.24928674517091,
        0.0414255378091868, 0.0531450498448169, 0.310352451033784,
        0.0414255378091868, 0.310352451033784, 0.0531450498448169,
        0.0414255378091868, 0.0531450498448169, 0.636502499121399,
        0.0414255378091868, 0.636502499121399, 0.0531450498448169,
        0.0414255378091868, 0.310352451033784, 0.636502499121399,
        0.0414255378091868, 0.636502499121399, 0.310352451033784
    };

    const double segment_P0[2] = {
        1., 0.5
    };

    const double segment_P2[4] = {
        0.5, 0.21132486540518708,
        0.5, 0.7886751345948129
    };

    Quadrature Quadrature::get_quadrature( int order, bool border )
    {
        double* pts = NULL;
        int nb_pts = 0;
        Quadrature Q;
        if ( order == 0 && !border ) {
            pts = const_cast<double*>(triangle_P0);
            nb_pts = 1;
        } else if ( order == 2 && !border ) {
            pts = const_cast<double*>(triangle_P2);
            nb_pts = 3;
        } else if ( order == 4 && !border ) {
            pts = const_cast<double*>(triangle_P4);
            nb_pts = 6;
        } else if ( order == 6 && !border ) {
            pts = const_cast<double*>(triangle_P6);
            nb_pts = 12;
        } else if ( order == 0 && border ) {
            pts = const_cast<double*>(segment_P0);
            nb_pts = 1;
        } else if ( order == 2 && border ) {
            pts = const_cast<double*>(segment_P2);
            nb_pts = 2;
        } else {
            std::cout << "Quadrature not implemented for order " << order << std::endl;
            assert( false );
        }
        Q.wxy_.resize(nb_pts * 3);
        for ( int i = 0; i < nb_pts; ++i ) {
            if ( !border ) {
                Q.wxy_[3*i+0] = pts[3*i+0];
                Q.wxy_[3*i+1] = pts[3*i+1];
                Q.wxy_[3*i+2] = pts[3*i+2];
            } else {
                Q.wxy_[3*i+0] = pts[2*i+0];
                Q.wxy_[3*i+1] = pts[2*i+1];
                Q.wxy_[3*i+2] = 0.;
            }
        }
        return Q;
    }

    /****************************************************************/
    /* Implementation of ElementMapping */
    /****************************************************************/
    ElementMapping::ElementMapping( const Mesh& M, bool border, int i )
        : border_( border )
    {
        //std::cout << "[ElementMapping] constructor for element " << i << " ";
        std::vector< vertex > vertices;
        if ( border ) {
        	for( int v=0; v<2; v++){
        		vertices.push_back( M.get_edge_vertex(i, v));
        	}
        	vertices_ = vertices; 
        	
        }
        else {
        	for( int v=0; v < 3; v++){
        		vertices.push_back( M.get_triangle_vertex(i, v));
        	}
        	vertices_ = vertices;
        }
        		
        //std::cout << '\n';
        // TODO
    }

    vertex ElementMapping::transform( vertex x_r ) const
    {
        //std::cout << "[ElementMapping] transform reference to world space" << '\n';
        // TODO
        vertex r ;
        double xi = x_r.x;
        double eta = x_r.y;
        if (border_){
        	r.x = (vertices_[0].x)*(1-xi) + (vertices_[1].x)*xi;
        	r.y = (vertices_[0].y)*(1-xi) + (vertices_[1].y)*xi;
        	
        }
        
        else {
        	r.x = (vertices_[0].x)*(1.0-xi-eta) + (vertices_[1].x)*xi + (vertices_[2].x)*eta;
        	r.y = (vertices_[0].y)*(1.0-xi-eta) + (vertices_[1].y)*xi + (vertices_[2].y)*eta;
        }
        
        return r ;
    }

    DenseMatrix ElementMapping::jacobian_matrix( vertex x_r ) const
    {
        //std::cout << "[ElementMapping] compute jacobian matrix" << '\n';
        // TODO
        DenseMatrix J ;
        if (border_){
        	J.set_size(1,2);
        	J.set(0, 0, vertices_[1].x-vertices_[0].x);
        	J.set(0, 1, vertices_[1].y-vertices_[0].y);
        }
        else {
        	J.set_size(2,2);
        	J.set(0, 0, vertices_[1].x-vertices_[0].x);
        	J.set(1, 0, vertices_[1].y-vertices_[0].y);
        	J.set(0, 1, vertices_[2].x-vertices_[0].x);
        	J.set(1, 1, vertices_[2].y-vertices_[0].y);
        }
        return J ;
    }

    double ElementMapping::jacobian( vertex x_r ) const
    {
        //std::cout << "[ElementMapping] compute jacobian determinant" << '\n';
        // TODO
        DenseMatrix J = jacobian_matrix( x_r );
        double determinant;
        if (border_){
        	determinant = pow(0.5, pow(2, J.get(0, 0)) + pow(2, J.get(0, 1)));
        }
        else {
        	determinant = J.det_2x2();
        }
        return determinant ;
    }

    /****************************************************************/
    /* Implementation of ShapeFunctions */
    /****************************************************************/
    ShapeFunctions::ShapeFunctions( int dim, int order )
        : dim_( dim ), order_( order )
    {
        if (dim != 1 && dim != 2 || order != 1){std::cout << "alerte valeur de dim ou order impossible" << "\n";}
    }

    int ShapeFunctions::nb_functions() const
    {
        return dim_+1 ;
    }

    double ShapeFunctions::evaluate( int i, vertex x_r ) const
    {
        if (dim_ == 2){
        	switch(i){
        		case 0: return 1-x_r.x-x_r.y;
        		case 1: return x_r.x;
        		case 2: return x_r.y;
        	}
        }
        if (dim_ == 1){
        	switch(i){
        		case 0: return 1-x_r.x;
        		case 1: return x_r.x;
        	}
        }
        return 0. ; // should not be reached
    }

    vec2 ShapeFunctions::evaluate_grad( int i, vertex x_r ) const
    {
        vec2 g ;
        if (dim_ == 2){
        	switch(i){
        		case 0:
        		g.x = -1;
        		g.y = -1;
        		break;
        		case 1:
        		g.x = 1;
        		g.y = 0;
        		break;
        		case 2:
        		g.x = 0;
        		g.y = 1;
        		break;
        	}
        }
        if (dim_ == 1){
        	switch(i){
        		case 0:
        		g.x = -1;
        		g.y = 0;
        		break;
        		case 1:
        		g.x = 1;
        		g.y = 0;
        		break;
        	}
        }
        return g ;
    }

    /****************************************************************/
    /* Implementation of Finite Element functions */
    /****************************************************************/
    void assemble_elementary_matrix(
        const ElementMapping& elt_mapping,
        const ShapeFunctions& reference_functions,
        const Quadrature& quadrature,
        double (*coefficient)(vertex),
        DenseMatrix& Ke )
    {
        //std::cout << "compute elementary matrix" << '\n';
        // TODO
        int max = reference_functions.nb_functions();
        int nbr_quad = quadrature.nb_points(); 
        
        //PARCOURS DES FONCTIONS DE REFERENCE + GRAD SUIVANT I OU J
        for (int i = 0; i < max; ++i){
        	for (int j = 0; j < max; ++j){
        	
        	// PARCOURS DES POINTS DE QUADRATURE
        		for (int q = 0; q < nbr_quad; ++q){
        			vertex quad = quadrature.point(q);
        			double w = quadrature.weight(q); 
        			vec2 gradi = reference_functions.evaluate_grad(i, quad);
        			vec2 gradj = reference_functions.evaluate_grad(j, quad);
        			DenseMatrix J = elt_mapping.jacobian_matrix(quad);
        			DenseMatrix JinvT = J.invert_2x2().transpose();
        			double det = elt_mapping.jacobian(quad);
        			vertex Me = elt_mapping.transform(quad);
        			vec2 pdti = JinvT.mult_2x2_2(gradi);
        			vec2 pdtj = JinvT.mult_2x2_2(gradj);
        			
        			// CALCUL DE LA i-j-ème COMPOSANTE DE Ke pour le q-ième point de quadrature 
        			double Kij = w*coefficient(Me)*(pdti.x*pdtj.x + pdti.y*pdtj.y)*det;
        			Ke.add(i, j, Kij);
        		} 
        	}
        }
        
        
    }

    

    void local_to_global_matrix(
        const Mesh& M,
        int t,
        const DenseMatrix& Ke,
        SparseMatrix& K )
    {    
        for (int i = 0; i < Ke.height(); i++){
        	for (int j = 0; j < Ke.width(); j++){
          		K.add(M.get_triangle_vertex_index(t,i), M.get_triangle_vertex_index(t,j), Ke.get(i,j));
        	}
        }
    }

    void assemble_elementary_vector(
        const ElementMapping& elt_mapping,
        const ShapeFunctions& reference_functions,
        const Quadrature& quadrature,
        double (*source)(vertex),
        std::vector< double >& Fe )
    {
        //std::cout << "compute elementary vector (source term)" << '\n';
        // TODO
        
        int nbr_quad = quadrature.nb_points(); 
        int max = reference_functions.nb_functions();
        
        //PARCOURS DES FONCTIONS DE REFERENCE
        for (int i = 0; i < max; ++i){
        	double Fe_i = 0;
        	
        	// PARCOURS DES POINTS DE QUADRATURE
        	for (int q = 0; q < nbr_quad; ++q){
        		vertex quad = quadrature.point(q);
        		double w = quadrature.weight(q); 
        		double shape_i = reference_functions.evaluate(i, quad);
        		double det = elt_mapping.jacobian(quad);
        		vertex Me = elt_mapping.transform(quad);
        		
        		// CALCUL DE LA i-ème COMPOSANTE DE Fe
        		Fe_i += w*shape_i*source(Me)*det;
        	}
    		Fe.push_back(Fe_i);
        }
    }

    void assemble_elementary_neumann_vector(
        const ElementMapping& elt_mapping_1D,
        const ShapeFunctions& reference_functions_1D,
        const Quadrature& quadrature_1D,
        double (*neumann)(vertex),
        std::vector< double >& Fe )
    {
        //std::cout << "compute elementary vector (neumann condition)" << '\n';
        // TODO
        
        int nbr_quad = quadrature_1D.nb_points(); 
        int max = reference_functions_1D.nb_functions();
        
        //PARCOURS DES FONCTIONS DE REFERENCE
        for (int i = 0; i < max; ++i){
        	double Fe_i = 0;
        	
        	// PARCOURS DES POINTS DE QUADRATURE
        	for (int q = 0; q < nbr_quad; ++q){
        		vertex quad = quadrature_1D.point(q);
        		double w = quadrature_1D.weight(q); 
        		double shape_i = reference_functions_1D.evaluate(i, quad);
        		double det = elt_mapping_1D.jacobian(quad);
        		vertex Me = elt_mapping_1D.transform(quad);
        		
        		// CALCUL DE LA i-ème COMPOSANTE DE Fe
        		Fe_i += w*shape_i*neumann(Me)*det;
        	}
    		Fe.push_back(Fe_i);
        }
        
        
    }

    void local_to_global_vector(
        const Mesh& M,
        bool border,
        int i,
        std::vector< double >& Fe,
        std::vector< double >& F )
    {
        //std::cout << "Fe -> F" << '\n';
        // TODO
        
        int indice;
        if (border){
        	for (int j = 0; j < Fe.size(); ++j){
        		indice = M.get_edge_vertex_index( i, j );
        		F[indice] += Fe[j];
        	}
        }
        else {
        	for (int j = 0; j < Fe.size(); ++j){
        		indice = M.get_triangle_vertex_index(i,j);
        		F[indice] += Fe[j];
        	}
        }
        
    }

    void apply_dirichlet_boundary_conditions(
        const Mesh& M,
        const std::vector< bool >& attribute_is_dirichlet, /* size: nb of attributes */
        const std::vector< double >& values, /* size: nb of DOFs */
        SparseMatrix& K,
        std::vector< double >& F )
    {
        // std::cout << "apply dirichlet boundary conditions" << '\n';
        
        std::vector < bool > processed_verticies(values.size(), false);
        double P = 10000;
        for( int edge = 0; edge < M.nb_edges(); edge++){
        	int edge_attribute = M.get_edge_attribute(edge);
        	if ( attribute_is_dirichlet[edge_attribute] ) {
        		for( int v = 0; v < 2; v++ ){
        			int vertex_index = M.get_edge_vertex_index(edge, v);
        			if( !processed_verticies[vertex_index] ){
        				processed_verticies[vertex_index] = true;
        				K.add(vertex_index, vertex_index, P);
        				F[vertex_index] += P*values[vertex_index];
        			}
        		}
        	}
        }
    }

/*
    void solve_poisson_problem(
            const Mesh& M,
            double (*diffusion_coef)(vertex),
            double (*source_term)(vertex),
            double (*dirichlet_fct)(vertex),
            double (*neumann_fct)(vertex),
            std::vector<double>& solution,
            int verbose )
    {
        std::cout << "solve poisson problem" << '\n';

	    
	    //M.set_attribute( Simu::unit_fct,  1, true);
	    //M.set_attribute( Simu::mug_eau_bouill,  0, true);

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
	    	assemble_elementary_matrix( element, shape, quad, diffusion_coef, Ke);
	    	local_to_global_matrix( M, i, Ke, K);
	    	
	    	//assemble_elementary_vector( element, shape, quad, Simu::source_fct, Fe);
	    	local_to_global_vector( M, false, i, Fe, F);
	    }
	    
	    std::vector < bool > attr_dirich;
	    std::vector < double > values(M.nb_vertices());
	    attr_dirich.push_back(true);
	    attr_dirich.push_back(false);
	    for (int i = 0; i < M.nb_vertices(); ++i){
	    	values[i] = dirichlet_fct(M.get_vertex(i));
	    }
	    for (int i = 0; i < M.nb_edges(); ++i){
			if( M.get_edge_attribute(i) == 1 ){
				ElementMapping element( M, true, i);
				std::vector < double > Fe;
				assemble_elementary_neumann_vector(element, shape_1D, quad_1D, neumann_fct, Fe);
				local_to_global_vector( M, true, i, Fe, F);
			}
		}	
			
		apply_dirichlet_boundary_conditions( M, attr_dirich, values, K, F);
		std::vector <double> U;
		solve(K, F, U);
		save_solution(U, "mug_0_5.bb");
		if ( verbose ) {
			for (int i = 0; i < U.size(); ++i){ 
				std::cout << U[i] << std::endl;
			}
		}
		
	}
*/

}
