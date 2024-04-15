#pragma once

#include "mesh.h"
#include "fem.h"
#include "solver.h"

#include <assert.h>
#include <iostream>
#include <iomanip>
#include <cmath>
#include <algorithm>
#include <stdlib.h>

namespace FEM2A {
    namespace Tests {

        bool test_load_mesh()
        {
            Mesh mesh;
            mesh.load("data/square.mesh");

            std::cout << "Vertices <x> <y> <att>" << std::endl;
            for( int v = 0; v < mesh.nb_vertices(); v++ ) {
                std::cout << mesh.get_vertex(v).x << " " << mesh.get_vertex(v).y
                    << " " << mesh.get_vertex_attribute(v) << std::endl;
            }

            std::cout << "Edges <v0> <v1> <att>" << std::endl ;
            for( int ed = 0; ed < mesh.nb_edges(); ed++ ) {
                std::cout << mesh.get_edge_vertex_index(ed, 0) << " "
                    << mesh.get_edge_vertex_index(ed, 1) << " "
                    << mesh.get_edge_attribute(ed) << std::endl;
            }

            std::cout << "Triangles <v0> <v1> <v2> <att>" << std::endl ;
            for( int tr = 0; tr < mesh.nb_triangles(); tr++ ) {
                std::cout << mesh.get_triangle_vertex_index(tr, 0) << " "
                    << mesh.get_triangle_vertex_index(tr, 1) << " "
                    << mesh.get_triangle_vertex_index(tr, 2) << " "
                    << mesh.get_triangle_attribute(tr) << std::endl;
            }

            return true;
        }

        bool test_load_save_mesh()
        {
            Mesh mesh;
            mesh.load("data/geothermie_4.mesh");
            mesh.save("data/geothermie_4.mesh");
            return true;
        }
        
        bool test_quadrature()
        {
        	std::cout << "Test Quadrature :" << std::endl;
        	double poids = 0;
        	Quadrature Total;
        	Total = Quadrature::get_quadrature(2);
        	for(int i = 0; i < Total.Quadrature::nb_points(); i++){
        		poids += Total.Quadrature::weight(i);
        	}
        	std::cout<<"Le poids total vaut "<< poids << std::endl;
        	return true;
        }
        
        bool test_element_edge()
        {
        	std::cout << "Test element edge :" << std::endl;
        	Mesh mesh;
           	mesh.load("data/square.mesh");
            ElementMapping element( mesh, true, 4);
        	vertex point;
        	point.x = 0.2;
        	point.y = 0;
        	vertex r;
        	r = element.transform( point );
        	std::cout << "transform : " << r.x << " " << r.y << std::endl;
        	return true;
        }
        
        bool test_element_triangle()
        {
        	std::cout << "Test element triangle :" << std::endl;
        	Mesh mesh;
           	mesh.load("data/square.mesh");
            ElementMapping element(mesh, false, 4);
        	vertex point;
        	point.x = 0.2;
        	point.y = 0.4;
        	vertex r = element.transform( point );
        	std::cout << "transform : " << r.x << " " << r.y << std::endl;
        	return true;
        }
        
        bool test_jacob_edge()
        {
        	std::cout << "Test Jacob edge :" << std::endl;
        	Mesh mesh;
           	mesh.load("data/square.mesh");
            ElementMapping element(mesh, true, 4);
        	vertex point;
        	point.x = 0.2;
        	point.y = 0.4;
        	DenseMatrix J;
        	J = element.jacobian_matrix(point);
        	std::cout << "J vaut :\n" << J.get(0, 0) << "\n" << J.get(0, 1) << std::endl;
        	std::cout << "determinant : " << element.jacobian(point) << std::endl; 
        	return true;
        }
        
        bool test_jacob_triangle()
        {
        	std::cout << "Test Jacob triangle :" << std::endl;
        	Mesh mesh;
           	mesh.load("data/square.mesh");
            	ElementMapping element(mesh, false, 4);
        	vertex point;
        	point.x = 0.2;
        	point.y = 0.4;
        	DenseMatrix J;
        	J = element.jacobian_matrix(point);
        	std::cout << "J vaut :\n" << J.get(0, 0) << "\n" << J.get(0, 1) << "\n" << J.get(1, 0) << "\n" << J.get(1, 1) <<std::endl; 
        	std::cout << "determinant : " << element.jacobian(point) << std::endl;
        	return true;
        }
        
        double unit_fct( vertex v )
        {
            return 1.;
        }
        
        bool test_assemble()
        {
        	std::cout << "Test assemble :" << std::endl;
        	Mesh mesh;
           	mesh.load("data/square.mesh");
            	ElementMapping element(mesh, false, 4);
            	DenseMatrix Ke;
            	Ke.set_size(3, 3);
            	for (int i = 0; i < 3; ++i){
            		for (int j = 0; j < 3; ++j){
            			Ke.set(i, j, 0);
            		}
            	}
            	ShapeFunctions shape(2, 1);
            	Quadrature Q;
            	Q = Quadrature::get_quadrature(2);
            	
            	assemble_elementary_matrix(element, shape, Q, unit_fct, Ke );
            	for (int i = 0; i < 3; ++i){
            		for (int j = 0; j < 3; ++j){
            			std::cout << Ke.get(i, j) << " ";
            		}
            		std::cout << std::endl;
            	}
            	
            	return true;
        }
        	
    }
}
