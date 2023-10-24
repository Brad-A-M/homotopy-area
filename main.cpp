// Copyright (c) 2021,2022,2023 Montana State (USA).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s)     : Bradley McCoy <bradleymccoy@montana.edu>
//



// Adapting the dual of an arrangement to a BGL graph.
#include <CGAL/config.h>
#include <boost/graph/breadth_first_search.hpp>
#include <boost/graph/visitors.hpp>
#include <CGAL/Arr_extended_dcel.h>
#include <CGAL/graph_traits_dual_arrangement_2.h>
#include <CGAL/Arr_face_index_map.h>
#include "Extended_face_property_map.h"
#include "arr_exact_construction_segments.h"
#include "arr_print.h"
#include <CGAL/HalfedgeDS_halfedge_base.h>
#include <CGAL/HalfedgeDS_decorator.h>
#include <CGAL/HalfedgeDS_face_base.h>
#include <CGAL/Arrangement_on_surface_2.h>
#include <CGAL/Arr_dcel_base.h>
#include <chrono>

#include <iostream>
#include <vector>
#include <CGAL/Exact_rational.h>

//class for edge data (cables, orientation)
//(vector of ints - cables, int orientation of curve)

namespace CGAL { namespace internal {
class Half_edge_DS {
public:
    std::vector<int> cables;
    int ort;
    
    Half_edge_DS(){
        ort = 0;
    }
    Half_edge_DS(std::vector<int> c, int r){
        cables = c;
        ort = r;
    }
    
    std::vector<int> get_cables() {return cables;}
    void update_cables(std::vector<int> cabs){cables = cabs;}
    void add_cable(int cb){cables.push_back(cb);}
    void insert_cable(int cb){cables.insert(cables.begin(),-cb);}
    int get_ort() const {return ort;}
    void  set_ort(int tation){ort = tation;}
    
    void  print() {
        std::cout << "Cables vector: ";
        for (auto i = cables.begin(); i != cables.end(); ++i)
            std::cout << *i << " ";
    }
    
};
}}

typedef CGAL::Arr_extended_dcel<Traits, int, CGAL::internal::Half_edge_DS, int> Dcel;
typedef CGAL::Arrangement_2<Traits, Dcel>                  Ex_arrangement;
typedef CGAL::Dual<Ex_arrangement>                         Dual_arrangement;
typedef CGAL::Arr_face_index_map<Ex_arrangement>           Face_index_map;
typedef Extended_face_property_map<Ex_arrangement,unsigned int>
                                                           Face_property_map;

//For the area of a polygon
//Potential problem: the faces are weakly simple polygons
//but not simple. Shoelace formula still works but I'm scared.
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Polygon_2.h>
#include <CGAL/Point_2.h>
#include <list>
typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef CGAL::Point_2<K> Points;
typedef CGAL::Polygon_2<K> Polygon_2;


//Class to hold face data a face is (id,area,edge,flag)
//f=(int - depth in a BFS search, double a - area, Ex_arrangement::Ccb_halfedge_circulator - edge to face discoverd earlier, bool flag for exploration)
namespace CGAL { namespace internal {
class Face_DS {
public:
    int id;
    double area;
    Ex_arrangement::Ccb_halfedge_circulator edge;
    bool explored;
    
    Face_DS(){
        id = 0;
        area = 0;
        explored = 0;
    }
    Face_DS(int f, double a, Ex_arrangement::Ccb_halfedge_circulator g, bool e){
        id = f;
        area = a;
        edge = g;
        explored = e;
    }
    
    int get_id() const {return id;}
    void set_id(int idd){id=idd;}
    double get_area() const {return area;}
    void  set_area(double ara){area = ara;}
    Ex_arrangement::Ccb_halfedge_circulator get_edge(){return edge;}
    void set_edge(Ex_arrangement::Ccb_halfedge_circulator edg){edge=edg;}
    bool is_explored() {return explored;}
    void set_explored(bool ex){explored=ex;}
    
    void  print() {
        
        std::cout << "face: "<< id << " has area: " << area<<
        " edge: ("<< edge->source()->point() <<") to ("<< edge->target()->point() <<") explored: "<<explored<<
        std::endl;
    }
    
};
}}

std::vector<int> insert_negative(std::vector<int> vec1, std::vector<int> vec2){
    for (auto i = vec2.begin(); i != vec2.end(); ++i)
        vec1.insert(vec1.begin(),-*i);
    return vec1;
}


int main() {
  
    //start timer
    auto start = std::chrono::high_resolution_clock::now();

  // Construct an curve as anarrangement
    /*
     //Example 1 the eight area should be 8
    Point p1(1, 0), p2(2, 2), p3(0,3), p4(-2, 2), p5(-1, 0), p6(1, 1), p7(-1, 2), p8(1, 2),p9(-1, 1),p10(1, 0);
  Traits traits;
  Ex_arrangement  arr(&traits);
  insert(arr, Segment(p1, p2));
  insert(arr, Segment(p2, p3));  insert(arr, Segment(p3, p4));
  insert(arr, Segment(p4, p5));  insert(arr, Segment(p5, p6));
  insert(arr, Segment(p6, p7));  insert(arr, Segment(p7, p8));
  insert(arr, Segment(p8, p9));  insert(arr, Segment(p9, p10));
     */
    //Example 2 double eight (area should be 18)
    Point p1(1, 0), p2(2, 1), p3(1,2), p4(2, 2), p5(1, 1), p6(2, 0), p7(3,0), p8(3, 3),p9(-3, 3),p10(-3, 0),p11(-2,0),p12(-1,1),p13(-2,2),p14(-1,2),p15(-2,1),p16(-1,0);
  Traits traits;
  Ex_arrangement  arr(&traits);
  insert(arr, Segment(p1, p2));
  insert(arr, Segment(p2, p3));  insert(arr, Segment(p3, p4));
  insert(arr, Segment(p4, p5));  insert(arr, Segment(p5, p6));
  insert(arr, Segment(p6, p7));  insert(arr, Segment(p7, p8));
  insert(arr, Segment(p8, p9));  insert(arr, Segment(p9, p10));
  insert(arr, Segment(p10, p11));  insert(arr, Segment(p11, p12));
  insert(arr, Segment(p12, p13));  insert(arr, Segment(p13, p14));
  insert(arr, Segment(p14, p15));  insert(arr, Segment(p15, p16));
    insert(arr, Segment(p16, p1));
    /*
    //Example 3 mouse, note unequal ear size, area should be 25.9
    Point p1(1, 0), p2(3, 0), p3(3,4), p4(-3, 4), p5(-3, 0), p6(-1, 0), p7(1,1), p8(1, 3),p9(2.2, 2),p10(-2, 2),p11(-1,3),p12(-1,1);
  Traits traits;
  Ex_arrangement  arr(&traits);
  insert(arr, Segment(p1, p2));
  insert(arr, Segment(p2, p3));  insert(arr, Segment(p3, p4));
  insert(arr, Segment(p4, p5));  insert(arr, Segment(p5, p6));
  insert(arr, Segment(p6, p7));  insert(arr, Segment(p7, p8));
  insert(arr, Segment(p8, p9));  insert(arr, Segment(p9, p10));
  insert(arr, Segment(p10, p11));  insert(arr, Segment(p11, p12));
  insert(arr, Segment(p12, p1));
     */
    
    
    
  //%%%%%%%%%%%%%%%%%%%%%%%%%%%
  //Begin: Build array of phaces
  //BFS to mape faces to discover time.
  // Create a mapping of the arrangement faces to indices.
  Face_index_map index_map(arr);
  // Perform breadth-first search from the unbounded face, using the event
  // visitor to associate each arrangement face with its discover time.
  int time = -1;
  boost::breadth_first_search(Dual_arrangement(arr), arr.unbounded_face(),
                              boost::vertex_index_map(index_map).visitor
                              (boost::make_bfs_visitor
                               (stamp_times(Face_property_map(), time,
                                            boost::on_discover_vertex()))));
    
    //%%%%%%%%%%%%%%%%%%%%%%%%%%%
    //End BFS
    //Begin: Build array of phaces and exterior edge
    //%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    CGAL::internal::Face_DS phaces[arr.number_of_faces()];
    CGAL::internal::Face_DS loader;
    Ex_arrangement::Ccb_halfedge_circulator loadedge;

    int i = 0;
    //std::cout << "We are here"<< std::endl;
    for (auto fit = arr.faces_begin(); fit != arr.faces_end(); ++fit) {
       // std::cout << "Face no. " << fit->data() << ": ";
        if (fit->is_unbounded()){
            // std::cout << "Face no. " << fit->data()<< " Unbounded.\n";
        }
            
        else {
            auto currr = fit->outer_ccb();
            int neighbor = 10e5;
            //Find edge with minimal neighbor id
            //Can improve by just find one less than current id.
            if(currr->twin()->face()->data()<neighbor)
              {
                  loadedge = currr;
                  neighbor = currr->twin()->face()->data();
              }
            
            while (++currr != fit->outer_ccb())
                if(currr->twin()->face()->data()<neighbor)
                  {
                      loadedge = currr;
                      neighbor = currr->twin()->face()->data();
                  }
            phaces[fit->data()-1] = CGAL::internal::Face_DS(fit->data(), 5.0, loadedge,0);
           // std::cout << loadedge->source()->point();
           // std::cout << std::endl;
        }
        
    }
    
   // for (int i=0;i<arr.number_of_faces()-1;i++) {
    //        phaces[i].print();
    //    }
    //std::cout << "I am here"<< std::endl;
    //%%%%%%%%%%%%%%%%%%%%%%%%%%%
    //Begin finding area of faces
    double x;
    double y;
    
    for (auto fit = arr.faces_begin(); fit != arr.faces_end(); ++fit) {
        Polygon_2 poly;
        
        if (fit->is_unbounded()){
            //std::cout << "Face no. " << fit->data()<< " Unbounded."<< std::endl;
            }
        else{
           // std::cout << "Face no. " << fit->data()<< std::endl;
            auto bdy = fit->outer_ccb();
            x =CGAL::to_double(bdy->source()->point().hx());
            y =CGAL::to_double(bdy->source()->point().hy());
                do{
                    x =CGAL::to_double(bdy->source()->point().hx());
                    y =CGAL::to_double(bdy->source()->point().hy());
                    //std::cout <<"("<< x << ","<< y<< ")"<< std::endl;
                    poly.push_back(Points(x,y));
                    ++bdy;
                }while (bdy != fit->outer_ccb());
            
            double Area = poly.area();
            //std::cout <<"(Area)"<< Area<<" face "<< fit->data()<<" is simple "<<poly.is_simple()<<std::endl;
            phaces[fit->data()-1].set_area(Area);
        }
    }
        
    for (int i=0;i<arr.number_of_faces()-1;i++) {
            phaces[i].print();
        }
    
    //%%%%%%%%%%%%%%%%%%%%%%%%%%%
    //End: Build array of phaces
    //Begin: Edges
    //%%%%%%%%%%%%%%%%%%%%%%%%%%%
    std::cout << "%%%%%%%%%%%%%%%%%%%%%%%%%%%"<< std::endl;
    
    
    //pick an edge to begin
    Ex_arrangement::Edge_iterator              eit;
    Ex_arrangement::Halfedge_handle crnt, upedge;
    
    eit=arr.edges_begin();
    
    //make sure tip has degree 2
    while(eit->target()->degree()!=2)
        {++eit;}
  //  std::cout << "eit: [" << eit->curve() <<"] "<<std::endl;
   // eit->data().print();
    eit->data().set_ort(1);
    crnt=eit->next();
   // std::cout << "edge into loop: [" << crnt->curve() <<"] "<<std::endl;
  //  crnt->data().print();
    crnt->data().set_ort(1);

    
    Ex_arrangement::Vertex_handle v,u;
    Ex_arrangement::Halfedge_around_vertex_circulator first1, crnt1;
    bool    eq;
    bool    co;
    
    while(crnt!=eit)
        {
            crnt->data().set_ort(1);
           // std::cout << "next edge: [" << crnt->curve() <<"] "<<std::endl;
            
            //crnt->data().print();

            v = crnt->target();
            first1 = crnt1 = v->incident_halfedges();
            //std::cout << "crnt1: [" << crnt1->curve() <<"] "<<std::endl;
            
            if(crnt->target()->degree()==2)
            {
                crnt=crnt->next();
            }
             else{
                do {
                   //iterate edges around target(v), find
                   //edge that is colinear but not the same
                   // std::cout << "inside crnt1: [" << crnt1->curve() <<"] "<<std::endl;
                    // Note that the current halfedge is directed from u to v:
                    u = crnt1->source();
                    
                    eq = (u->point()==crnt->source()->point());
                    co = CGAL::collinear(crnt->source()->point(), v->point(), u->point());
                    
                    if(!eq && co){
                        
                        upedge=crnt1->twin();
                    }
                    
                } while (++crnt1 != first1);
                 crnt=upedge;
            }
        }
    
    std::cout << "%%%%%%%%%%%%%%%%%%%%%%%%%%%"<< std::endl;
    
    //%%%%%%%%%%%%%%%%%%%%%%%%%%%
    //End: Traverse cuver (assign orientation)
    //Begin: Laying cables.
    //%%%%%%%%%%%%%%%%%%%%%%%%%%%
    

    std::vector<int> newcables, edgecables;
    int j;
    
    //phaces[0].print();
    
    for (int i=arr.number_of_faces()-2;i>-1;i--) {
        std::cout<< i <<std::endl;
        j=i;
        newcables.clear();
        if(phaces[i].is_explored()==0){
            
            phaces[i].set_explored(1);
            if(phaces[i].get_edge()->data().get_ort()==0){
                phaces[i].get_edge()->twin()->data().insert_cable(i+1);
                
            }
            else{
                phaces[i].get_edge()->data().add_cable(i+1);
                
            }
            newcables.push_back(i+1);
            j=phaces[i].get_edge()->twin()->face()->data()-1;
            std::cout<< j <<std::endl;
            
            while(j>=0){
                //std::cout<<"inside j "<< j <<std::endl;
                
                if(phaces[j].is_explored()==0){
                    newcables.push_back(j+1);
                    phaces[j].set_explored(1);
                }
                
                
                if(phaces[j].get_edge()->data().get_ort()==0){
                    //edgecables=phaces[j].get_edge()->twin()->data().get_cables();
                    phaces[j].get_edge()->twin()->data().update_cables(insert_negative(phaces[j].get_edge()->twin()->data().get_cables(),newcables));
                    //phaces[j].get_edge()->twin()->data().insert_cable(i+1);
                    
                }
                else{
                    edgecables=phaces[j].get_edge()->data().get_cables();
                    edgecables.insert(edgecables.end(), newcables.begin(), newcables.end());
                }
                j=phaces[j].get_edge()->twin()->face()->data()-1;

               // std::cout<<"j bottom"<< j <<std::endl;
            }
            
            
        }
            phaces[i].print();
    }
    
    std::cout << "%%%%%%%%%%%%%%%%%%%%%%%%%%%"<< std::endl;
    
    
    std::vector<int> word_vector, sub_word_vector;
    eit=arr.edges_begin();
    
    //make sure tip has degree 2
    while(eit->target()->degree()!=2)
        {++eit;}
    //std::cout << "eit: [" << eit->curve() <<"] "<<std::endl;
   // eit->data().print();
    eit->data().set_ort(1);
    crnt=eit->next();
   // std::cout << "edge into loop: [" << crnt->curve() <<"] "<<std::endl;
   // crnt->data().print();
    //crnt->data().set_ort(1);

    
    //Ex_arrangement::Vertex_handle v,u;
   // Ex_arrangement::Halfedge_around_vertex_circulator first1, crnt1;
   // bool    eq;
    //bool    co;
    
    while(crnt!=eit)
        {
            
            crnt->data().set_ort(1);
            //std::cout << "next edge: [" << crnt->curve() <<"] "<<std::endl;
            
            //crnt->data().print();
            sub_word_vector.clear();
            sub_word_vector= crnt->data().get_cables();
            for (auto i = sub_word_vector.begin(); i != sub_word_vector.end(); ++i){
               word_vector.push_back(*i);
            }


            v = crnt->target();
            first1 = crnt1 = v->incident_halfedges();
            //std::cout << "crnt1: [" << crnt1->curve() <<"] "<<std::endl;
            
            if(crnt->target()->degree()==2)
            {
                crnt=crnt->next();
            }
             else{
                do {
                   //iterate edges around target(v), find
                   //edge that is colinear but not the same
                   // std::cout << "inside crnt1: [" << crnt1->curve() <<"] "<<std::endl;
                    // Note that the current halfedge is directed from u to v:
                    u = crnt1->source();
                    
                    eq = (u->point()==crnt->source()->point());
                    co = CGAL::collinear(crnt->source()->point(), v->point(), u->point());
                    
                    if(!eq && co){
                        
                        upedge=crnt1->twin();
                    }
                    
                } while (++crnt1 != first1);
                 crnt=upedge;
            }
        }
    
    for (auto i = word_vector.begin(); i != word_vector.end(); ++i)
        std::cout << *i << " ";
    
    
    std::cout << "%%%%%%%%%%%%%%%%%%%%%%%% "<<std::endl;
    /*
    std::vector<int> vec3, vec4;
    
    vec3.push_back(4);
    vec3.push_back(2);
    vec4.push_back(7);
    vec4.push_back(9);
    
    for (auto i = vec3.begin(); i != vec3.end(); ++i)
        std::cout << *i << " ";
    std::cout << "now add "<<std::endl;
    vec3=insert_negative(vec3,vec4);
    
    for (auto i = vec3.begin(); i != vec3.end(); ++i)
        std::cout << *i << " ";
    */
    
    /////////////////////////////////
    //begin compute the norm
    /////////////////////////////////
    ///
    int wl=word_vector.size();
    double word[wl][2];
    
    
    for(int i=0;i<wl;i++){
        word[i][0]=word_vector.at(i);
        word[i][1]=phaces[abs(word_vector.at(i))-1].get_area();
        std::cout<< "Face letter " << word[i][0]<< " and area: " <<word[i][1]<<std::endl;
        
    }
    double dynamic[wl][wl];
    
    //load up the diagonal
    for (int i=0; i<wl; i++){
        dynamic[i][i]=word[i][1];
    }
    
    //fill up dynamic matrix
    double l,r,t,w;
    for(int j=1; j<wl; j++){
        for(int i=0; i<(wl-j); i++){
            
            //double w;//initial weight ||w'||+wt(l_n)
            w=dynamic[i][i+j-1]+dynamic[i+j][i+j];
            
            for(int k=i;k<(i+j-1);k++){
                if(word[k][0]==-word[i+j][0]){
                    //double l;
                    //double r;
                    if(k==i){l=0;}
                    else{l=dynamic[i][k-1];}
                    if(k==(i+j-1)){r=0;}
                    else{r=dynamic[k+1][i+j-1];}
                    
                   // double t;
                    t=l+r;
                    
                    if(t<w){w=t;}
                }
                
            }
            dynamic[i][i+j]=w;
            
        }
    }//end fill up loops
    
    std::cout << "The homotopy area is: " << dynamic[0][wl-1] << std::endl;
    
    //stop timer
    auto stop = std::chrono::high_resolution_clock::now();
    
    auto duration = std::chrono::duration_cast<std::chrono::microseconds>(stop - start);
    std::cout << "This took: " << duration.count() << std::endl;
  return 0;
}
