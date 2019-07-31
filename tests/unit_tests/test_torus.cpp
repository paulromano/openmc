#include <iostream>
#include "openmc/surface.h"
#include "pugixml.hpp"

void ray_through_middle() {
    pugi::xml_document doc;

    pugi::xml_node xml_node = doc.append_child("surface");

    xml_node.set_name("surface");
    xml_node.append_attribute("id") = "1";
    xml_node.append_attribute("type") = "z-torus";
    xml_node.append_attribute("coeffs") = "0 0 0 20 5 5";

    openmc::SurfaceZTorus *surf = new openmc::SurfaceZTorus(xml_node);
    
    openmc::Position r = {-200,0,0};
    openmc::Direction w = {1,0,0};
    double d = surf->distance(r,w,false);
    std::cout << d << std::endl;
    r.x += w.x*d;
    d = surf->distance(r,w,false);
    std::cout << d << std::endl;
    r.x += w.x*d;
    d = surf->distance(r,w,false);
    std::cout << d << std::endl;
    r.x += w.x*d;
    d = surf->distance(r,w,false);
    std::cout << d << std::endl;
    r.x += w.x*d;
    d = surf->distance(r,w,false);
    std::cout << d << std::endl;
    return;
}

void ray_across_the_top() {
    pugi::xml_document doc;

    pugi::xml_node xml_node = doc.append_child("surface");

    xml_node.set_name("surface");
    xml_node.append_attribute("id") = "1";
    xml_node.append_attribute("type") = "z-torus";
    xml_node.append_attribute("coeffs") = "0 0 0 20 5 5";

    openmc::SurfaceZTorus *surf = new openmc::SurfaceZTorus(xml_node);
    
    openmc::Position r = {-200,0,5};
    openmc::Direction w = {1,0,0};

    double d = surf->distance(r,w,false);
    std::cout << d << std::endl;
    r.x += w.x*d;
    d = surf->distance(r,w,false);
    std::cout << d << std::endl;
    r.x += w.x*d;
    d = surf->distance(r,w,false);
    std::cout << d << std::endl;

}

void ray_outboard_middle() {
    pugi::xml_document doc;

    pugi::xml_node xml_node = doc.append_child("surface");

    xml_node.set_name("surface");
    xml_node.append_attribute("id") = "1";
    xml_node.append_attribute("type") = "z-torus";
    xml_node.append_attribute("coeffs") = "0 0 0 20 5 5";

    openmc::SurfaceZTorus *surf = new openmc::SurfaceZTorus(xml_node);
    
    openmc::Position r = { 20,0,0};
    openmc::Direction w = {1,0,0};

    double d = surf->distance(r,w,false);
    std::cout << d << std::endl;
    r.x += w.x*d;
    d = surf->distance(r,w,false);
    std::cout << d << std::endl;
}

void ray_normals() {
  pugi::xml_document doc;

  pugi::xml_node xml_node = doc.append_child("surface");

  xml_node.set_name("surface");
  xml_node.append_attribute("id") = "1";
  xml_node.append_attribute("type") = "z-torus";
  xml_node.append_attribute("coeffs") = "0 0 0 20 5 5";
  openmc::SurfaceZTorus *surf = new openmc::SurfaceZTorus(xml_node);
  
  openmc::Position inside = {-20,0,0};
  openmc::Direction normal = surf->normal(inside);
  std::cout << normal.x << " " << normal.y << " " << normal.z << std::endl;

  inside = {-25,0,0};
  normal = surf->normal(inside);
  std::cout << normal.x << " " << normal.y << " " << normal.z << std::endl;

  inside = {-19,0,0};
  normal = surf->normal(inside);
  std::cout << normal.x << " " << normal.y << " " << normal.z << std::endl;


  inside = {-15,0,0};
  normal = surf->normal(inside);
  std::cout << normal.x << " " << normal.y << " " << normal.z << std::endl;

  inside = {-10,0,0};
  normal = surf->normal(inside);
  std::cout << normal.x << " " << normal.y << " " << normal.z << std::endl;


  inside.x = -80;
  normal = surf->normal(inside);
  std::cout << "out " << normal.x << " " << normal.y << " " << normal.z << std::endl;

  inside.x = 20;
  normal = surf->normal(inside);
  std::cout << normal.x << " " << normal.y << " " << normal.z << std::endl;

  inside.x = 80;
  normal = surf->normal(inside);
  std::cout << normal.x << " " << normal.y << " " << normal.z << std::endl;

  inside.x = -20;
  inside.y = 2;
  inside.z = 1;
  normal = surf->normal(inside);
  std::cout << normal.x << " " << normal.y << " " << normal.z << std::endl;

}

int main(int argc, char* argv[]) {

    std::cout << "-------" << std::endl;
    ray_through_middle();
    std::cout << "-------" << std::endl;
    ray_across_the_top();
    std::cout << "-------" << std::endl;
    ray_outboard_middle();
    std::cout << "-------" << std::endl;
    ray_normals();
    return 0;
}