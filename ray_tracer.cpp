/**
 * ray_tracer.cpp
 * CS230 Assignment 2, Winter 2012
 * -------------------------------
 * Implement ray tracer here.
 */

#define SET_RED(P, C)   (P = (((P) & 0x00ffffff) | ((C) << 24)))
#define SET_GREEN(P, C)  (P = (((P) & 0xff00ffff) | ((C) << 16)))
#define SET_BLUE(P, C) (P = (((P) & 0xffff00ff) | ((C) << 8)))

#include "ray_tracer.h"

using namespace std;

const double Object::small_t=1e-6;
//--------------------------------------------------------------------------------
// utility functions
//--------------------------------------------------------------------------------
double sqr(const double x)
{
    return x*x;
}

Pixel Pixel_Color(const Vector_3D<double>& color)
{
    Pixel pixel=0;
    SET_RED(pixel,(unsigned char)(min(color.x,1.0)*255));
    SET_GREEN(pixel,(unsigned char)(min(color.y,1.0)*255));
    SET_BLUE(pixel,(unsigned char)(min(color.z,1.0)*255));
    return pixel;
}

//--------------------------------------------------------------------------------
// Shader
//--------------------------------------------------------------------------------
Vector_3D<double> Phong_Shader::
Shade_Surface(const Ray& ray,const Object& intersection_object,const Vector_3D<double>& intersection_point,const Vector_3D<double>& same_side_normal) const
{

    Vector_3D<double> color;

    for( int i = 0; i < world.lights.size(); i++){
	    Vector_3D<double> light = world.lights[i]->position - intersection_point;
	    light.Normalize();
	    double diffuse = max(0.0, Vector_3D<double>::Dot_Product(same_side_normal, light));
	    
	    Vector_3D<double> h = light - ray.direction;
	    h.Normalize();
	    double spec = max(0.0, Vector_3D<double>::Dot_Product(same_side_normal, h));
	    double specular = pow(spec,specular_power);
	    
	    if( world.enable_shadows ){
	   	Ray shadow;
	    	shadow.endpoint = intersection_point + same_side_normal * intersection_object.small_t; 
	    	shadow.direction = light;
	    
	    	const Object *obj = world.Closest_Intersection(shadow);
	    
	    	if( shadow.semi_infinite == false ) continue;
	    }
	    
	    Vector_3D<double> I = world.lights[i]->Emitted_Light(ray);
	    color = color + color_diffuse * I * diffuse + color_specular * I * specular;

    } 
    
    return color;
}

Vector_3D<double> Reflective_Shader::
Shade_Surface(const Ray& ray,const Object& intersection_object,const Vector_3D<double>& intersection_point,const Vector_3D<double>& same_side_normal) const
{

    Vector_3D<double> color = Phong_Shader::Shade_Surface(ray, intersection_object, intersection_point, same_side_normal);
    
    Ray reflection;
    reflection.endpoint = intersection_point + same_side_normal * intersection_object.small_t;
    reflection.direction = ray.direction - same_side_normal * 2 * Vector_3D<double>::Dot_Product( ray.direction, same_side_normal);

    color = color + world.Cast_Ray(reflection, ray) * reflectivity;

    return color;

}

Vector_3D<double> Flat_Shader::
Shade_Surface(const Ray& ray,const Object& intersection_object,const Vector_3D<double>& intersection_point,const Vector_3D<double>& same_side_normal) const
{
    return color;
}

//--------------------------------------------------------------------------------
// Objects
//--------------------------------------------------------------------------------
// determine if the ray intersects with the sphere
// if there is an intersection, set t_max, current_object, and semi_infinite as appropriate and return true
bool Sphere::
Intersection(Ray& ray) const
{
    // TODO
    
    double a = Vector_3D<double>::Dot_Product(ray.direction, ray.direction);
    double b = 2 * Vector_3D<double>::Dot_Product(ray.direction, ray.endpoint - center);
    double c = Vector_3D<double>::Dot_Product(ray.endpoint - center, ray.endpoint - center) - radius * radius;
    
    double intersection = b * b - 4 * a * c;
    
    if( intersection < 0 ){
    	return false;
    }
    else if( intersection >= 0 ){
  
    	double t = 0;
    	double t_plus =(-b + sqrt( intersection) ) / (2.0 * a);
        double t_minus = (-b - sqrt( intersection) ) / (2.0 * a);
        
    	if( t_plus >= 0 && t_minus >= 0 ) t = min(t_plus, t_minus);
    	else if( t_plus >= 0 && t_minus < 0 ) t = t_plus;
    	else if( t_minus >= 0 && t_plus < 0 ) t = t_minus;
    	
	
    	if( ray.semi_infinite == true && t > small_t ){
    		ray.semi_infinite = false;
    		ray.t_max = t;
    		ray.current_object = this;
    		return true;
    	}
    	else if(ray.semi_infinite == false && t > small_t && t < ray.t_max){
    		ray.t_max = t;
    		ray.current_object = this;
    		return true;
    	}
    	else{
    		return false;
    	}

    }

}

Vector_3D<double> Sphere::
Normal(const Vector_3D<double>& location) const
{

    Vector_3D<double> normal = ( location - center);
    normal.Normalize();
    return normal;
}

// determine if the ray intersects with the sphere
// if there is an intersection, set t_max, current_object, and semi_infinite as appropriate and return true
bool Plane::
Intersection(Ray& ray) const
{

    double numerator = Vector_3D<double>::Dot_Product( normal, x1 - ray.endpoint);
    double denominator = Vector_3D<double>::Dot_Product( normal, ray.direction );
    double t = numerator / denominator;
    
    if( denominator == 0){
    	return false;
    }
    else if( t > small_t && ray.semi_infinite == true){
    	ray.semi_infinite = false;
    	ray.t_max = t;
    	ray.current_object = this;
    	return true;
    }
    else if(t > small_t && t < ray.t_max && ray.semi_infinite == false)
    {
    	ray.t_max = t;
    	ray.current_object = this;
    	return true;
    }
    else{
    	return false;
    }

}

Vector_3D<double> Plane::
Normal(const Vector_3D<double>& location) const
{
    return normal;
}

//--------------------------------------------------------------------------------
// Camera
//--------------------------------------------------------------------------------
// Find the world position of the input pixel
Vector_3D<double> Camera::
World_Position(const Vector_2D<int>& pixel_index)
{
    Vector_2D<double> index = film.pixel_grid.X(pixel_index);
    Vector_3D<double> result = horizontal_vector * index.x + vertical_vector * index.y + focal_point;
    return result;
}

//--------------------------------------------------------------------------------
// Render_World
//--------------------------------------------------------------------------------
// Find the closest object of intersection and return a pointer to it
//   if the ray intersects with an object, then ray.t_max, ray.current_object, and ray.semi_infinite will be set appropriately
//   if there is no intersection do not modify the ray and return 0
const Object* Render_World::
Closest_Intersection(Ray& ray)
{
    // TODO
    for( int i = 0; i < objects.size(); i++)
    {
    	objects[i]->Intersection(ray);
    }
    return ray.current_object;
}

// set up the initial view ray and call 
void Render_World::
Render_Pixel(const Vector_2D<int>& pixel_index)
{
    // TODO
    Ray ray = Ray( camera.position, camera.World_Position(pixel_index) - camera.position);
    Ray dummy_root;
    
    Vector_3D<double> color=Cast_Ray(ray,dummy_root);
    camera.film.Set_Pixel(pixel_index,Pixel_Color(color));
}

// cast ray and return the color of the closest intersected surface point, 
// or the background color if there is no object intersection
Vector_3D<double> Render_World::
Cast_Ray(Ray& ray,const Ray& parent_ray)
{
    // TODO

    Vector_3D<double> color;
    const Object * obj = Closest_Intersection(ray);
  
    if( ray.semi_infinite == false ){ 
    	Vector_3D<double> intersection_point = ray.Point(ray.t_max);
    	Vector_3D<double> normal = obj->Normal(intersection_point);
    	return obj->material_shader->Shade_Surface(ray, *obj, intersection_point, normal);
    }
    else{
        return color;
    }
}
