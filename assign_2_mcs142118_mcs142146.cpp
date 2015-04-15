#include <GL/glut.h>
#include <iostream>
#include <math.h>
#include <vector>

float point_size;

int is_affine;
float sc_x = 1.0, sc_y = 1.0, sc_z = 1.0, sh_xy = 0.0, sh_yx = 0.0;

float boundary_t_inc, inner_t_inc, theta_inc, fi_inc;

float eyex, eyey , eyez, lookAt_x, lookAt_y, lookAt_z, up_x, up_y, up_z;
float near_plane, far_plane;
float width, height;

float ambient, diffuse, specular;

float reflected_red, reflected_green, reflected_blue;
float refracted_red, refracted_green, refracted_blue;
float rec_reflected_red, rec_reflected_green, rec_reflected_blue;
float rec_refracted_red, rec_refracted_green, rec_refracted_blue;

bool is_rec_calculating;

bool is_reflected, is_refracted, is_rec_reflected, is_rec_refracted, in_shadow;

float constant_atten, linear_atten, quadratic_atten;

bool behind = false, right = true;

enum COLOR
{
    RED,
    GREEN,
    BLUE
};

enum SHAPE_TYPE
{
    SPHERE,
    TRIANGLE
};

enum INTERSECTION_TYPE
{
    REFLECTION,
    REFRACTION
};

// these are the identities of the object which have intersected the secondary rays
// and are responsible for recursive ray tracing
int which_recursive_number;
SHAPE_TYPE which_recursive_shape;


// vertex or position vector
class vertex
{
public:
  float x,y,z;

  vertex(){}

  vertex(float x, float y, float z)
  {
    this->x = x;
    this->y = y;
    this->z = z;
  }

  vertex operator-( vertex v )
  {
    vertex w;
    w.x = x - v.x;
    w.y = y - v.y;
    w.z = z - v.z;
    return w;
  }

  vertex operator+( vertex v )
  {
    vertex w;
    w.x = x + v.x;
    w.y = y + v.y;
    w.z = z + v.z;
    return w;
  }

  // scalar product or dot product
  float operator*( vertex &v )
  {
    return x*v.x + y*v.y + z*v.z;
  }

  // multiplication with scalar
  vertex operator *( float i )
  {
      vertex temp;
      temp.x = x*i;
      temp.y = y*i;
      temp.z = z*i;
      return temp;
  }

};

class light
{
public:
    vertex position;
    float intensity;
    light()
    {

    }
    light( vertex position, float intensity)
    {
        this->position = position;
        this->intensity = intensity;
    }
};

// distance between 2 vertices
float distance(vertex v1, vertex v2)
{
    return sqrt( pow((v1.x-v2.x),2) +pow((v1.y-v2.y),2) +pow((v1.z-v2.z),2) );
}

// x2 + y2 + z2 = 1 and (x,y,z) is direction
class direction_unit_vector
{
public:
  float x,y,z;

  direction_unit_vector(){}

  // line from u to v is the direction and we are computing unit vector in the same direction
  direction_unit_vector( vertex u, vertex v)
  {
      float dist = distance( u, v);

      x = v.x - u.x;
      y = v.y - u.y;
      z = v.z - u.z;

      x = x / dist;
      y = y / dist;
      z = z / dist;
  }

  // line from u to v is the direction and we are computing unit vector in the same direction
  direction_unit_vector( vertex *u, vertex *v)
  {
      float dist = distance( *u, *v);

      x = v->x - u->x;
      y = v->y - u->y;
      z = v->z - u->z;

      x = x / dist;
      y = y / dist;
      z = z / dist;
  }

  // u -> v is the direction, we are dividing it with its length to get unit vector
  void unit_vectorize()
  {
    vertex u = vertex(0.0 ,0.0 ,0.0);
    vertex v = vertex(x   ,y   ,z  );
    float dist = distance( u, v);
    x /= dist;
    y /= dist;
    z /= dist;
  }

  direction_unit_vector operator-( direction_unit_vector &v )
  {
    direction_unit_vector w;
    w.x = x - v.x;
    w.y = y - v.y;
    w.z = z - v.z;
    return w;
  }

  direction_unit_vector operator+( direction_unit_vector &v )
  {
    direction_unit_vector w;
    w.x = x + v.x;
    w.y = y + v.y;
    w.z = z + v.z;
    return w;
  }

  // multiply with scalar
  direction_unit_vector operator*( float i )
  {
    direction_unit_vector w;
    w.x = x * i;
    w.y = y * i;
    w.z = z * i;
    return w;
  }

  // inverse the direction
  void operator-()
  {
      x = -x;
      y = -y;
      z = -z;
  }

  float check_mag()
  {
      return x*x + y*y + z*z;
  }
};

bool check_intersection( vertex surface_position, SHAPE_TYPE s, int i, direction_unit_vector transmit_dir, enum INTERSECTION_TYPE in );

class sphere
{
    public:

    float radius;
    vertex center;

    float cr, cg, cb;
    float ka, kd, ks, n;
    float refl_cof, refr_cof, refr_ind;

    sphere()
    {

    }
    sphere( float radius, vertex center, float cr, float cg, float cb, float ka, float kd, float ks, float n, float refl_cof, float refr_cof, float refr_ind )
    {
        this->radius = radius;
        this->center = center;
        this->cr = cr;
        this->cg = cg;
        this->cb = cb;
        this->ka = ka;
        this->kd = kd;
        this->ks = ks;
        this->n = n;
        this->refl_cof = refl_cof;
        this->refr_cof = refr_cof;
        this->refr_ind = refr_ind;
    }

};

class triangle
{
    public:

    vertex v1, v2, v3;
    float eqa,eqb,eqc,eqd;

    direction_unit_vector surface_normal;

    float cr, cg, cb;
    float ka, kd, ks, n;
    float refl_cof, refr_cof, refr_ind;

    triangle() {}

    triangle( float eqa, float eqb, float eqc, float eqd, direction_unit_vector surface_normal, vertex *v1, vertex *v2, vertex *v3, float cr, float cg, float cb, float ka, float kd, float ks, float n, float refl_cof, float refr_cof, float refr_ind )
    {
        this->eqa = eqa;
        this->eqb = eqb;
        this->eqc = eqc;
        this->eqd = eqd;
        this->surface_normal = surface_normal;
        this->v1 = *v1;
        this->v2 = *v2;
        this->v3 = *v3;
        this->cr = cr;
        this->cg = cg;
        this->cb = cb;
        this->ka = ka;
        this->kd = kd;
        this->ks = ks;
        this->n = n;
        this->refl_cof = refl_cof;
        this->refr_cof = refr_cof;
        this->refr_ind = refr_ind;
    }
};


// u.v cos (theta) where theta is the angle between 2 directions
float dot_product( direction_unit_vector u, direction_unit_vector v)
{
    return u.x*v.x + u.y*v.y + u.z*v.z;
}

direction_unit_vector cross_product( direction_unit_vector u, direction_unit_vector v )
{
  direction_unit_vector w;
  w.x = u.y*v.z - u.z*v.y;
  w.y = u.z*v.x - u.x*v.z;
  w.z = u.x*v.y - u.y*v.x;
  return w;
}

vertex cross_product_vertex(vertex u, vertex v)
{
  vertex w;
  w.x = u.y*v.z - u.z*v.y;
  w.y = u.z*v.x - u.x*v.z;
  w.z = u.x*v.y - u.y*v.x;
  return w;

}

int no_of_spheres, no_of_triangles, no_of_light_source;
std::vector<sphere> spheres;
std::vector<triangle> triangles;
std::vector<light> lights;

// returns false if obect is obstructed by other object
bool is_visible( vertex surface_position, SHAPE_TYPE s, int i )
{
    vertex *eye = new vertex( eyex, eyey, eyez);
    direction_unit_vector *dir = new direction_unit_vector( surface_position, *eye);

    vertex start = surface_position;
    vertex *direction = new vertex( dir->x, dir->y, dir->z);

    float t, eye_t;

    // this is the direction of ray from surface to eye
    vertex ray = start + *direction * t;

    // any object with t < eye_t is an obstacle and will make this surface under non visible
    eye_t = (*eye-start).x / direction->x;

    // checking for each sphere, method from sir's slides
    for( int j = 0; j < no_of_spheres; j++ )
    {
        if( s == SPHERE && i == j ) continue;

        float local_t;

        float A = pow( dir->x, 2) + pow( dir->y, 2) + pow( dir->z, 2);
        float B = 2 * ( dir->x*(start.x - spheres.at(j).center.x) + dir->y*(start.y - spheres.at(j).center.y) + dir->z*(start.z - spheres.at(j).center.z) );
        float C = pow( start.x - spheres.at(j).center.x, 2) + pow( start.y - spheres.at(j).center.y, 2) + pow( start.z - spheres.at(j).center.z, 2) - pow( spheres.at(j).radius, 2);
        float D = B*B - 4*A*C;

        if( D < 0 ) continue;

        float t0 = (-B - sqrt(D))/2*A;
        float t1 = (-B + sqrt(D))/2*A;
        if( t0 < t1 ) local_t = t0; else local_t = t1;

        // this sphere is in between surface and eye
        if( local_t > 0 && local_t < eye_t ) return false;

    }

    // checking for each triangle, method from sir's slides
    for( int j = 0; j < no_of_triangles; j++ )
    {
        if( s == TRIANGLE && i == j ) continue;

        float local_t;

        direction_unit_vector pn = triangles.at(j).surface_normal;
        float eqa = triangles.at(j).eqa;
        float eqb = triangles.at(j).eqb;
        float eqc = triangles.at(j).eqc;
        float eqd = triangles.at(j).eqd;

        if( eqa*dir->x + eqb*dir->y + eqc*dir->z == 0 ) continue;

        local_t =  eqa*start.x + eqb*start.y + eqc*start.z + eqd;
        local_t /= eqa*dir->x + eqb*dir->y + eqc*dir->z;
        local_t = -local_t;

        direction_unit_vector surface_normal = pn;

        // ray is the point which is possibly obstructing start and eye
        ray = start + *direction * local_t;

        // containment test
        vertex v1 = triangles.at(j).v1, v2 = triangles.at(j).v2, v3 = triangles.at(j).v3;

        vertex cp1 = cross_product_vertex( v1-v2,  v3-v2);
        vertex cp2 = cross_product_vertex( v1-v2, ray-v2);

        if( cp1*cp2 < 0 ) continue;

        vertex cp3 = cross_product_vertex( v2-v3,  v1-v3);
        vertex cp4 = cross_product_vertex( v2-v3, ray-v3);

        if( cp3*cp4 < 0 ) continue;

        vertex cp5 = cross_product_vertex( v3-v1,  v2-v1);
        vertex cp6 = cross_product_vertex( v3-v1, ray-v1);

        if( cp5*cp6 < 0 ) continue;

        // this triangle is in between surface and eye
        if( local_t > 0 && local_t < eye_t ) return false;
    }

    return true;

}

bool is_shadow( vertex surface_position, vertex light_source_position, SHAPE_TYPE s, int i )
{
    // unit vector from surface to light source
    direction_unit_vector *dir = new direction_unit_vector(surface_position, light_source_position);

    vertex start = surface_position;
    vertex *direction = new vertex( dir->x, dir->y, dir->z);

    float t, light_t;

    // this is the direction of ray from surface to light source
    vertex ray = start + *direction * t;

    // any object with t < light_t is an obstacle and will make this surface under shadow
    light_t = (light_source_position-start).x / direction->x;

    // checking for each sphere, method from sir's slides
    for( int j = 0; j < no_of_spheres; j++ )
    {
        if( s == SPHERE && i == j ) continue;

        float local_t;

        float A = pow( dir->x, 2) + pow( dir->y, 2) + pow( dir->z, 2);
        float B = 2 * ( dir->x*(start.x - spheres.at(j).center.x) + dir->y*(start.y - spheres.at(j).center.y) + dir->z*(start.z - spheres.at(j).center.z) );
        float C = pow( start.x - spheres.at(j).center.x, 2) + pow( start.y - spheres.at(j).center.y, 2) + pow( start.z - spheres.at(j).center.z, 2) - pow( spheres.at(j).radius, 2);
        float D = B*B - 4*A*C;

        if( D < 0 ) continue;

        float t0 = (-B - sqrt(D))/2*A;
        float t1 = (-B + sqrt(D))/2*A;
        if( t0 < t1 ) local_t = t0; else local_t = t1;

        // this sphere is in between surface and light source
        if( local_t > 0 && local_t < light_t ) return true;

    }

    // checking for each triangle, method from sir's slides
    for( int j = 0; j < no_of_triangles; j++ )
    {
        if( s == TRIANGLE && i == j ) continue;

        float local_t;

        direction_unit_vector pn = triangles.at(j).surface_normal;
        float eqa = triangles.at(j).eqa;
        float eqb = triangles.at(j).eqb;
        float eqc = triangles.at(j).eqc;
        float eqd = triangles.at(j).eqd;

        if( eqa*dir->x + eqb*dir->y + eqc*dir->z == 0 ) continue;

        local_t =  eqa*start.x + eqb*start.y + eqc*start.z + eqd;
        local_t /= eqa*dir->x + eqb*dir->y + eqc*dir->z;
        local_t = -local_t;

        direction_unit_vector surface_normal = pn;

        // ray is the point which is possibly obstructing start and light source
        ray = start + *direction * local_t;

        // containment test
        vertex v1 = triangles.at(j).v1, v2 = triangles.at(j).v2, v3 = triangles.at(j).v3;

        vertex cp1 = cross_product_vertex( v1-v2,  v3-v2);
        vertex cp2 = cross_product_vertex( v1-v2, ray-v2);

        if( cp1*cp2 < 0 ) continue;

        vertex cp3 = cross_product_vertex( v2-v3,  v1-v3);
        vertex cp4 = cross_product_vertex( v2-v3, ray-v3);

        if( cp3*cp4 < 0 ) continue;

        vertex cp5 = cross_product_vertex( v3-v1,  v2-v1);
        vertex cp6 = cross_product_vertex( v3-v1, ray-v1);

        if( cp5*cp6 < 0 ) continue;

        // this triangle is in between surface and light source
        if( local_t > 0 && local_t < light_t ) return true;
    }

    return false;
}

// last argument stands for -> from original surface to point which is to be reflected/refracted, it is to find the recursive reflected/refracted ray
// last argument stands useless if "is_rec_calculating" is TRUE because there we don't have to find the reflection/refraction
float calculate_light_intensity( COLOR c, SHAPE_TYPE s, int i, vertex surface_position, direction_unit_vector surface_normal, direction_unit_vector view_dir, direction_unit_vector transmitted_direction )
{
    // check for recursive calculation if not doing it
    if( !is_rec_calculating )
    {
        is_rec_calculating = true;

        // this is from point to be reflected to original surface
        direction_unit_vector inverse_transmitted_direction = transmitted_direction;
        -inverse_transmitted_direction;

        if( s == TRIANGLE )
        {

            // original direction is itself the direction for refraction
            direction_unit_vector refract_dir = transmitted_direction;
            if( triangles.at(i).refr_cof > 0.01 ) check_intersection( surface_position, s, i, refract_dir, REFRACTION );

            // this is the reflection direction
            direction_unit_vector reflect_dir;

            direction_unit_vector temp = surface_normal * (2 * dot_product( inverse_transmitted_direction, surface_normal));
            reflect_dir = temp - inverse_transmitted_direction;

            if( triangles.at(i).refl_cof > 0.01 ) check_intersection( surface_position, s, i, reflect_dir, REFLECTION );

        }

        if( s == SPHERE )
        {

            // REFRACTION THROUGH SPHERE - complicated one
            // this is the refraction direction

            // to find direction of refracted ray, method from sir'slides
            float cos_theta_i    = dot_product( inverse_transmitted_direction, surface_normal);
            float sin_theta_i    = 1 - pow(cos_theta_i, 2);
            float refr_ind_ratio = 1/(spheres.at(i).refr_ind);
            float sin_theta_t    = refr_ind_ratio * sin_theta_i;
            float cos_theta_t    = 1 - pow(sin_theta_t, 2);

            direction_unit_vector temp_I = transmitted_direction;

            vertex *I = new vertex( temp_I.x         , temp_I.y         , temp_I.z          );
            vertex *N = new vertex( surface_normal.x, surface_normal.y, surface_normal.z );

            // this is the refraction direction from surface to in the direction of refracted ray
            vertex *direction = new vertex();
            *direction = ( *I + *N * cos_theta_i ) * refr_ind_ratio - *N * cos_theta_t;

            direction_unit_vector *refract_dir = new direction_unit_vector( new vertex(0.0,0.0,0.0), direction);

            // t is at infinity i.e. the object to be refracted is at infinity
            float t = 1000;
            vertex start = surface_position;
            direction = new vertex( refract_dir->x, refract_dir->y, refract_dir->z);

            // this is the direction of ray from (tx,ty,tz) in the direction of refraction i.e. going inside the sphere
            vertex ray = start + *direction * t;

            // now, we have to find the position from where this ray will come out of the sphere
            float local_t;

            float A = pow( refract_dir->x, 2) + pow( refract_dir->y, 2) + pow( refract_dir->z, 2);
            float B = 2 * ( refract_dir->x*(start.x - spheres.at(i).center.x) + refract_dir->y*(start.y - spheres.at(i).center.y) + refract_dir->z*(start.z - spheres.at(i).center.z) );
            float C = pow( start.x - spheres.at(i).center.x, 2) + pow( start.y - spheres.at(i).center.y, 2) + pow( start.z - spheres.at(i).center.z, 2) - pow( spheres.at(i).radius, 2);
            float D = B*B - 4*A*C;

            float t0 = (-B - sqrt(D))/2*A;
            float t1 = (-B + sqrt(D))/2*A;
            if( t0 < t1 ) local_t = t1; else local_t = t0;

            t = local_t;

            // ray is the point which is refracted at (tx, ty, tz) and intensity at (ray.x, ray.y, ray.z) will be computed now
            ray = start + *direction * t;

            vertex *cen = new vertex( spheres.at(i).center.x, spheres.at(i).center.y, spheres.at(i).center.z);
            direction_unit_vector *surface_normal_1 = new direction_unit_vector(ray, *cen);


            direction_unit_vector *view_dir_1 = new direction_unit_vector( ray, surface_position);

            // to find direction of refracted ray, method from sir'slides
            cos_theta_i    = dot_product( *view_dir_1, *surface_normal_1);
            sin_theta_i    = 1 - pow(cos_theta_i, 2);
            refr_ind_ratio = spheres.at(i).refr_ind;
            sin_theta_t    = refr_ind_ratio * sin_theta_i;
            cos_theta_t    = 1 - pow(sin_theta_t, 2);

            temp_I = *view_dir_1;
            -temp_I;

            I = new vertex( temp_I.x         , temp_I.y         , temp_I.z          );
            N = new vertex( surface_normal_1->x, surface_normal_1->y, surface_normal_1->z );

            // this is the refraction direction from surface to in the direction of refracted ray
            *direction = ( *I + *N * cos_theta_i ) * refr_ind_ratio - *N * cos_theta_t;

            refract_dir = new direction_unit_vector( new vertex(0.0,0.0,0.0), direction);

            if( spheres.at(i).refr_cof > 0.01 ) check_intersection( surface_position, s, i, *refract_dir, REFRACTION );

            // this is the reflection direction
            direction_unit_vector reflect_dir;

            direction_unit_vector temp = surface_normal * (2 * dot_product( inverse_transmitted_direction, surface_normal));
            reflect_dir = temp - inverse_transmitted_direction;

            if( spheres.at(i).refl_cof > 0.01 ) check_intersection( surface_position, s, i, reflect_dir, REFLECTION );

        }

        is_rec_calculating = false;
    }

    // NOW, go for the local light sources
    float temp_diffuse = 0.0, temp_specular = 0.0;
    for( int j = 0; j < no_of_light_source; j++ )
    {

        // diffuse and specular contributed by this particular light source
        float local_diffuse = 0.0, local_specular = 0.0;

        // light source position
        vertex *light_source_position = new vertex( lights.at( j ).position.x, lights.at( j ).position.y, lights.at( j ).position.z);

         // light direction unit vector is negated to match sir's conventions
        direction_unit_vector *light_dir = new direction_unit_vector( *light_source_position, surface_position);
        -(*light_dir);

        // SHADOW, if in shadow, check for the next light source

        if( is_shadow(surface_position, *light_source_position, s, i) ) continue;

        // DIFFUSE

        if( s == SPHERE )
        {
            float diffuse_dot_prod = dot_product( *light_dir, surface_normal);

            // if (tx,ty,tz) is not visible by the light source, no need to proceed for this light source
            if( diffuse_dot_prod < 0 ) continue;

            local_diffuse = spheres.at( i ).kd * lights.at( j ).intensity * diffuse_dot_prod;
        }
        else
        {
            float diffuse_dot_prod;

            if( dot_product( view_dir, surface_normal) > 0 && dot_product( *light_dir, surface_normal) > 0 )
            {
                diffuse_dot_prod = dot_product( *light_dir, surface_normal);
                local_diffuse = triangles.at( i ).kd * lights.at( j ).intensity * diffuse_dot_prod;
            }
            else if( dot_product( view_dir, surface_normal) < 0 && dot_product( *light_dir, surface_normal) < 0 )
            {
                diffuse_dot_prod = -dot_product( *light_dir, surface_normal);
                local_diffuse = triangles.at( i ).kd * lights.at( j ).intensity * diffuse_dot_prod;
            }
            else
                continue;
        }

        // SPECULAR

        // compute Half Way Vector H = (L+V)/|L+V|, NOTE half_way_vector is not a pointer
        direction_unit_vector half_way_vector = *light_dir + view_dir;
        half_way_vector.unit_vectorize();

        if( s == SPHERE )
        {
            float specular_dot_prod = pow ( dot_product( surface_normal, half_way_vector), spheres.at( i ).n);
            local_specular += spheres.at( i ).ks * lights.at( j ).intensity * specular_dot_prod;
        }
        else
        {
            float specular_dot_prod = pow ( dot_product( surface_normal, half_way_vector), triangles.at( i ).n);
            local_specular += triangles.at( i ).ks * lights.at( j ).intensity * specular_dot_prod;
        }
        //ATTENUATE THE INTENSITIES

        local_diffuse  /= constant_atten + (linear_atten)*distance(*light_source_position, surface_position) + (quadratic_atten)*distance(*light_source_position, surface_position)*distance(*light_source_position, surface_position);
        local_specular /= constant_atten + (linear_atten)*distance(*light_source_position, surface_position) + (quadratic_atten)*distance(*light_source_position, surface_position)*distance(*light_source_position, surface_position);

        temp_diffuse += local_diffuse;
        temp_specular += local_specular;
    }

    // if TRUE, it means no need to check rec prefixed booleans as it will mess up the values to be returned
    if( is_rec_calculating )
    {
        if( s == SPHERE )
        {
            if( c == RED   ) return (ambient * spheres.at( i ).ka + temp_diffuse) * spheres.at( i ).cr + temp_specular;
            if( c == GREEN ) return (ambient * spheres.at( i ).ka + temp_diffuse) * spheres.at( i ).cg + temp_specular;
            if( c == BLUE  ) return (ambient * spheres.at( i ).ka + temp_diffuse) * spheres.at( i ).cb + temp_specular;
        }
        else
        {
            if( c == RED   ) return (ambient * triangles.at( i ).ka + temp_diffuse) * triangles.at( i ).cr + temp_specular;
            if( c == GREEN ) return (ambient * triangles.at( i ).ka + temp_diffuse) * triangles.at( i ).cg + temp_specular;
            if( c == BLUE  ) return (ambient * triangles.at( i ).ka + temp_diffuse) * triangles.at( i ).cb + temp_specular;
        }
    }

    if( is_rec_reflected )
    {
        if( is_rec_refracted )
        {
            if( s == SPHERE )
            {
                if( spheres.at(i).refl_cof + spheres.at(i).refr_cof > 1 )
                {
                    spheres.at(i).refl_cof /= spheres.at(i).refl_cof + spheres.at(i).refr_cof;
                    spheres.at(i).refr_cof /= spheres.at(i).refl_cof + spheres.at(i).refr_cof;
                }

                if( c == RED   ) return (ambient * spheres.at( i ).ka + temp_diffuse) * spheres.at( i ).cr * (1 - spheres.at(i).refl_cof - spheres.at(i).refr_cof ) + temp_specular + rec_reflected_red   * ( spheres.at(i).refl_cof ) + rec_refracted_red   * ( spheres.at(i).refr_cof );
                if( c == GREEN ) return (ambient * spheres.at( i ).ka + temp_diffuse) * spheres.at( i ).cg * (1 - spheres.at(i).refl_cof - spheres.at(i).refr_cof ) + temp_specular + rec_reflected_green * ( spheres.at(i).refl_cof ) + rec_refracted_green * ( spheres.at(i).refr_cof );
                if( c == BLUE  ) return (ambient * spheres.at( i ).ka + temp_diffuse) * spheres.at( i ).cb * (1 - spheres.at(i).refl_cof - spheres.at(i).refr_cof ) + temp_specular + rec_reflected_blue  * ( spheres.at(i).refl_cof ) + rec_refracted_blue  * ( spheres.at(i).refr_cof );
            }
            else
            {
                if( triangles.at(i).refl_cof + triangles.at(i).refr_cof > 1 )
                {
                    triangles.at(i).refl_cof /= triangles.at(i).refl_cof + triangles.at(i).refr_cof;
                    triangles.at(i).refr_cof /= triangles.at(i).refl_cof + triangles.at(i).refr_cof;
                }
                if( c == RED   ) return (ambient * triangles.at( i ).ka + temp_diffuse) * triangles.at( i ).cr * (1 - triangles.at(i).refl_cof - triangles.at(i).refr_cof ) + temp_specular + rec_reflected_red   * ( triangles.at(i).refl_cof ) + rec_refracted_red   * ( triangles.at(i).refr_cof );
                if( c == GREEN ) return (ambient * triangles.at( i ).ka + temp_diffuse) * triangles.at( i ).cg * (1 - triangles.at(i).refl_cof - triangles.at(i).refr_cof ) + temp_specular + rec_reflected_green * ( triangles.at(i).refl_cof ) + rec_refracted_green * ( triangles.at(i).refr_cof );
                if( c == BLUE  ) return (ambient * triangles.at( i ).ka + temp_diffuse) * triangles.at( i ).cb * (1 - triangles.at(i).refl_cof - triangles.at(i).refr_cof ) + temp_specular + rec_reflected_blue  * ( triangles.at(i).refl_cof ) + rec_refracted_blue  * ( triangles.at(i).refr_cof );

            }
        }
        else
        {

            if( s == SPHERE )
            {
                if( c == RED   ) return (ambient * spheres.at( i ).ka + temp_diffuse) * spheres.at( i ).cr * (1 - spheres.at(i).refl_cof ) + temp_specular + rec_reflected_red   * ( spheres.at(i).refl_cof );
                if( c == GREEN ) return (ambient * spheres.at( i ).ka + temp_diffuse) * spheres.at( i ).cg * (1 - spheres.at(i).refl_cof ) + temp_specular + rec_reflected_green * ( spheres.at(i).refl_cof );
                if( c == BLUE  ) return (ambient * spheres.at( i ).ka + temp_diffuse) * spheres.at( i ).cb * (1 - spheres.at(i).refl_cof ) + temp_specular + rec_reflected_blue  * ( spheres.at(i).refl_cof );
            }
            else
            {
                if( c == RED   ) return (ambient * triangles.at( i ).ka + temp_diffuse) * triangles.at( i ).cr * (1 - triangles.at(i).refl_cof ) + temp_specular + rec_reflected_red   * ( triangles.at(i).refl_cof );
                if( c == GREEN ) return (ambient * triangles.at( i ).ka + temp_diffuse) * triangles.at( i ).cg * (1 - triangles.at(i).refl_cof ) + temp_specular + rec_reflected_green * ( triangles.at(i).refl_cof );
                if( c == BLUE  ) return (ambient * triangles.at( i ).ka + temp_diffuse) * triangles.at( i ).cb * (1 - triangles.at(i).refl_cof ) + temp_specular + rec_reflected_blue  * ( triangles.at(i).refl_cof );
            }
        }
    }
    else if( is_rec_refracted )
    {
        if( s == SPHERE )
        {
            if( c == RED   ) return (ambient * spheres.at( i ).ka + temp_diffuse) * spheres.at( i ).cr * (1 - spheres.at(i).refr_cof ) + temp_specular + rec_refracted_red   * ( spheres.at(i).refr_cof );
            if( c == GREEN ) return (ambient * spheres.at( i ).ka + temp_diffuse) * spheres.at( i ).cg * (1 - spheres.at(i).refr_cof ) + temp_specular + rec_refracted_green * ( spheres.at(i).refr_cof );
            if( c == BLUE  ) return (ambient * spheres.at( i ).ka + temp_diffuse) * spheres.at( i ).cb * (1 - spheres.at(i).refr_cof ) + temp_specular + rec_refracted_blue  * ( spheres.at(i).refr_cof );
        }
        else
        {
            if( c == RED   ) return (ambient * triangles.at( i ).ka + temp_diffuse) * triangles.at( i ).cr * (1 - triangles.at(i).refr_cof ) + temp_specular + rec_refracted_red   * ( triangles.at(i).refr_cof );
            if( c == GREEN ) return (ambient * triangles.at( i ).ka + temp_diffuse) * triangles.at( i ).cg * (1 - triangles.at(i).refr_cof ) + temp_specular + rec_refracted_green * ( triangles.at(i).refr_cof );
            if( c == BLUE  ) return (ambient * triangles.at( i ).ka + temp_diffuse) * triangles.at( i ).cb * (1 - triangles.at(i).refr_cof ) + temp_specular + rec_refracted_blue  * ( triangles.at(i).refr_cof );
        }
    }
    else
    {
        if( s == SPHERE )
        {
            if( c == RED   ) return (ambient * spheres.at( i ).ka + temp_diffuse) * spheres.at( i ).cr + temp_specular;
            if( c == GREEN ) return (ambient * spheres.at( i ).ka + temp_diffuse) * spheres.at( i ).cg + temp_specular;
            if( c == BLUE  ) return (ambient * spheres.at( i ).ka + temp_diffuse) * spheres.at( i ).cb + temp_specular;
        }
        else
        {
            if( c == RED   ) return (ambient * triangles.at( i ).ka + temp_diffuse) * triangles.at( i ).cr + temp_specular;
            if( c == GREEN ) return (ambient * triangles.at( i ).ka + temp_diffuse) * triangles.at( i ).cg + temp_specular;
            if( c == BLUE  ) return (ambient * triangles.at( i ).ka + temp_diffuse) * triangles.at( i ).cb + temp_specular;
        }
    }
}

bool check_intersection( vertex surface_position, SHAPE_TYPE s, int i, direction_unit_vector transmit_dir, enum INTERSECTION_TYPE in )
{

    vertex start = surface_position;
    vertex *eye_position = new vertex(eyex, eyey, eyez);

    vertex *direction = new vertex( transmit_dir.x, transmit_dir.y, transmit_dir.z );

    float t = 1000;
    vertex ray = start + *direction * t;

    // checking for each sphere, method from sir's slides
    for( int j = 0; j < no_of_spheres; j++ )
    {
        if( s == SPHERE && i == j ) continue;

        float local_t;

        float A = pow( transmit_dir.x, 2) + pow( transmit_dir.y, 2) + pow( transmit_dir.z, 2);
        float B = 2 * ( transmit_dir.x*(start.x - spheres.at(j).center.x) + transmit_dir.y*(start.y - spheres.at(j).center.y) + transmit_dir.z*(start.z - spheres.at(j).center.z) );
        float C = pow( start.x - spheres.at(j).center.x, 2) + pow( start.y - spheres.at(j).center.y, 2) + pow( start.z - spheres.at(j).center.z, 2) - pow( spheres.at(j).radius, 2);
        float D = B*B - 4*A*C;

        if( D < 0 ) continue;

        float t0 = (-B - sqrt(D))/2*A;
        float t1 = (-B + sqrt(D))/2*A;
        if( t0 < t1 ) local_t = t0; else local_t = t1;

        if( local_t < 0 ) continue;

        // farther from some other object, will be rejected
        if( local_t > t ) continue;

        t = local_t;

        ray = start + *direction * t;

        vertex *cen = new vertex( spheres.at(j).center.x, spheres.at(j).center.y, spheres.at(j).center.z);
        direction_unit_vector *surface_normal = new direction_unit_vector(*cen, ray);

        // this is from surface to eye to match sir's convention
        direction_unit_vector *view_dir = new direction_unit_vector( ray, *eye_position);

        is_rec_calculating = true;

        if( in == REFLECTION )
        {
            rec_reflected_red   = calculate_light_intensity( RED,   SPHERE, j, ray, *surface_normal, *view_dir, *surface_normal );
            rec_reflected_green = calculate_light_intensity( GREEN, SPHERE, j, ray, *surface_normal, *view_dir, *surface_normal );
            rec_reflected_blue  = calculate_light_intensity( BLUE,  SPHERE, j, ray, *surface_normal, *view_dir, *surface_normal );

            is_rec_reflected = true;
        }
        else
        {
            rec_refracted_red   = calculate_light_intensity( RED,   SPHERE, j, ray, *surface_normal, *view_dir, *surface_normal );
            rec_refracted_green = calculate_light_intensity( GREEN, SPHERE, j, ray, *surface_normal, *view_dir, *surface_normal );
            rec_refracted_blue  = calculate_light_intensity( BLUE,  SPHERE, j, ray, *surface_normal, *view_dir, *surface_normal );

            is_rec_refracted = true;
        }
        which_recursive_number = j;
        which_recursive_shape = SPHERE;
    }

    // checking for each triangle, method from sir's slides
    for( int j = 0; j < no_of_triangles; j++ )
    {
        if( s == TRIANGLE && i == j ) continue;

        float local_t;

        direction_unit_vector pn = triangles.at(j).surface_normal;
        float eqa = triangles.at(j).eqa;
        float eqb = triangles.at(j).eqb;
        float eqc = triangles.at(j).eqc;
        float eqd = triangles.at(j).eqd;

        if( eqa*transmit_dir.x + eqb*transmit_dir.y + eqc*transmit_dir.z == 0 ) continue;

        local_t =  eqa*start.x + eqb*start.y + eqc*start.z + eqd;
        local_t /= eqa*transmit_dir.x + eqb*transmit_dir.y + eqc*transmit_dir.z;
        local_t = -local_t;

        if( local_t < 0 ) continue;

        // farther from some other object, will be rejected
        if( local_t > t ) continue;

        direction_unit_vector surface_normal = pn;

        ray = start + *direction * local_t;

        // containment test
        vertex v1 = triangles.at(j).v1, v2 = triangles.at(j).v2, v3 = triangles.at(j).v3;

        vertex cp1 = cross_product_vertex( v1-v2,  v3-v2);
        vertex cp2 = cross_product_vertex( v1-v2, ray-v2);

        if( cp1*cp2 < 0 ) continue;

        vertex cp3 = cross_product_vertex( v2-v3,  v1-v3);
        vertex cp4 = cross_product_vertex( v2-v3, ray-v3);

        if( cp3*cp4 < 0 ) continue;

        vertex cp5 = cross_product_vertex( v3-v1,  v2-v1);
        vertex cp6 = cross_product_vertex( v3-v1, ray-v1);

        if( cp5*cp6 < 0 ) continue;

        t = local_t;

        // this is from surface to eye to match sir's convention
        direction_unit_vector *view_dir = new direction_unit_vector( ray, *eye_position);

        is_rec_calculating = true;

        if( in == REFLECTION )
        {
            rec_reflected_red   = calculate_light_intensity( RED,   TRIANGLE, j, ray,  surface_normal, *view_dir,  surface_normal );
            rec_reflected_green = calculate_light_intensity( GREEN, TRIANGLE, j, ray,  surface_normal, *view_dir,  surface_normal );
            rec_reflected_blue  = calculate_light_intensity( BLUE,  TRIANGLE, j, ray,  surface_normal, *view_dir,  surface_normal );

            is_rec_reflected = true;
        }
        else
        {
            rec_refracted_red   = calculate_light_intensity( RED,   TRIANGLE, j, ray,  surface_normal, *view_dir,  surface_normal );
            rec_refracted_green = calculate_light_intensity( GREEN, TRIANGLE, j, ray,  surface_normal, *view_dir,  surface_normal );
            rec_refracted_blue  = calculate_light_intensity( BLUE,  TRIANGLE, j, ray,  surface_normal, *view_dir,  surface_normal );

            is_rec_refracted = true;
        }

        which_recursive_number = j;
        which_recursive_shape = TRIANGLE;
    }

}

void init(void)
{
  glEnable(GL_DEPTH_TEST);
  glClearColor (0.0, 0.0, 0.0, 0.0);
  glShadeModel (GL_FLAT);
  glPointSize( point_size );
  glLineWidth( 25.0 );
}

void display(void)
{
  glClear (GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
  glLoadIdentity();
  gluLookAt(eyex, eyey, eyez, lookAt_x, lookAt_y, lookAt_z, up_x, up_y, up_z);

  // DRAW COORDINATE AXES

  glColor3f( 0.3, 0.3, 0.3);
  glBegin(GL_LINES);
    glVertex3f(0.0, 0.0, 0.0);
    glVertex3f(4.0, 0.0, 0.0);
    glVertex3f(0.0, 0.0, 0.0);
    glVertex3f(0.0, 4.0, 0.0);
    glVertex3f(0.0, 0.0, 0.0);
    glVertex3f(0.0, 0.0, 4.0);
  glEnd();

  // DRAW LIGHT BULBS
  for( int i = 0; i < no_of_light_source; i++ )
  {
    glPushMatrix();
    glTranslatef( lights.at(i).position.x, lights.at(i).position.y, lights.at(i).position.z);
    glColor3f( 1.0, 1.0, 1.0);
    glutSolidSphere( 0.3, 10, 10 );
    glPopMatrix();
  }

  // DRAW SPHERES
  for( int i = 0; i < no_of_spheres; i++ )
  {
    std::cout << "S" << i << "\n";

    float tx , ty, tz;
    float theta, fi;

    // center and radius of sphere
    vertex *cen = new vertex( spheres.at( i ).center.x, spheres.at( i ).center.y, spheres.at( i ).center.z);
    float radius = spheres.at( i ).radius;

    glPushMatrix();
    glBegin(GL_POINTS);

    int count = 0;
    for( theta = 0; theta < 180; theta+= theta_inc )
    {
        for( fi = 0; fi < 360; fi+= fi_inc )
        {
            // sphere vertex at this point is tx,ty,tz centered at ( cen->x, cen->y, cen->z )
            tx =  cen->x + radius * sin(theta) * cos( fi );
            ty =  cen->y + radius * sin(theta) * sin( fi );
            tz =  cen->z + radius * cos(theta);

            // find normal to this point
            vertex *surface_position = new vertex( tx, ty, tz);

            // this point is obstructed by some other object, hence not visible
            if( !is_visible( *surface_position, SPHERE, i) ) continue;

            direction_unit_vector *surface_normal = new direction_unit_vector(*cen, *surface_position);

            // this is eye position and direction from surface to eye to match sir's convention
            vertex *eye_position = new vertex( eyex, eyey, eyez);
            direction_unit_vector *view_dir = new direction_unit_vector( *surface_position, *eye_position);

            // (tx,ty,tz) is not visible from eye, hence no need to draw this vertex
            if( dot_product( *view_dir, *surface_normal) < 0 ) continue;

            // initially, no vertex is illuminated either by diffuse or specular intensity at ( tx, ty, tz)
            // initially, no vertex is reflecting any other object
            diffuse         = 0.0;
            specular        = 0.0;
            reflected_red   = 0.0;
            reflected_green = 0.0;
            reflected_blue  = 0.0;
            refracted_red   = 0.0;
            refracted_green = 0.0;
            refracted_blue  = 0.0;
            rec_reflected_red   = 0.0;
            rec_reflected_green = 0.0;
            rec_reflected_blue  = 0.0;
            rec_refracted_red   = 0.0;
            rec_refracted_green = 0.0;
            rec_refracted_blue  = 0.0;


            is_reflected = false;
            is_refracted = false;
            is_rec_reflected = false;
            is_rec_refracted = false;
            is_rec_calculating = false;

            in_shadow    = false;

            // LOCAL INTENSITY DUE TO LIGHT SOURCES

            // calculate intensities contributed by each light source
            for( int j = 0; j < no_of_light_source; j++ )
            {

                // diffuse and specular contributed by this particular light source
                float local_diffuse = 0.0, local_specular = 0.0;

                // light source position
                vertex *light_source_position = new vertex( lights.at( j ).position.x, lights.at( j ).position.y, lights.at( j ).position.z);

                 // light direction unit vector is negated to match sir's conventions
                direction_unit_vector *light_dir = new direction_unit_vector( *light_source_position, *surface_position);
                -(*light_dir);

                // SHADOW, if in shadow, check for the next light source

                if( is_shadow(*surface_position, *light_source_position, SPHERE, i) ) continue;

                // DIFFUSE

                float diffuse_dot_prod = dot_product( *light_dir, *surface_normal);

                // if (tx,ty,tz) is not visible by the light source, no need to proceed for this light source
                if( diffuse_dot_prod < 0 ) continue;

                local_diffuse = spheres.at( i ).kd * lights.at( j ).intensity * diffuse_dot_prod;

                // SPECULAR

                // compute Half Way Vector H = (L+V)/|L+V|, NOTE half_way_vector is not a pointer
                direction_unit_vector half_way_vector = *light_dir + *view_dir;
                half_way_vector.unit_vectorize();

                float specular_dot_prod = pow ( dot_product( *surface_normal, half_way_vector), spheres.at( i ).n);
                local_specular += spheres.at( i ).ks * lights.at( j ).intensity * specular_dot_prod;

                //ATTENUATE THE INTENSITIES

                local_diffuse  /= constant_atten + (linear_atten)*distance(*light_source_position, *surface_position) + (quadratic_atten)*distance(*light_source_position, *surface_position)*distance(*light_source_position, *surface_position);
                local_specular /= constant_atten + (linear_atten)*distance(*light_source_position, *surface_position) + (quadratic_atten)*distance(*light_source_position, *surface_position)*distance(*light_source_position, *surface_position);

                diffuse += local_diffuse;
                specular += local_specular;
            }


            // REFLECTION OF OTHER OBJECTS AT (tx, ty, tz)

            // this is the reflection direction
            direction_unit_vector reflect_dir;
            direction_unit_vector temp = *surface_normal * (2 * dot_product( *view_dir, *surface_normal));
            reflect_dir = temp - *view_dir;
            reflect_dir.unit_vectorize();

            // t is at infinity i.e. the object to be reflected is at infinity
            float t = 1000;
            vertex *start = surface_position;
            vertex *direction = new vertex( reflect_dir.x, reflect_dir.y, reflect_dir.z);

            // this is the direction of ray from (tx,ty,tz) in the direction of reflection
            vertex ray = *start + *direction * t;

            // checking for each sphere, method from sir's slides
            for( int j = 0; j < no_of_spheres; j++ )
            {
                // if negligible reflection coefficient, skip reflection calculations
                if( spheres.at(i).refl_cof <= 0.01 ) break;

                if( i == j ) continue;

                float local_t;

                float A = pow( reflect_dir.x, 2) + pow( reflect_dir.y, 2) + pow( reflect_dir.z, 2);
                float B = 2 * ( reflect_dir.x*(start->x - spheres.at(j).center.x) + reflect_dir.y*(start->y - spheres.at(j).center.y) + reflect_dir.z*(start->z - spheres.at(j).center.z) );
                float C = pow( start->x - spheres.at(j).center.x, 2) + pow( start->y - spheres.at(j).center.y, 2) + pow( start->z - spheres.at(j).center.z, 2) - pow( spheres.at(j).radius, 2);
                float D = B*B - 4*A*C;

                if( D < 0 ) continue;

                float t0 = (-B - sqrt(D))/2*A;
                float t1 = (-B + sqrt(D))/2*A;
                if( t0 < t1 ) local_t = t0; else local_t = t1;

                if( local_t < 0 ) continue;

                // farther from some other object, will not be reflected
                if( local_t > t ) continue;

                t = local_t;

                // ray is the point which is reflected at (tx, ty, tz) and intensity at (ray.x, ray.y, ray.z) will be computed now
                ray = *start + *direction * t;

                vertex *cen = new vertex( spheres.at(j).center.x, spheres.at(j).center.y, spheres.at(j).center.z);
                direction_unit_vector *surface_normal = new direction_unit_vector(*cen, ray);

                // this is from surface to eye to match sir's convention
                direction_unit_vector *view_dir = new direction_unit_vector( ray, *eye_position);

                reflected_red   = calculate_light_intensity( RED,   SPHERE, j, ray, *surface_normal, *view_dir, reflect_dir );
                reflected_green = calculate_light_intensity( GREEN, SPHERE, j, ray, *surface_normal, *view_dir, reflect_dir );
                reflected_blue  = calculate_light_intensity( BLUE,  SPHERE, j, ray, *surface_normal, *view_dir, reflect_dir );

                is_reflected = true;

            }

            // checking for each triangle, method from sir's slides
            for( int j = 0; j < no_of_triangles; j++ )
            {
                // if negligible reflection coefficient, skip reflection calculations
                if( spheres.at(i).refl_cof <= 0.01 ) break;


                float local_t;

                direction_unit_vector pn = triangles.at(j).surface_normal;
                float eqa = triangles.at(j).eqa;
                float eqb = triangles.at(j).eqb;
                float eqc = triangles.at(j).eqc;
                float eqd = triangles.at(j).eqd;

                if( eqa*reflect_dir.x + eqb*reflect_dir.y + eqc*reflect_dir.z == 0 ) continue;

                local_t =  eqa*start->x + eqb*start->y + eqc*start->z + eqd;
                local_t /= eqa*reflect_dir.x + eqb*reflect_dir.y + eqc*reflect_dir.z;
                local_t = -local_t;

                if( local_t < 0 ) continue;

                // farther from some other object, will not be reflected
                if( local_t > t ) continue;

                direction_unit_vector surface_normal = pn;

                // ray is the point which is reflected at (tx, ty, tz) and intensity at (ray.x, ray.y, ray.z) will be computed now
                ray = *start + *direction * local_t;

                // containment test
                vertex v1 = triangles.at(j).v1, v2 = triangles.at(j).v2, v3 = triangles.at(j).v3;

                vertex cp1 = cross_product_vertex( v1-v2,  v3-v2);
                vertex cp2 = cross_product_vertex( v1-v2, ray-v2);

                if( cp1*cp2 < 0 ) continue;

                vertex cp3 = cross_product_vertex( v2-v3,  v1-v3);
                vertex cp4 = cross_product_vertex( v2-v3, ray-v3);

                if( cp3*cp4 < 0 ) continue;

                vertex cp5 = cross_product_vertex( v3-v1,  v2-v1);
                vertex cp6 = cross_product_vertex( v3-v1, ray-v1);

                if( cp5*cp6 < 0 ) continue;

                t = local_t;

                // this is from surface to eye to match sir's convention
                direction_unit_vector *view_dir = new direction_unit_vector( ray, *eye_position);

                reflected_red   = calculate_light_intensity( RED,   TRIANGLE, j, ray, surface_normal, *view_dir, reflect_dir );
                reflected_green = calculate_light_intensity( GREEN, TRIANGLE, j, ray, surface_normal, *view_dir, reflect_dir );
                reflected_blue  = calculate_light_intensity( BLUE,  TRIANGLE, j, ray, surface_normal, *view_dir, reflect_dir );

                is_reflected = true;

            }

            // REFLECTION PART OVER

            // REFRACTION OF OTHER OBJECTS AT (tx, ty, tz)
            // NOTE : THIS IS THE COMPLEX PART BECAUSE IT IS REFRACTION THROUGH SPHERE

            // to find direction of refracted ray, method from sir'slides
            float cos_theta_i    = dot_product( *view_dir, *surface_normal);
            float sin_theta_i    = 1 - pow(cos_theta_i, 2);
            float refr_ind_ratio = 1/(spheres.at(i).refr_ind);
            float sin_theta_t    = refr_ind_ratio * sin_theta_i;
            float cos_theta_t    = 1 - pow(sin_theta_t, 2);

            direction_unit_vector temp_I = *view_dir;
            -temp_I;

            vertex *I = new vertex( temp_I.x         , temp_I.y         , temp_I.z          );
            vertex *N = new vertex( surface_normal->x, surface_normal->y, surface_normal->z );

            // this is the refraction direction from surface to in the direction of refracted ray
            *direction = ( *I + *N * cos_theta_i ) * refr_ind_ratio - *N * cos_theta_t;

            direction_unit_vector *refract_dir = new direction_unit_vector( new vertex(0.0,0.0,0.0), direction);

            // t is at infinity i.e. the object to be refracted is at infinity
            t = 1000;
            start = surface_position;
            direction = new vertex( refract_dir->x, refract_dir->y, refract_dir->z);

            // this is the direction of ray from (tx,ty,tz) in the direction of refraction i.e. going inside the sphere
            ray = *start + *direction * t;

            // now, we have to find the position from where this ray will come out of the sphere
            float local_t;

            float A = pow( refract_dir->x, 2) + pow( refract_dir->y, 2) + pow( refract_dir->z, 2);
            float B = 2 * ( refract_dir->x*(start->x - spheres.at(i).center.x) + refract_dir->y*(start->y - spheres.at(i).center.y) + refract_dir->z*(start->z - spheres.at(i).center.z) );
            float C = pow( start->x - spheres.at(i).center.x, 2) + pow( start->y - spheres.at(i).center.y, 2) + pow( start->z - spheres.at(i).center.z, 2) - pow( spheres.at(i).radius, 2);
            float D = B*B - 4*A*C;

            if( D < 0 ) continue;

            float t0 = (-B - sqrt(D))/2*A;
            float t1 = (-B + sqrt(D))/2*A;
            if( t0 < t1 ) local_t = t1; else local_t = t0;

            if( local_t < 0 ) continue;

            t = local_t;

            // ray is the point which is refracted at (tx, ty, tz) and intensity at (ray.x, ray.y, ray.z) will be computed now
            ray = *start + *direction * t;

            vertex *cen = new vertex( spheres.at(i).center.x, spheres.at(i).center.y, spheres.at(i).center.z);
            surface_normal = new direction_unit_vector(ray, *cen);

            // this is from surface to eye to match sir's convention
            view_dir = new direction_unit_vector( ray, *surface_position);

            // to find direction of refracted ray, method from sir'slides
            cos_theta_i    = dot_product( *view_dir, *surface_normal);
            sin_theta_i    = 1 - pow(cos_theta_i, 2);
            refr_ind_ratio = spheres.at(i).refr_ind;
            sin_theta_t    = refr_ind_ratio * sin_theta_i;
            cos_theta_t    = 1 - pow(sin_theta_t, 2);

            temp_I = *view_dir;
            -temp_I;

            I = new vertex( temp_I.x         , temp_I.y         , temp_I.z          );
            N = new vertex( surface_normal->x, surface_normal->y, surface_normal->z );

            // this is the refraction direction from surface to in the direction of refracted ray
            *direction = ( *I + *N * cos_theta_i ) * refr_ind_ratio - *N * cos_theta_t;

            refract_dir = new direction_unit_vector( new vertex(0.0,0.0,0.0), direction);

            // t is at infinity i.e. the object to be refracted is at infinity
            t = 1000;
            start = surface_position;
            direction = new vertex( refract_dir->x, refract_dir->y, refract_dir->z);

            // this is the direction of ray from (tx,ty,tz) in the direction of refraction i.e. coming outside of the sphere
            ray = *start + *direction * t;

            // checking for each sphere, method from sir's slides
            for( int j = 0; j < no_of_spheres; j++ )
            {
                // if negligible refraction coefficient, skip refraction calculations
                if( spheres.at(i).refr_cof <= 0.01 ) break;


                if( i == j ) continue;

                float local_t;

                float A = pow( refract_dir->x, 2) + pow( refract_dir->y, 2) + pow( refract_dir->z, 2);
                float B = 2 * ( refract_dir->x*(start->x - spheres.at(j).center.x) + refract_dir->y*(start->y - spheres.at(j).center.y) + refract_dir->z*(start->z - spheres.at(j).center.z) );
                float C = pow( start->x - spheres.at(j).center.x, 2) + pow( start->y - spheres.at(j).center.y, 2) + pow( start->z - spheres.at(j).center.z, 2) - pow( spheres.at(j).radius, 2);
                float D = B*B - 4*A*C;

                if( D < 0 ) continue;

                float t0 = (-B - sqrt(D))/2*A;
                float t1 = (-B + sqrt(D))/2*A;
                if( t0 < t1 ) local_t = t0; else local_t = t1;

                if( local_t < 0 ) continue;

                // farther from some other object, will not be refracted
                if( local_t > t ) continue;

                t = local_t;

                // ray is the point which is refracted at (tx, ty, tz) and intensity at (ray.x, ray.y, ray.z) will be computed now
                ray = *start + *direction * t;

                vertex *cen = new vertex( spheres.at(j).center.x, spheres.at(j).center.y, spheres.at(j).center.z);
                direction_unit_vector *surface_normal = new direction_unit_vector(*cen, ray);

                // this is from surface to eye to match sir's convention
                direction_unit_vector *view_dir = new direction_unit_vector( ray, *eye_position);

                refracted_red   = calculate_light_intensity( RED,   SPHERE, j, ray, *surface_normal, *view_dir, *refract_dir );
                refracted_green = calculate_light_intensity( GREEN, SPHERE, j, ray, *surface_normal, *view_dir, *refract_dir );
                refracted_blue  = calculate_light_intensity( BLUE,  SPHERE, j, ray, *surface_normal, *view_dir, *refract_dir );

                is_refracted = true;

            }

            // checking for each triangle, method from sir's slides
            for( int j = 0; j < no_of_triangles; j++ )
            {
                // if negligible refraction coefficient, skip refraction calculations
                if( spheres.at(i).refr_cof <= 0.01 ) break;

                float local_t;

                direction_unit_vector pn = triangles.at(j).surface_normal;
                float eqa = triangles.at(j).eqa;
                float eqb = triangles.at(j).eqb;
                float eqc = triangles.at(j).eqc;
                float eqd = triangles.at(j).eqd;

                if( eqa*refract_dir->x + eqb*refract_dir->y + eqc*refract_dir->z == 0 ) continue;

                local_t =  eqa*start->x + eqb*start->y + eqc*start->z + eqd;
                local_t /= eqa*refract_dir->x + eqb*refract_dir->y + eqc*refract_dir->z;
                local_t = -local_t;

                if( local_t < 0 ) continue;

                // farther from some other object, will not be refracted
                if( local_t > t ) continue;

                direction_unit_vector surface_normal = pn;

                // ray is the point which is refracted at (tx, ty, tz) and intensity at (ray.x, ray.y, ray.z) will be computed now
                ray = *start + *direction * local_t;

                // containment test
                vertex v1 = triangles.at(j).v1, v2 = triangles.at(j).v2, v3 = triangles.at(j).v3;

                vertex cp1 = cross_product_vertex( v1-v2,  v3-v2);
                vertex cp2 = cross_product_vertex( v1-v2, ray-v2);

                if( cp1*cp2 < 0 ) continue;

                vertex cp3 = cross_product_vertex( v2-v3,  v1-v3);
                vertex cp4 = cross_product_vertex( v2-v3, ray-v3);

                if( cp3*cp4 < 0 ) continue;

                vertex cp5 = cross_product_vertex( v3-v1,  v2-v1);
                vertex cp6 = cross_product_vertex( v3-v1, ray-v1);

                if( cp5*cp6 < 0 ) continue;

                t = local_t;

                // this is from surface to eye to match sir's convention
                direction_unit_vector *view_dir = new direction_unit_vector( ray, *eye_position);

                refracted_red   = calculate_light_intensity( RED,   TRIANGLE, j, ray, surface_normal, *view_dir, *refract_dir );
                refracted_green = calculate_light_intensity( GREEN, TRIANGLE, j, ray, surface_normal, *view_dir, *refract_dir );
                refracted_blue  = calculate_light_intensity( BLUE,  TRIANGLE, j, ray, surface_normal, *view_dir, *refract_dir );

                is_refracted = true;

            }

            // REFRACTION PART OVER

            if( is_reflected )
            {
                if( is_refracted )
                {
                    if( spheres.at(i).refl_cof + spheres.at(i).refr_cof > 1 )
                    {
                        spheres.at(i).refl_cof /= spheres.at(i).refl_cof + spheres.at(i).refr_cof;
                        spheres.at(i).refr_cof /= spheres.at(i).refl_cof + spheres.at(i).refr_cof;
                    }

                    glColor3f((ambient * spheres.at( i ).ka + diffuse) * spheres.at( i ).cr * (1 - spheres.at(i).refl_cof - spheres.at(i).refr_cof ) + specular + reflected_red   * ( spheres.at(i).refl_cof ) + refracted_red   * ( spheres.at(i).refr_cof ),
                              (ambient * spheres.at( i ).ka + diffuse) * spheres.at( i ).cg * (1 - spheres.at(i).refl_cof - spheres.at(i).refr_cof ) + specular + reflected_green * ( spheres.at(i).refl_cof ) + refracted_green * ( spheres.at(i).refr_cof ),
                              (ambient * spheres.at( i ).ka + diffuse) * spheres.at( i ).cb * (1 - spheres.at(i).refl_cof - spheres.at(i).refr_cof ) + specular + reflected_blue  * ( spheres.at(i).refl_cof ) + refracted_blue  * ( spheres.at(i).refr_cof ));
                }
                else
                    glColor3f((ambient * spheres.at( i ).ka + diffuse) * spheres.at( i ).cr * (1 - spheres.at(i).refl_cof ) + specular + reflected_red   * ( spheres.at(i).refl_cof ),
                              (ambient * spheres.at( i ).ka + diffuse) * spheres.at( i ).cg * (1 - spheres.at(i).refl_cof ) + specular + reflected_green * ( spheres.at(i).refl_cof ),
                              (ambient * spheres.at( i ).ka + diffuse) * spheres.at( i ).cb * (1 - spheres.at(i).refl_cof ) + specular + reflected_blue  * ( spheres.at(i).refl_cof ));
            }
            else if( is_refracted )
                glColor3f((ambient * spheres.at( i ).ka + diffuse) * spheres.at( i ).cr * (1 - spheres.at(i).refr_cof ) + specular + refracted_red   * ( spheres.at(i).refr_cof ),
                          (ambient * spheres.at( i ).ka + diffuse) * spheres.at( i ).cg * (1 - spheres.at(i).refr_cof ) + specular + refracted_green   * ( spheres.at(i).refr_cof ),
                          (ambient * spheres.at( i ).ka + diffuse) * spheres.at( i ).cb * (1 - spheres.at(i).refr_cof ) + specular + refracted_blue   * ( spheres.at(i).refr_cof ));
            else
            {
                glPointSize( 3.0 );
                glColor3f((ambient * spheres.at( i ).ka + diffuse) * spheres.at( i ).cr + specular,
                          (ambient * spheres.at( i ).ka + diffuse) * spheres.at( i ).cg + specular,
                          (ambient * spheres.at( i ).ka + diffuse) * spheres.at( i ).cb + specular);
            }
            tx = tx * sc_x + ty * sh_xy;
            ty *= sc_y;
            tz *= sc_z;

            glVertex3f( tx, ty, tz);
        }
    }

    glEnd();
    glPopMatrix();

  }

  // DRAW TRIANGLES
  for( int i = 0; i < no_of_triangles; i++ )
  {
    std::cout << "T" << i << "\n";

    float tx , ty, tz;

    glPushMatrix();

    // end points of triangle
    vertex v1 = triangles.at( i ).v1;
    vertex v2 = triangles.at( i ).v2;
    vertex v3 = triangles.at( i ).v3;

    direction_unit_vector surface_normal = triangles.at( i ).surface_normal;

    glBegin(GL_POINTS);

    int count = 0;
    for( float boundary_t = 0; boundary_t <= 1; boundary_t += boundary_t_inc )
    {
        float x_begin = boundary_t * (v3.x - v1.x) + v1.x ;
        float y_begin = boundary_t * (v3.y - v1.y) + v1.y ;
        float z_begin = boundary_t * (v3.z - v1.z) + v1.z ;

        float x_end = boundary_t * (v3.x - v2.x) + v2.x ;
        float y_end = boundary_t * (v3.y - v2.y) + v2.y ;
        float z_end = boundary_t * (v3.z - v2.z) + v2.z ;

        // draw line from v_begin to v_end
        for( float inner_t = 0; inner_t <= 1; inner_t += inner_t_inc )
        {
            glPointSize( point_size );

            // triangle vertex at this point is tx,ty,tz
            tx = inner_t * (x_end - x_begin) + x_begin ;
            ty = inner_t * (y_end - y_begin) + y_begin ;
            tz = inner_t * (z_end - z_begin) + z_begin ;

            // this is the vertex to be illuminated
            vertex *surface_position = new vertex( tx, ty, tz);

            // this point is obstructed by some other object, hence not visible
            if( !is_visible( *surface_position, TRIANGLE, i) ) continue;

            // this is eye position and direction from surface to eye to match sir's convention
            vertex *eye_position = new vertex( eyex, eyey, eyez);
            direction_unit_vector *view_dir = new direction_unit_vector( *surface_position, *eye_position);

            // initially, no vertex is illuminated either by diffuse or specular intensity at ( tx, ty, tz)
            // initially, no vertex is reflecting any other object
            diffuse         = 0.0;
            specular        = 0.0;
            reflected_red   = 0.0;
            reflected_green = 0.0;
            reflected_blue  = 0.0;
            refracted_red   = 0.0;
            refracted_green = 0.0;
            refracted_blue  = 0.0;
            rec_reflected_red   = 0.0;
            rec_reflected_green = 0.0;
            rec_reflected_blue  = 0.0;
            rec_refracted_red   = 0.0;
            rec_refracted_green = 0.0;
            rec_refracted_blue  = 0.0;


            is_reflected = false;
            is_refracted = false;
            is_rec_reflected = false;
            is_rec_refracted = false;
            is_rec_calculating = false;

            in_shadow    = false;

            // calculate intensities contributed by each light source
            for( int j = 0; j < no_of_light_source; j++ )
            {

                // diffuse and specular contributed by this particular light source
                float local_diffuse = 0.0, local_specular = 0.0;

               // light source position
                vertex *light_source_position = new vertex( lights.at( j ).position.x, lights.at( j ).position.y, lights.at( j ).position.z);

                 // light direction unit vector is negated to match sir's conventions
                direction_unit_vector *light_dir = new direction_unit_vector( *light_source_position, *surface_position);
                -(*light_dir);

                // SHADOW, if in shadow, check for the next light source

                if( is_shadow(*surface_position, *light_source_position, TRIANGLE, i) ) continue;

                // DIFFUSE

                float diffuse_dot_prod;

                if( dot_product( *view_dir, surface_normal) > 0 && dot_product( *light_dir, surface_normal) > 0 )
                {
                    diffuse_dot_prod = dot_product( *light_dir, surface_normal);
                    local_diffuse = triangles.at( i ).kd * lights.at( j ).intensity * diffuse_dot_prod;
                }
                else if( dot_product( *view_dir, surface_normal) < 0 && dot_product( *light_dir, surface_normal) < 0 )
                {
                    diffuse_dot_prod = -dot_product( *light_dir, surface_normal);
                    local_diffuse = triangles.at( i ).kd * lights.at( j ).intensity * diffuse_dot_prod;
                }
                else
                    continue;

                // SPECULAR

                // compute Half Way Vector H = (L+V)/|L+V|, NOTE half_way_vector is not a pointer
                direction_unit_vector half_way_vector = *light_dir + *view_dir;
                half_way_vector.unit_vectorize();

                float specular_dot_prod = pow ( dot_product( surface_normal, half_way_vector), triangles.at( i ).n);
                local_specular = triangles.at( i ).ks * lights.at( j ).intensity * specular_dot_prod;

                //ATTENUATE THE INTENSITIES

                local_diffuse  /= constant_atten + (linear_atten)*distance(*light_source_position, *surface_position) + (quadratic_atten)*distance(*light_source_position, *surface_position)*distance(*light_source_position, *surface_position);
                local_specular /= constant_atten + (linear_atten)*distance(*light_source_position, *surface_position) + (quadratic_atten)*distance(*light_source_position, *surface_position)*distance(*light_source_position, *surface_position);

                diffuse += local_diffuse;
                specular += local_specular;
            }


            // REFLECTION OF OTHER OBJECTS AT (tx, ty, tz)

            // this is the reflection direction
            direction_unit_vector reflect_dir;
            direction_unit_vector temp = surface_normal * (2 * dot_product( *view_dir, surface_normal));
            reflect_dir = temp - *view_dir;

            // t is at infinity i.e. the object to be reflected is at infinity
            float t = 1000;
            vertex *start = surface_position;
            vertex *direction = new vertex( reflect_dir.x, reflect_dir.y, reflect_dir.z);

            // this is the direction of ray from (tx,ty,tz) in the direction of reflection
            vertex ray = *start + *direction * t;

            // checking for each sphere, method from sir's slides
            for( int j = 0; j < no_of_spheres; j++ )
            {
                // if negligible reflection coefficient, skip reflection calculations
                if( triangles.at(i).refl_cof <= 0.01 ) break;

                float local_t;

                float A = pow( reflect_dir.x, 2) + pow( reflect_dir.y, 2) + pow( reflect_dir.z, 2);
                float B = 2 * ( reflect_dir.x*(start->x - spheres.at(j).center.x) + reflect_dir.y*(start->y - spheres.at(j).center.y) + reflect_dir.z*(start->z - spheres.at(j).center.z) );
                float C = pow( start->x - spheres.at(j).center.x, 2) + pow( start->y - spheres.at(j).center.y, 2) + pow( start->z - spheres.at(j).center.z, 2) - pow( spheres.at(j).radius, 2);
                float D = B*B - 4*A*C;

                if( D < 0 ) continue;

                float t0 = (-B - sqrt(D))/2*A;
                float t1 = (-B + sqrt(D))/2*A;
                if( t0 < t1 ) local_t = t0; else local_t = t1;

                if( local_t < 0 ) continue;

                // farther from some other object, will not be reflected
                if( local_t > t ) continue;

                t = local_t;

                // ray is the point which is reflected at (tx, ty, tz) and intensity at (ray.x, ray.y, ray.z) will be computed now
                ray = *start + *direction * t;

                vertex *cen = new vertex( spheres.at(j).center.x, spheres.at(j).center.y, spheres.at(j).center.z);
                direction_unit_vector *surface_normal = new direction_unit_vector(*cen, ray);

                // this is from surface to eye to match sir's convention
                direction_unit_vector *view_dir = new direction_unit_vector( ray, *eye_position);

                reflected_red   = calculate_light_intensity( RED,   SPHERE, j, ray, *surface_normal, *view_dir, reflect_dir );
                reflected_green = calculate_light_intensity( GREEN, SPHERE, j, ray, *surface_normal, *view_dir, reflect_dir  );
                reflected_blue  = calculate_light_intensity( BLUE,  SPHERE, j, ray, *surface_normal, *view_dir, reflect_dir  );

                is_reflected = true;

            }

            // checking for each triangle, method from sir's slides
            for( int j = 0; j < no_of_triangles; j++ )
            {
                // if negligible reflection coefficient, skip reflection calculations
                if( triangles.at(i).refl_cof <= 0.01 ) break;

                if( i == j ) continue;

                float local_t;

                direction_unit_vector pn = triangles.at(j).surface_normal;
                float eqa = triangles.at(j).eqa;
                float eqb = triangles.at(j).eqb;
                float eqc = triangles.at(j).eqc;
                float eqd = triangles.at(j).eqd;

                if( eqa*reflect_dir.x + eqb*reflect_dir.y + eqc*reflect_dir.z == 0 ) continue;

                local_t =  eqa*start->x + eqb*start->y + eqc*start->z + eqd;
                local_t /= eqa*reflect_dir.x + eqb*reflect_dir.y + eqc*reflect_dir.z;
                local_t = -local_t;

                if( local_t < 0 ) continue;

                // farther from some other object, will not be reflected
                if( local_t > t ) continue;

                direction_unit_vector surface_normal = pn;

                // ray is the point which is reflected at (tx, ty, tz) and intensity at (ray.x, ray.y, ray.z) will be computed now
                ray = *start + *direction * local_t;

                // containment test
                vertex v1 = triangles.at(j).v1, v2 = triangles.at(j).v2, v3 = triangles.at(j).v3;

                vertex cp1 = cross_product_vertex( v1-v2,  v3-v2);
                vertex cp2 = cross_product_vertex( v1-v2, ray-v2);

                if( cp1*cp2 < 0 ) continue;

                vertex cp3 = cross_product_vertex( v2-v3,  v1-v3);
                vertex cp4 = cross_product_vertex( v2-v3, ray-v3);

                if( cp3*cp4 < 0 ) continue;

                vertex cp5 = cross_product_vertex( v3-v1,  v2-v1);
                vertex cp6 = cross_product_vertex( v3-v1, ray-v1);

                if( cp5*cp6 < 0 ) continue;

                t = local_t;

                // this is from surface to eye to match sir's convention
                direction_unit_vector *view_dir = new direction_unit_vector( ray, *eye_position);

                reflected_red   = calculate_light_intensity( RED,   TRIANGLE, j, ray, surface_normal, *view_dir, reflect_dir );
                reflected_green = calculate_light_intensity( GREEN, TRIANGLE, j, ray, surface_normal, *view_dir, reflect_dir );
                reflected_blue  = calculate_light_intensity( BLUE,  TRIANGLE, j, ray, surface_normal, *view_dir, reflect_dir );

                is_reflected = true;

            }

            // REFLECTION PART OVER

            // REFRACTION OF OTHER OBJECTS AT (tx, ty, tz)


            // this is the refraction direction from surface to in the direction of refracted ray
            direction_unit_vector refract_dir = *view_dir;
            -refract_dir;

            // t is at infinity i.e. the object to be refracted is at infinity
            t = 1000;
            start = surface_position;
            direction = new vertex( refract_dir.x, refract_dir.y, refract_dir.z);

            // this is the direction of ray from (tx,ty,tz) in the direction of refraction
            ray = *start + *direction * t;

            // checking for each sphere, method from sir's slides
            for( int j = 0; j < no_of_spheres; j++ )
            {
                // if negligible refraction coefficient, skip refraction calculations
                if( triangles.at(i).refr_cof <= 0.01 ) break;

                float local_t;

                float A = pow( refract_dir.x, 2) + pow( refract_dir.y, 2) + pow( refract_dir.z, 2);
                float B = 2 * ( refract_dir.x*(start->x - spheres.at(j).center.x) + refract_dir.y*(start->y - spheres.at(j).center.y) + refract_dir.z*(start->z - spheres.at(j).center.z) );
                float C = pow( start->x - spheres.at(j).center.x, 2) + pow( start->y - spheres.at(j).center.y, 2) + pow( start->z - spheres.at(j).center.z, 2) - pow( spheres.at(j).radius, 2);
                float D = B*B - 4*A*C;

                if( D < 0 ) continue;

                float t0 = (-B - sqrt(D))/2*A;
                float t1 = (-B + sqrt(D))/2*A;
                if( t0 < t1 ) local_t = t0; else local_t = t1;

                if( local_t < 0 ) continue;

                // farther from some other object, will not be refracted
                if( local_t > t ) continue;

                t = local_t;

                // ray is the point which is refracted at (tx, ty, tz) and intensity at (ray.x, ray.y, ray.z) will be computed now
                ray = *start + *direction * t;

                vertex *cen = new vertex( spheres.at(j).center.x, spheres.at(j).center.y, spheres.at(j).center.z);
                direction_unit_vector *surface_normal = new direction_unit_vector(*cen, ray);

                // this is from surface to eye to match sir's convention
                direction_unit_vector *view_dir = new direction_unit_vector( ray, *eye_position);

                refracted_red   = calculate_light_intensity( RED,   SPHERE, j, ray, *surface_normal, *view_dir, refract_dir );
                refracted_green = calculate_light_intensity( GREEN, SPHERE, j, ray, *surface_normal, *view_dir, refract_dir );
                refracted_blue  = calculate_light_intensity( BLUE,  SPHERE, j, ray, *surface_normal, *view_dir, refract_dir );

                is_refracted = true;
            }

            // checking for each triangle, method from sir's slides
            for( int j = 0; j < no_of_triangles; j++ )
            {
                // if negligible refraction coefficient, skip refraction calculations
                if( triangles.at(i).refr_cof <= 0.01 ) break;

                if( i == j ) continue;

                float local_t;

                direction_unit_vector pn = triangles.at(j).surface_normal;
                float eqa = triangles.at(j).eqa;
                float eqb = triangles.at(j).eqb;
                float eqc = triangles.at(j).eqc;
                float eqd = triangles.at(j).eqd;

                if( eqa*refract_dir.x + eqb*refract_dir.y + eqc*refract_dir.z == 0 ) continue;

                local_t =  eqa*start->x + eqb*start->y + eqc*start->z + eqd;
                local_t /= eqa*refract_dir.x + eqb*refract_dir.y + eqc*refract_dir.z;
                local_t = -local_t;

                if( local_t < 0 ) continue;

                // farther from some other object, will not be refracted
                if( local_t > t ) continue;

                direction_unit_vector surface_normal = pn;

                // ray is the point which is refracted at (tx, ty, tz) and intensity at (ray.x, ray.y, ray.z) will be computed now
                ray = *start + *direction * local_t;

                // containment test
                vertex v1 = triangles.at(j).v1, v2 = triangles.at(j).v2, v3 = triangles.at(j).v3;

                vertex cp1 = cross_product_vertex( v1-v2,  v3-v2);
                vertex cp2 = cross_product_vertex( v1-v2, ray-v2);

                if( cp1*cp2 < 0 ) continue;

                vertex cp3 = cross_product_vertex( v2-v3,  v1-v3);
                vertex cp4 = cross_product_vertex( v2-v3, ray-v3);

                if( cp3*cp4 < 0 ) continue;

                vertex cp5 = cross_product_vertex( v3-v1,  v2-v1);
                vertex cp6 = cross_product_vertex( v3-v1, ray-v1);

                if( cp5*cp6 < 0 ) continue;

                t = local_t;

                // this is from surface to eye to match sir's convention
                direction_unit_vector *view_dir = new direction_unit_vector( ray, *eye_position);

                refracted_red   = calculate_light_intensity( RED,   TRIANGLE, j, ray, surface_normal, *view_dir, refract_dir );
                refracted_green = calculate_light_intensity( GREEN, TRIANGLE, j, ray, surface_normal, *view_dir, refract_dir );
                refracted_blue  = calculate_light_intensity( BLUE,  TRIANGLE, j, ray, surface_normal, *view_dir, refract_dir );

                is_refracted = true;

            }

            // REFRACTION PART OVER

            if( is_reflected )
            {
                if( is_refracted )
                {
                    if( triangles.at(i).refl_cof + triangles.at(i).refr_cof > 1 )
                    {
                        triangles.at(i).refl_cof /= triangles.at(i).refl_cof + triangles.at(i).refr_cof;
                        triangles.at(i).refr_cof /= triangles.at(i).refl_cof + triangles.at(i).refr_cof;
                    }

                    glColor3f((ambient * triangles.at( i ).ka + diffuse) * triangles.at( i ).cr * (1 - triangles.at(i).refl_cof - triangles.at(i).refr_cof ) + specular + reflected_red   * ( triangles.at(i).refl_cof ) + refracted_red   * ( triangles.at(i).refr_cof ),
                              (ambient * triangles.at( i ).ka + diffuse) * triangles.at( i ).cg * (1 - triangles.at(i).refl_cof - triangles.at(i).refr_cof ) + specular + reflected_green * ( triangles.at(i).refl_cof ) + refracted_green * ( triangles.at(i).refr_cof ),
                              (ambient * triangles.at( i ).ka + diffuse) * triangles.at( i ).cb * (1 - triangles.at(i).refl_cof - triangles.at(i).refr_cof ) + specular + reflected_blue  * ( triangles.at(i).refl_cof ) + refracted_blue  * ( triangles.at(i).refr_cof ));
                }
                else
                    glColor3f((ambient * triangles.at( i ).ka + diffuse) * triangles.at( i ).cr * (1 - triangles.at(i).refl_cof ) + specular + reflected_red   * ( triangles.at(i).refl_cof ),
                              (ambient * triangles.at( i ).ka + diffuse) * triangles.at( i ).cg * (1 - triangles.at(i).refl_cof ) + specular + reflected_green * ( triangles.at(i).refl_cof ),
                              (ambient * triangles.at( i ).ka + diffuse) * triangles.at( i ).cb * (1 - triangles.at(i).refl_cof ) + specular + reflected_blue  * ( triangles.at(i).refl_cof ));
            }
            else if( is_refracted )
                glColor3f((ambient * triangles.at( i ).ka + diffuse) * triangles.at( i ).cr * (1 - triangles.at(i).refr_cof ) + specular + refracted_red   * ( triangles.at(i).refr_cof ),
                          (ambient * triangles.at( i ).ka + diffuse) * triangles.at( i ).cg * (1 - triangles.at(i).refr_cof ) + specular + refracted_green   * ( triangles.at(i).refr_cof ),
                          (ambient * triangles.at( i ).ka + diffuse) * triangles.at( i ).cb * (1 - triangles.at(i).refr_cof ) + specular + refracted_blue   * ( triangles.at(i).refr_cof ));
            else
            {
                glColor3f((ambient * triangles.at( i ).ka + diffuse) * triangles.at( i ).cr + specular,
                          (ambient * triangles.at( i ).ka + diffuse) * triangles.at( i ).cg + specular,
                          (ambient * triangles.at( i ).ka + diffuse) * triangles.at( i ).cb + specular);
            }
            glVertex3f( tx, ty, tz);
        }
    }

    glEnd();
    glPopMatrix();

  }

  glutSwapBuffers();
}

void reshape (int w, int h)
{
  glViewport (0, 0, (GLsizei) w, (GLsizei) h);
  glMatrixMode (GL_PROJECTION);
  glLoadIdentity ();
  glFrustum (-width/2, width/2, -height/2, height/2, near_plane, far_plane);
  glMatrixMode(GL_MODELVIEW);
}

void keyboard(unsigned char key, int x, int y)
{
      switch( key )
        {

        case 'a':
          eyex -= 0.5;
          break;
        case 'd':
          eyex += 0.5;
          break;
        case 'w':
          eyez += 0.5;
          break;
        case 's':
          eyez -= 0.5;
          break;

        case 'p':
          for( int i  = 0; i < no_of_light_source; i++ )
            lights.at(i).intensity += 2.0;
          break;
        case 'o':
          for( int i  = 0; i < no_of_light_source; i++ )
            lights.at(i).intensity -= 2.0;
          break;

        case 'l':
          if( eyex == 20.0 )
        behind = true;
          if( eyex == -20.0 )
        behind = false;
          if( eyez == 20.0 )
        right = true;
          if( eyez == -20.0 )
        right = false;

          if( !behind && right )
        {
          eyex += 0.5;
          eyez -= 0.5;
        }
          else if( behind && right )
        {
          eyex -= 0.5;
          eyez -= 0.5;
        }
          else if( behind && !right )
        {
          eyex -= 0.5;
          eyez += 0.5;
        }
          else if( !behind && !right )
        {
          eyex += 0.5;
          eyez += 0.5;
        }
          break;
        case 'j':
          if( eyex == 20.0 )
        behind = false;
          if( eyex == -20.0 )
        behind = true;
          if( eyez == 20.0 )
        right = false;
          if( eyez == -20.0 )
        right = true;

          if( !behind && right )
        {
          eyex -= 0.5;
          eyez += 0.5;
        }
          else if( behind && right )
        {
          eyex += 0.5;
          eyez += 0.5;
        }
          else if( behind && !right )
        {
          eyex += 0.5;
          eyez -= 0.5;
        }
          else if( !behind && !right )
        {
          eyex -= 0.5;
          eyez -= 0.5;
        }
          break;

    }
    glutPostRedisplay();
}

int main(int argc, char** argv)
{
    std::cin >> is_affine;
    if( is_affine == 1 )
    {
        std::cin >> sc_x >> sc_y >> sc_z;
        std::cin >> sh_xy >> sh_yx;
    }

    std::cin >> point_size;

    std::cin >> boundary_t_inc >> inner_t_inc >> theta_inc >> fi_inc;

    std::cin >> no_of_spheres;

    for( int i = 0; i < no_of_spheres; i++ )
    {
        float radius, x, y, z, cr, cg, cb, ka, kd, ks, n, refl_cof, refr_cof, refr_ind;

        std::cin >> radius >> x >> y >> z >> cr >> cg >> cb >> ka >> kd >> ks >> n >> refl_cof >> refr_cof >> refr_ind;

        vertex *center = new vertex( x, y, z);
        sphere *temp = new sphere( radius, *center, cr, cg, cb, ka, kd, ks, n, refl_cof, refr_cof, refr_ind );
        spheres.push_back( *temp );
    }

    std::cin >> no_of_triangles;

    for( int i = 0; i < no_of_triangles; i++ )
    {
        // enter the vertices of the triangle
        float x1, y1, z1, x2, y2, z2, x3, y3, z3;
        std::cin >> x1 >> y1 >> z1 >> x2 >> y2 >> z2 >> x3 >> y3 >> z3;

        // enter other attributes
        float cr, cg, cb, ka, kd, ks, n, refl_cof, refr_cof, refr_ind;

        std::cin >> cr >> cg >> cb >> ka >> kd >> ks >> n >> refl_cof >> refr_cof >> refr_ind;

        vertex *v1 = new vertex(x1,y1,z1), *v2 = new vertex(x2,y2,z2), *v3 = new vertex(x3,y3,z3);

        // find normal to this triangle
        direction_unit_vector *plane_vector_1 = new direction_unit_vector( *v1, *v2);
        direction_unit_vector *plane_vector_2 = new direction_unit_vector( *v1, *v3);
        direction_unit_vector surface_normal = cross_product(*plane_vector_2, *plane_vector_1);
        surface_normal.unit_vectorize();

        float eqa = surface_normal.x, eqb = surface_normal.y, eqc = surface_normal.z, eqd;
        eqd = - ( eqa*x1 + eqb*y1 + eqc*z1 );

        triangle *temp = new triangle( eqa, eqb, eqc, eqd, surface_normal, v1, v2, v3, cr, cg, cb, ka, kd, ks, n, refl_cof, refr_cof, refr_ind );

        triangles.push_back( *temp );
    }

    std::cin >> eyex >> eyey >> eyez;
    std::cin >> lookAt_x >> lookAt_y >> lookAt_z;
    std::cin >> up_x >> up_y >> up_z;
    std::cin >> near_plane >> far_plane;
    std::cin >> width >> height;

    std::cin >> no_of_light_source;

    for( int i = 0; i < no_of_light_source; i++ )
    {
        float x ,y ,z, intensity;
        std::cin >> x >> y >> z >> intensity;
        vertex *vtemp = new vertex( x, y, z );
        light *temp = new light( *vtemp, intensity );
        lights.push_back( *temp );
    }

    std::cin >> constant_atten >> linear_atten >> quadratic_atten;
    std::cin >> ambient;

    glutInit(&argc, argv);
    glutInitDisplayMode (GLUT_DEPTH | GLUT_RGB);
    glutInitWindowSize (1600, 900);
    glutInitWindowPosition (100, 100);
    glutCreateWindow (argv[0]);
    init();
    glutDisplayFunc(display);
    glutReshapeFunc(reshape);
    glutKeyboardFunc(keyboard);
    glutMainLoop();
    return 0;
}
