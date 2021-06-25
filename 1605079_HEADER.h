#include <iostream>
#include<fstream>
#include<stack>
#include<vector>
#include <cmath>
#include<string>
#include<string.h>
#include<sstream>
#include <ctime>
#include <cstdlib>
#include<limits>
#include <iomanip>
#include <windows.h>
#include <GL/glut.h>

#define pi (2*acos(0.0))

using namespace std;

struct Vector3D {
    double x, y, z;

    inline Vector3D operator*(double v) {
        return {x * v, y * v, z * v};
    }

    inline Vector3D operator+(Vector3D p) {
        return {x + p.x, y + p.y, z + p.z};
    }

    inline Vector3D operator-(Vector3D p) {
        return {x - p.x, y - p.y, z - p.z};
    }

};

typedef Vector3D Vector3D;

void printVector3D(Vector3D p) {
    cout << "x = " << p.x << " y = " << p.y << " z = " << p.z << endl;
}

double getVectorMagnitude(Vector3D v) {
    return sqrt(v.x * v.x + v.y * v.y + v.z * v.z);
}

double getDotProduct(Vector3D a, Vector3D b) {
    return a.x * b.x + a.y * b.y + a.z * b.z;
}

Vector3D getUnitVector(Vector3D v) {
    double norm = getVectorMagnitude(v);
    Vector3D normalized = {0, 0, 0};

    normalized.x = v.x / norm;
    normalized.y = v.y / norm;
    normalized.z = v.z / norm;
    return normalized;
}

double convertDegreeToRadian(double theta) {
    return theta * pi / 180;
}

Vector3D getCrossProduct(Vector3D a, Vector3D b) {
    Vector3D result = {0, 0, 0};

    result.x = a.y * b.z - a.z * b.y;
    result.y = a.z * b.x - a.x * b.z;
    result.z = a.x * b.y - b.x * a.y;
    return result;
}


class Ray {
public:
    Vector3D start;
    Vector3D dir;

    // normalize for easier calculations
    // write appropriate constructor
    Ray(Vector3D R_0, Vector3D R_d) {
        start = R_0;
        dir = getUnitVector(R_d);
    }
};

class Object {
public:
    Vector3D reference_point;
    double height, width, length;
    double color[3];
    double coefficients[4]; // reflection coefficients
    int shine; // exponent term of specular component
    object() {}

    virtual void draw() {}

    void setColor(double color_r, double color_g, double color_b) {
        cout << "In set color" << endl;
        cout << color_r << color_g << color_b << endl << "\n";
        color[0] = color_r;
        color[1] = color_g;
        color[2] = color_b;
    }

    void setShine(int shine_value) {
        cout << "In shine value" << endl;
        cout << shine_value << endl << "\n";
        shine = shine_value;
    }

    void setCoEfficients(double ambient, double diffuse, double specular, double recur_reflect) {
        cout << "In set coeff" << endl;
        cout << ambient << diffuse << specular << recur_reflect << endl << "\n";
        coefficients[0] = ambient;
        coefficients[1] = diffuse;
        coefficients[2] = specular;
        coefficients[3] = recur_reflect;
    }

    virtual double intersect(Ray& r, vector<double>&final_color, int level) {
        return -1.0;
    }

    virtual Vector3D get_normal_vector() {}

    virtual double get_t_value() {}
};

class Sphere : public Object {

public:
    Sphere(Vector3D center, double radius) {
        reference_point = center;
        length = radius;
    }

    void draw() {
    int slices = 50;
    int stacks = 30;
    Vector3D points[100][100];
    int i, j;
    double h, r;
    //printVector3D(reference_point);
    //generate points
    for (i = 0; i <= stacks; i++) {
        h = length * sin(((double) i / (double) stacks) * (pi / 2));
        r = length * cos(((double) i / (double) stacks) * (pi / 2));
        for (j = 0; j <= slices; j++) {
            points[i][j].x = r * cos(((double) j / (double) slices) * 2 * pi) + reference_point.x;
            points[i][j].y = r * sin(((double) j / (double) slices) * 2 * pi) + reference_point.y;
            points[i][j].z = h;
        }
    }

    //glRotatef(90, 0, 1, 0);

    //glTranslatef(reference_point.x, reference_point.y, reference_point.z)
    //draw quads using generated points
    for (i = 0; i < stacks; i++) {

        for (j = 0; j < slices; j++) {
            glColor3f(color[0], color[1], color[2]);
            glBegin(GL_QUADS);
            {
                //lower hemisphere
                glVertex3f(points[i][j].x , points[i][j].y , -points[i][j].z + reference_point.z);
                glVertex3f(points[i][j + 1].x , points[i][j + 1].y , -points[i][j + 1].z + reference_point.z);
                glVertex3f(points[i + 1][j + 1].x , points[i + 1][j + 1].y , -points[i + 1][j + 1].z + reference_point.z);
                glVertex3f(points[i + 1][j].x, points[i + 1][j].y, -points[i + 1][j].z + reference_point.z);
            }
            glEnd();
        }
    }

    for (i = 0; i < stacks; i++) {
        for (j = 0; j < slices; j++) {
            glColor3f(color[0], color[1], color[2]);
            glBegin(GL_QUADS);
            {
                //upper hemisphere
                glVertex3f(points[i][j].x, points[i][j].y , points[i][j].z + reference_point.z);
                glVertex3f(points[i][j + 1].x , points[i][j + 1].y , points[i][j + 1].z + reference_point.z);
                glVertex3f(points[i + 1][j + 1].x , points[i + 1][j + 1].y, points[i + 1][j + 1].z +  reference_point.z);
                glVertex3f(points[i + 1][j].x , points[i + 1][j].y, points[i + 1][j].z + reference_point.z);
            }
            glEnd();
        }

    }
    }
    virtual double intersect(Ray& r, vector<double>&final_color, int level) {

        //printVector3D(r.dir);
        Vector3D R_start = {r.start.x - reference_point.x, r.start.y - reference_point.y, r.start.z - reference_point.z};
        //cout<<"R_START "<<endl;
        //printVector3D(R_start);
        //cout<<"R_DIR "<<endl;
        //printVector3D(r.dir);
        double a = 1;
        double b = 2 * getDotProduct(r.dir, R_start);
        //cout<<"b "<<b<<endl;
        double c = getDotProduct(R_start, R_start) - length*length;
        //cout<<"c "<<c<<endl;

        double temp = b*b - 4 * a * c;

        if(temp<0) {
            //cout<<"d<0"<<endl;
            return 1000000;
        }

        double d = sqrt(b*b - 4 * a * c);

        //cout<<"d "<<d<<endl;

        double t_pos = (-b+d)/(2*a);
        double t_neg = (-b-d)/(2*a);

        double final_t;

        if(t_pos > 0 && t_neg> 0){

            final_t = min(t_pos, t_neg);

        }else if(t_pos > 0 && t_neg < 0){

            final_t = t_pos;

        }else if(t_pos< 0 && t_neg > 0){

            final_t = t_neg;

        }else{
            return 1000000;
        }

        //cout<<"In sphere intersect"<<endl;

        // setting the color of the pixel of intersection
        final_color[0]=color[0]*1*coefficients[0];
        final_color[1]=color[1]*1*coefficients[0];
        final_color[2]=color[2]*1*coefficients[0];

        //cout<<final_color[0]<<final_color[1]<<final_color[2]<<endl;

        // returning the t_neg for now(not sure here
        return final_t;
    }

};

class Triangle : public Object {

public:
    Vector3D first_point;
    Vector3D second_point;
    Vector3D third_point;

    Triangle(Vector3D a, Vector3D b, Vector3D c) {
        first_point = a;
        second_point = b;
        third_point = c;
    }

    void draw() {
        glColor3f(color[0], color[1], color[2]);
        glBegin(GL_TRIANGLES);
            glVertex3f(first_point.x, first_point.y, first_point.z);
            glVertex3f(second_point.x, second_point.y, second_point.z);
            glVertex3f(third_point.x, third_point.y, third_point.z);
        glEnd();
    }
};

class GeneralQuadrates : public Object {
public:
    vector<double> eqn_coefficients;

    GeneralQuadrates(vector<double> coefficients, Vector3D cube_reference_point, double cube_l, double cube_w,
                     double cube_h) {
        eqn_coefficients = coefficients;
        reference_point = cube_reference_point;
        length = cube_l;
        width = cube_w;
        height = cube_h;
        printVector3D(reference_point);

    }
};

class Floor : public Object {
public:
    double num_of_tiles;
    Floor(double floorWidth, double tileWidth) {
        Vector3D reference_p;
        reference_p.x = -floorWidth / 2;
        reference_p.y = -floorWidth / 2;
        reference_p.z = 0;

        reference_point = Vector3D{-floorWidth / 2, -floorWidth / 2, 0};
        length = tileWidth;
        num_of_tiles = (int) floorWidth/tileWidth;
    }

    void draw() {


// write codes for drawing a checkerboard-like
// floor with alternate colors on adjacent tiles
    double x,y;
    int color = 0;
    for(int i =0; i<num_of_tiles;i++){

        y = reference_point.y + 20*i;
        //cout<<"y = "<<y<<endl;

        for(int j = 0;j<num_of_tiles;j++){

            x = reference_point.x + 20*j;

            //cout<<color<<endl;
            glColor3f(color,color,color);
            glBegin(GL_QUADS);
            {
                glVertex3f(x, y, 0);
                glVertex3f(x + 20, y, 0);
                glVertex3f(x+20, y+20, 0);
                glVertex3f(x, y+20, 0);
            }
            glEnd();

            color= 1- color;
        }
        color= 1- color;
        x = reference_point.x;
    }

    }

    virtual double intersect(Ray& r, vector<double>&final_color, int level) {

        Vector3D normal= {0,0,1};

        double t = getDotProduct(normal, r.start)/getDotProduct(normal, r.dir);
        t = t*-1;

        if(t < 0) return -1;

        return t;
    }


};

class Light {
public:
    Vector3D light_pos;
    double color[3];

    Light(Vector3D light_coord, double color_x, double color_y, double color_z) {
        light_pos = light_coord;
        color[0] = color_x;
        color[1] = color_y;
        color[2] = color_z;
    }
};
