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

    virtual double intersect(Ray *r, double *color, int level) {
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

    void draw() {}
};

class Floor : public Object {
public:
    Floor(double floorWidth, double tileWidth) {
        Vector3D reference_p;
        reference_p.x = -floorWidth / 2;
        reference_p.y = -floorWidth / 2;
        reference_p.z = 0;

        reference_point = Vector3D{-floorWidth / 2, -floorWidth / 2, 0};
        length = tileWidth;
    }

    void draw() {
// write codes for drawing a checkerboard-like
// floor with alternate colors on adjacent tiles
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
