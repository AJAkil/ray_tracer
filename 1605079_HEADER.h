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

double get_phong_intensity(double Ip, double kd, double ks, double object_color_comp,
                           int shine, Vector3D L, Vector3D N, Vector3D R, Vector3D V) {

    double diffused_comp = Ip * kd * object_color_comp * max(getDotProduct(L, N), 0.0);
    double specular_comp = Ip * ks  * max(pow(getDotProduct(R, V), shine), 0.0);

    //cout<<getDotProduct(L, N)<<" "<<getDotProduct(R, V);
    //cout<<specular_comp<<endl;

    return diffused_comp + specular_comp;

}

class Light {
public:
    Vector3D light_pos{};
    double color[3]{};

    Light(Vector3D light_coord, double color_x, double color_y, double color_z) {
        light_pos = light_coord;
        color[0] = color_x;
        color[1] = color_y;
        color[2] = color_z;
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
            h = 2 * sin(((double) i / (double) stacks) * (pi / 2));
            r = 2 * cos(((double) i / (double) stacks) * (pi / 2));
            for (j = 0; j <= slices; j++) {
                points[i][j].x = r * cos(((double) j / (double) slices) * 2 * pi) + light_pos.x;
                points[i][j].y = r * sin(((double) j / (double) slices) * 2 * pi) + light_pos.y;
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
                    glVertex3f(points[i][j].x, points[i][j].y, -points[i][j].z + light_pos.z);
                    glVertex3f(points[i][j + 1].x, points[i][j + 1].y, -points[i][j + 1].z + light_pos.z);
                    glVertex3f(points[i + 1][j + 1].x, points[i + 1][j + 1].y,
                               -points[i + 1][j + 1].z + light_pos.z);
                    glVertex3f(points[i + 1][j].x, points[i + 1][j].y, -points[i + 1][j].z + light_pos.z);
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
                    glVertex3f(points[i][j].x, points[i][j].y, points[i][j].z + light_pos.z);
                    glVertex3f(points[i][j + 1].x, points[i][j + 1].y, points[i][j + 1].z + light_pos.z);
                    glVertex3f(points[i + 1][j + 1].x, points[i + 1][j + 1].y,
                               points[i + 1][j + 1].z + light_pos.z);
                    glVertex3f(points[i + 1][j].x, points[i + 1][j].y, points[i + 1][j].z + light_pos.z);
                }
                glEnd();
            }

        }
    }
};


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
    Object() {}

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

    virtual double intersect(Ray &r, vector<double> &final_color, int level, int self_index) {
        return -1.0;
    }

    virtual Vector3D get_normal_vector(Vector3D poi) {}

    virtual double get_t_value() {}
};

extern vector<Object *> objects;
extern vector<Light> lights;

Vector3D calculate_color(double red, double green, double blue, Vector3D N, Vector3D V, Vector3D poi, int self_index,
                         int shine, double ka, double kd, double ks) {

    // by default we consider that the point is in shadows
    double Ir = red * 1 * ka;
    double Ig = green * 1 * ka;
    double Ib = blue * 1 * ka;

    // We loop over the light object
    for (int i = 0; i < lights.size(); i++) {

        Vector3D l_source = lights[i].light_pos;
        //cout<<"after getting light source"<<endl;

        // we form the light vector from light source to point of intersection and we get the length
        Vector3D L = {l_source.x - poi.x, l_source.y - poi.y, l_source.z - poi.z};
        double LR_length = getVectorMagnitude(L);

        //cout<<"after gettinng L and L.len()"<<endl;

        // we form the light ray
        Ray light_ray = Ray(poi, L);

        //cout<<"after getting light ray"<<endl;

        // we et a flag to check whether we have an object in between light ray and light source or not
        bool is_obstructed = false;

        /**
         * We run a for loop to check for each of the objects whether there is another object infront of the object
         * we are trying to set the color for. We recursively call the intersect function with level = 0. If
         * we obtain a t < LR_length, then it means the main object we are trying to color is being obstructed
         * by another object. Then we break the loop and set the is_obstructed flag to true
         */
        for (int j = 0; j < objects.size(); j++) {

            vector<double> dummy_color;
            dummy_color.push_back(0);
            dummy_color.push_back(0);
            dummy_color.push_back(0);

            //cout<<" before intersect call inside"<<endl;
            double t_test = objects[j]->intersect(light_ray, dummy_color, 0, j);
            //cout<<" after intersect call inside"<<endl;

            // we check the condition here
            if (j != self_index) {
                //cout<<" j!= self_index "<<endl;
                if (t_test < LR_length) {

                    //cout<<"obstruction condition"<<endl;
                    is_obstructed = true;
                    break;
                }
            }
        }

        // outside for loop of objects we set the colors as needed here

        if (!is_obstructed) {

            // we first normalize L
            L = getUnitVector(L);

            // calculate R
            double scaler = 2 * getDotProduct(L, N);
            Vector3D R = {N.x * scaler - L.x, N.y * scaler - L.y, N.z * scaler - L.z};
            R = getUnitVector(R); // normalize R

            //ambient, diffused, specular, recursive
            Ir += get_phong_intensity(lights[i].color[0], kd, ks,
                                      red, shine, L, N, R, V * -1);

            Ig += get_phong_intensity(lights[i].color[1], kd, ks,
                                      green, shine, L, N, R, V * -1);

            Ib += get_phong_intensity(lights[i].color[2], kd, ks,
                                      blue, shine, L, N, R, V * -1);

        }

    }

    Vector3D light_comps = {Ir, Ig, Ib};

    return light_comps;
}


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
                    glVertex3f(points[i][j].x, points[i][j].y, -points[i][j].z + reference_point.z);
                    glVertex3f(points[i][j + 1].x, points[i][j + 1].y, -points[i][j + 1].z + reference_point.z);
                    glVertex3f(points[i + 1][j + 1].x, points[i + 1][j + 1].y,
                               -points[i + 1][j + 1].z + reference_point.z);
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
                    glVertex3f(points[i][j].x, points[i][j].y, points[i][j].z + reference_point.z);
                    glVertex3f(points[i][j + 1].x, points[i][j + 1].y, points[i][j + 1].z + reference_point.z);
                    glVertex3f(points[i + 1][j + 1].x, points[i + 1][j + 1].y,
                               points[i + 1][j + 1].z + reference_point.z);
                    glVertex3f(points[i + 1][j].x, points[i + 1][j].y, points[i + 1][j].z + reference_point.z);
                }
                glEnd();
            }

        }
    }

    virtual double intersect(Ray &r, vector<double> &final_color, int level, int self_index) {

        //printVector3D(r.dir);
        Vector3D R_start = {r.start.x - reference_point.x, r.start.y - reference_point.y,
                            r.start.z - reference_point.z};
        //cout<<"R_START "<<endl;
        //printVector3D(R_start);
        //cout<<"R_DIR "<<endl;
        //printVector3D(r.dir);
        double a = 1;
        double b = 2 * getDotProduct(r.dir, R_start);
        //cout<<"b "<<b<<endl;
        double c = getDotProduct(R_start, R_start) - length * length;
        //cout<<"c "<<c<<endl;

        double temp = b * b - 4 * a * c;

        if (temp < 0) {
            //cout<<"d<0"<<endl;
            return 1000000;
        }

        double d = sqrt(b * b - 4 * a * c);

        //cout<<"d "<<d<<endl;

        double t_pos = (-b + d) / (2 * a);
        double t_neg = (-b - d) / (2 * a);

        double final_t;

        if (t_pos > 0 && t_neg > 0) {

            final_t = min(t_pos, t_neg);

        } else if (t_pos > 0 && t_neg < 0) {

            final_t = t_pos;

        } else if (t_pos < 0 && t_neg > 0) {

            final_t = t_neg;

        } else {
            return 1000000;
        }

        //cout<<"In sphere intersect"<<endl;

        if(level == 0) return final_t;

        Vector3D poi = {r.start.x + final_t * r.dir.x, r.start.y + final_t * r.dir.y, r.start.z + final_t * r.dir.z};

        Vector3D N = get_normal_vector(poi);
        N = getUnitVector(N);
//
        Vector3D light_components = calculate_color(color[0], color[1], color[2], N, r.dir, poi, self_index, shine,
                                                    coefficients[0],
                                                    coefficients[1], coefficients[2]);
//
        final_color[0] = light_components.x;
        final_color[1] = light_components.y;
        final_color[2] = light_components.z;



        // setting the color of the pixel of intersection
//        final_color[0] = color[0] * 1 * coefficients[0];
//        final_color[1] = color[1] * 1 * coefficients[0];
//        final_color[2] = color[2] * 1 * coefficients[0];

        //cout<<final_color[0]<<final_color[1]<<final_color[2]<<endl;

        return final_t;
    }

    virtual Vector3D get_normal_vector(Vector3D poi) {

        Vector3D normal = {poi.x - reference_point.x, poi.y - reference_point.y, poi.z - reference_point.z};

        return normal;

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

    virtual double intersect(Ray &r, vector<double> &final_color, int level, int self_index) {

        double beta, gamma, t, A, D1, D2, D3, x_o, y_o, z_o, x_d, y_d, z_d;

        // setting up the ray components
        x_o = r.start.x;
        y_o = r.start.y;
        z_o = r.start.z;

        x_d = r.dir.x;
        y_d = r.dir.y;
        z_d = r.dir.z;

        Vector3D a = first_point;
        Vector3D b = second_point;
        Vector3D c = third_point;

        // calculate A
        A = (a.x - b.x) * ((a.y - c.y) * z_d - y_d * (a.z - c.z)) -
            (a.x - c.x) * ((a.y - b.y) * z_d - (a.z - b.z) * y_d) +
            x_d * ((a.y - b.y) * (a.z - c.z) - (a.y - c.y) * (a.z - b.z));

        // calculate D1
        D1 = (a.x - x_o) * ((a.y - c.y) * z_d - y_d * (a.z - c.z)) -
             (a.x - c.x) * ((a.y - y_o) * z_d - (a.z - z_o) * y_d) +
             x_d * ((a.y - y_o) * (a.z - c.z) - (a.y - c.y) * (a.z - z_o));

        // calculate D2
        D2 = (a.x - b.x) * ((a.y - y_o) * z_d - y_d * (a.z - z_o)) -
             (a.x - x_o) * ((a.y - b.y) * z_d - (a.z - b.z) * y_d) +
             x_d * ((a.y - b.y) * (a.z - z_o) - (a.y - y_o) * (a.z - b.z));

        // calculate D3
        D3 = (a.x - b.x) * ((a.y - c.y) * (a.z - z_o) - (a.y - y_o) * (a.z - c.z)) -
             (a.x - c.x) * ((a.y - b.y) * (a.z - z_o) - (a.z - b.z) * (a.y - y_o)) +
             (a.x - x_o) * ((a.y - b.y) * (a.z - c.z) - (a.y - c.y) * (a.z - b.z));

        beta = D1 / A;
        gamma = D2 / A;
        t = D3 / A;

        if (beta + gamma < 1 && beta > 0 && gamma > 0 && t > 0) {

            if(level == 0 ) return t;

            Vector3D poi = {r.start.x + t * r.dir.x, r.start.y + t * r.dir.y, r.start.z + t * r.dir.z};

            //cout<<"before normal calculation"<<endl;
            Vector3D N = get_normal_vector(poi);
            N = getUnitVector(N);
            //cout<<"after normal calculation"<<endl;
//
            Vector3D light_components = calculate_color(color[0], color[1], color[2], N, r.dir, poi, self_index, shine,
                                                        coefficients[0],
                                                        coefficients[1], coefficients[2]);
//
            final_color[0] = light_components.x;
            final_color[1] = light_components.y;
            final_color[2] = light_components.z;

//            final_color[0] = color[0] * 1 * coefficients[0];
//            final_color[1] = color[1] * 1 * coefficients[0];
//            final_color[2] = color[2] * 1 * coefficients[0];

            // returning the final_t
            return t;

        } else {

            return 1000000;
        }

    }

    virtual Vector3D get_normal_vector(Vector3D poi) {

        Vector3D side_a = {second_point.x - first_point.x, second_point.y - first_point.y,
                           second_point.z - first_point.z};
        Vector3D side_b = {third_point.x - first_point.x, third_point.y - first_point.y,
                           third_point.z - first_point.z};

        return getCrossProduct(side_a, side_b);

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

        cout << "coeff" << eqn_coefficients.size() << endl;

    }

    bool is_inside_cube(Vector3D point) {

        bool flag = true;

        if (length != 0) {

            if (reference_point.x <= point.x && point.x <= reference_point.x + length) {
                flag = true;
            } else {
                return false;
            }

        }

        if (width != 0) {

            if (reference_point.y <= point.y && point.y <= reference_point.y + width) {
                flag = true;
            } else {
                return false;
            }

        }

        if (height != 0) {

            if (reference_point.z <= point.z && point.z <= reference_point.z + height) {
                flag = true;
            } else {
                return false;
            }

        }

        return flag;

    }


    virtual double intersect(Ray &r, vector<double> &final_color, int level, int self_index) {

        double x_o, y_o, z_o, x_d, y_d, z_d;

        // setting up the coefficients
        double A = eqn_coefficients[0];
        double B = eqn_coefficients[1];
        double C = eqn_coefficients[2];
        double D = eqn_coefficients[3];
        double E = eqn_coefficients[4];
        double F = eqn_coefficients[5];
        double G = eqn_coefficients[6];
        double H = eqn_coefficients[7];
        double I = eqn_coefficients[8];
        double J = eqn_coefficients[9];

        //cout<<A<<" "<<B<<" "<<C<<" "<<D<<" "<<E<<" "<<F<<" "<<G<<" "<<H<<" "<<I<<" "<<J<<endl;

        // setting up the ray components
        x_o = r.start.x;
        y_o = r.start.y;
        z_o = r.start.z;

        x_d = r.dir.x;
        y_d = r.dir.y;
        z_d = r.dir.z;


        // we calculate a,b,c and find d
        double a = A * x_d * x_d + B * y_d * y_d + C * z_d * z_d + D * x_d * y_d + E * x_d * z_d + F * y_d * z_d;

        double b = 2 * A * x_o * x_d + 2 * B * y_o * y_d + 2 * C * z_o * z_d + D * x_o * y_d + D * y_o * x_d +
                   E * x_o * z_d + E * z_o * x_d + F * y_o * z_d + F * z_o * y_d + G * x_d + H * y_d + I * z_d;

        double c = A * x_o * x_o + B * y_o * y_o + C * z_o * z_o + D * x_o * y_o + E * x_o * z_o + F * y_o * z_o +
                   G * x_o + H * y_o + I * z_o + J;

        // we first calculate the discriminant to check for non intersection
        double temp = b * b - 4 * a * c;

        if (temp < 0) {
            //cout<<"d<0"<<endl;
            return 1000000;
        }

        double d = sqrt(b * b - 4 * a * c);

        //cout<<"d "<<d<<endl;

        // we compute the value of t here
        double t_pos = (-b + d) / (2 * a);
        double t_neg = (-b - d) / (2 * a);

        // we check if the point of intersection is within the bounding cube or not

        // we first find the point of intersections
        Vector3D poi_pos = {r.start.x + t_pos * r.dir.x, r.start.y + t_pos * r.dir.y, r.start.z + t_pos * r.dir.z};
        Vector3D poi_neg = {r.start.x + t_neg * r.dir.x, r.start.y + t_neg * r.dir.y, r.start.z + t_neg * r.dir.z};

        // if t_neg > 0 and point 2 is valid return t_neg --> else same for t_pos > 0 and point 1 is valid ---> non intersection

        if (t_neg > 0 && is_inside_cube(poi_neg)) {

            if(level == 0 ) return t_neg;

            Vector3D poi = {r.start.x + t_neg * r.dir.x, r.start.y + t_neg * r.dir.y, r.start.z + t_neg * r.dir.z};

            //cout<<"before normal calculation"<<endl;
            Vector3D N = get_normal_vector(poi);
            N = getUnitVector(N);
            //cout<<"after normal calculation"<<endl;
//
            Vector3D light_components = calculate_color(color[0], color[1], color[2], N, r.dir, poi, self_index, shine,
                                                        coefficients[0],
                                                        coefficients[1], coefficients[2]);
//
            final_color[0] = light_components.x;
            final_color[1] = light_components.y;
            final_color[2] = light_components.z;

            // setting the color of the pixel of intersection
//            final_color[0] = color[0] * 1 * coefficients[0];
//            final_color[1] = color[1] * 1 * coefficients[0];
//            final_color[2] = color[2] * 1 * coefficients[0];

            return t_neg;

        } else if (t_pos > 0 && is_inside_cube(poi_pos)) {

            if(level == 0) return t_pos;


            Vector3D poi = {r.start.x + t_pos * r.dir.x, r.start.y + t_pos * r.dir.y, r.start.z + t_pos * r.dir.z};

            //cout<<"before normal calculation"<<endl;
            Vector3D N = get_normal_vector(poi);
            N = getUnitVector(N);
            //cout<<"after normal calculation"<<endl;
//
            Vector3D light_components = calculate_color(color[0], color[1], color[2], N, r.dir, poi, self_index, shine,
                                                        coefficients[0],
                                                        coefficients[1], coefficients[2]);
//
            final_color[0] = light_components.x;
            final_color[1] = light_components.y;
            final_color[2] = light_components.z;

            // setting the color of the pixel of intersection
//            final_color[0] = color[0] * 1 * coefficients[0];
//            final_color[1] = color[1] * 1 * coefficients[0];
//            final_color[2] = color[2] * 1 * coefficients[0];

            return t_pos;

        } else {
            return 10000000;
        }

        // then we check if the poi2 and poi2 is within cube
        /*if ((reference_point.x <= poi1.x && poi1.x <= reference_point.x + length) &&
            (reference_point.y <= poi1.y && poi1.y <= reference_point.y + width) &&
            (reference_point.z <= poi1.z && poi1.z <= reference_point.z + height)
            && (reference_point.x <= poi2.x && poi2.x <= reference_point.x + length) &&
            (reference_point.y <= poi2.y && poi2.y <= reference_point.y + width) &&
            (reference_point.z <= poi2.z && poi2.z <= reference_point.z + height)) {

            // we choose the appropriate t from here and then we set the color here
            cout << "FIRST CASE" << endl;
            double final_t;

            if (t_pos > 0 && t_neg > 0) {

                final_t = min(t_pos, t_neg);

            } else if (t_pos > 0 && t_neg < 0) {

                final_t = t_pos;

            } else if (t_pos < 0 && t_neg > 0) {

                final_t = t_neg;

            } else {
                cout << "No two point" << endl;
                return 1000000;
            }


            // we set the color and then we return the final t

            // setting the color of the pixel of intersection
            final_color[0] = color[0] * 1 * coefficients[0];
            final_color[1] = color[1] * 1 * coefficients[0];
            final_color[2] = color[2] * 1 * coefficients[0];

            //cout<<final_color[0]<<final_color[1]<<final_color[2]<<endl;

            // returning the final_t
            cout << "In general " << final_t << endl;
            return final_t;

        } else if ((reference_point.x <= poi1.x && poi1.x <= reference_point.x + length) &&
                   (reference_point.y <= poi1.y && poi1.y <= reference_point.y + width) &&
                   (reference_point.z <= poi1.z && poi1.z <= reference_point.z + height)) {

            // we check for the t_pos here
            cout << "SECOND CASE" << endl;

            if (t_pos > 0) {


                // setting the color of the pixel of intersection
                final_color[0] = color[0] * 1 * coefficients[0];
                final_color[1] = color[1] * 1 * coefficients[0];
                final_color[2] = color[2] * 1 * coefficients[0];

                //cout<<final_color[0]<<final_color[1]<<final_color[2]<<endl;

                // returning the final_t
                cout << "in general " << t_pos << endl;
                return t_pos;

            } else {

                cout << "No t_pos" << endl;
                return 1000000;

            }

        } else if ((reference_point.x <= poi2.x && poi2.x <= reference_point.x + length) &&
                   (reference_point.y <= poi2.y && poi2.y <= reference_point.y + width) &&
                   (reference_point.z <= poi2.z && poi2.z <= reference_point.z + height)) {

            cout << "THIRD CASE" << endl;
            if (t_neg > 0) {


                // setting the color of the pixel of intersection
                final_color[0] = color[0] * 1 * coefficients[0];
                final_color[1] = color[1] * 1 * coefficients[0];
                final_color[2] = color[2] * 1 * coefficients[0];

                //cout<<final_color[0]<<final_color[1]<<final_color[2]<<endl;

                // returning the final_t
                cout << "in general " << t_neg;
                return t_neg;

            } else {

                cout << "No t_neg" << endl;
                return 1000000;

            }

        } else {
            //cout<<"FOURTH CASE"<<endl;
            //cout<<"No ultimate"<<endl;
            return 1000000;
        }*/

    }

    virtual Vector3D get_normal_vector(Vector3D poi) {

        Vector3D normal = {0.0, 0.0, 0.0};

        // setting up the coefficients
        double A = eqn_coefficients[0];
        double B = eqn_coefficients[1];
        double C = eqn_coefficients[2];
        double D = eqn_coefficients[3];
        double E = eqn_coefficients[4];
        double F = eqn_coefficients[5];
        double G = eqn_coefficients[6];
        double H = eqn_coefficients[7];
        double I = eqn_coefficients[8];
        double J = eqn_coefficients[9];

        double x = 2 * A * poi.x + D * poi.y + E * poi.z + G;
        double y = 2 * B * poi.y + D * poi.x + F * poi.z + H;
        double z = 2 * C * poi.z + E * poi.x + F * poi.y + I;

        normal = {x, y, z};

        return normal;

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
        num_of_tiles = (int) floorWidth / tileWidth;
    }

    void draw() {
// write codes for drawing a checkerboard-like
// floor with alternate colors on adjacent tiles
        double x, y;
        int color = 0;
        for (int i = 0; i < num_of_tiles; i++) {

            y = reference_point.y + 20 * i;
            //cout<<"y = "<<y<<endl;

            for (int j = 0; j < num_of_tiles; j++) {

                x = reference_point.x + 20 * j;

                //cout<<color<<endl;
                glColor3f(color, color, color);
                glBegin(GL_QUADS);
                {
                    glVertex3f(x, y, 0);
                    glVertex3f(x + 20, y, 0);
                    glVertex3f(x + 20, y + 20, 0);
                    glVertex3f(x, y + 20, 0);
                }
                glEnd();

                color = 1 - color;
            }
            color = 1 - color;
            x = reference_point.x;
        }

    }


    virtual double intersect(Ray &r, vector<double> &final_color, int level, int self_index) {

        Vector3D normal = {0, 0, 1};

        double t = getDotProduct(normal, r.start) / getDotProduct(normal, r.dir);
        t = t * -1;

        if (t < 0) return 1000000;

        Vector3D poi = {r.start.x + t * r.dir.x, r.start.y + t * r.dir.y, r.start.z + t * r.dir.z};

        // checking to see if the poi is within the floor itself
        //cout<<reference_point.x<<" " <<reference_point.y<<endl;
        if ((reference_point.x <= poi.x && poi.x <= -1 * reference_point.x) &&
            (reference_point.y <= poi.y && poi.y <= -1 * reference_point.y)) {

            // setting the color of the pixel of intersection

            if (level == 0) return t;

            // if level is not 0 then calculate the required lighting components here

            // first we calculate the normal function here and normalize the vector
            Vector3D N = getUnitVector(normal);
            double Ir, Ig, Ib;

            int x = (poi.x - reference_point.x) / length;
            int y = (poi.y - reference_point.y) / length;

            double red = 0, green = 0, blue = 0;

            if ((x + y) % 2 == 0) {

                red = 0;
                green = 0;
                blue = 0;

            } else {

                red = 1;
                green = 1;
                blue = 1;

            }

            // by default we consider that the point is in shadows
            /*Ir = red * 1 * coefficients[0];
            Ig = green * 1 * coefficients[0];
            Ib = blue * 1 * coefficients[0];

            // We loop over the light object
            for (int i = 0; i < lights.size(); i++) {

                Vector3D l_source = lights[i].light_pos;

                // we form the light vector from light source to point of intersection and we get the length
                Vector3D L = {l_source.x - poi.x, l_source.y - poi.y, l_source.z - poi.z};
                double LR_length = getVectorMagnitude(L);

                // we form the light ray
                Ray light_ray = Ray(poi, L);

                // we et a flag to check whether we have an object in between light ray and light source or not
                bool is_obstructed = false;

                *//**
                 * We run a for loop to check for each of the objects whether there is another object infront of the object
                 * we are trying to set the color for. We recursively call the intersect function with level = 0. If
                 * we obtain a t < LR_length, then it means the main object we are trying to color is being obstructed
                 * by another object. Then we break the loop and set the is_obstructed flag to true
                 *//*
                for (int j = 0; j < objects.size(); j++) {

                    vector<double> dummy_color;
                    dummy_color.push_back(0);
                    dummy_color.push_back(0);
                    dummy_color.push_back(0);

                    double t_test = objects[j]->intersect(light_ray, dummy_color, 0, j);

                    // we check the condition here
                    if (j != self_index) {
                        if (t_test < LR_length) {
                            is_obstructed = true;
                            break;
                        }
                    }
                }

                // outside for loop of objects we set the colors as needed here

                if (!is_obstructed) {

                    // we first normalize L
                    L = getUnitVector(L);

                    // calculate R
                    double scaler = 2 * getDotProduct(L, N);
                    Vector3D R = {N.x * scaler - L.x, N.y * scaler - L.y, N.z * scaler - L.z};
                    R = getUnitVector(R); // normalize R

                    //ambient, diffused, specular, recursive
                    Ir += get_phong_intensity(lights[i].color[0], coefficients[1], coefficients[2],
                                              red, shine, L, N, R, r.dir * -1);

                    Ig += get_phong_intensity(lights[i].color[1], coefficients[1], coefficients[2],
                                              green, shine, L, N, R, r.dir * -1);

                    Ib += get_phong_intensity(lights[i].color[2], coefficients[1], coefficients[2],
                                              blue, shine, L, N, R, r.dir * -1);

                }

            }


            final_color[0] = Ir;
            final_color[1] = Ig;
            final_color[2] = Ib;*/

            Vector3D light_components = calculate_color(red, green, blue, N, r.dir, poi, self_index, shine,
                                                        coefficients[0],
                                                        coefficients[1], coefficients[2]);

            final_color[0] = light_components.x;
            final_color[1] = light_components.y;
            final_color[2] = light_components.z;

            return t;


        } else return 1000000;


    }


};
