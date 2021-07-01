#include "1605079_Header.h"
#include "bitmap_image.hpp"


using namespace std;

#define pi (2*acos(0.0))

double cameraHeight;
double cameraAngle;
int drawaxes;
double angle;
double camTheta = 0.75;


Vector3D pos = {100, 100, 0};
Vector3D u = {0, 0, 1};
Vector3D l = {-1 / sqrt(2), -1 / sqrt(2), 0};
Vector3D r = {-1 / sqrt(2), 1 / sqrt(2), 0};
double k = 2;
vector<string> fileLines;
vector<Object *> objects;
vector<Light> lights;
int recursion_level, image_width, image_height;
double windowHeight = 500, windowWidth = 500;
int num_objects, num_light_sources;
double view_angle = 80;

vector<string> tokenizeString(string s) {
    vector<string> tokens;
    //cout << s << endl;

    // create the stream
    stringstream tokenizer(s);
    string token;

    while (getline(tokenizer, token, ' ')) {
        tokens.push_back(token);
    }

    //cout << tokens[0] << endl;

    return tokens;
}

Vector3D parseLine(string line) {
    //cout << "In here" << endl;
    vector<string> lines;
    Vector3D a = {0, 0, 0};

    lines = tokenizeString(line);

    istringstream(lines[0]) >> a.x;
    istringstream(lines[1]) >> a.y;
    istringstream(lines[2]) >> a.z;

    return a;
}

void readFile(const char *fileName) {

    fileLines.clear();
    fstream newfile;

    newfile.open(fileName, ios::in);
    if (newfile.is_open()) {
        string line;
        while (getline(newfile, line)) { //read data from file object and put it into string.
            fileLines.push_back(line);
        }
        newfile.close(); //close the file object.
    }

//    for (int i = 0; i < fileLines.size(); i++) {
//        cout << fileLines[i] << endl;
//    }


}

double convertDegreeToRadian(int modifier) {
    return modifier * camTheta * pi / 180;
}

void camLookLr(int modifier) {
    Vector3D temp_vec;

    temp_vec = r;
    double theta = convertDegreeToRadian(modifier);
    r = r * cos(theta) + l * sin(theta);
    l = l * cos(theta) + (temp_vec * -1) * sin(theta);
}

void camLookUpDown(int modifier) {
    Vector3D temp_vec;

    temp_vec = l;
    double theta = convertDegreeToRadian(modifier);
    l = l * cos(theta) + u * sin(theta);
    u = u * cos(theta) + (temp_vec * -1) * sin(theta);
}

void camTilt(int modifier) {
    Vector3D temp_vec;

    temp_vec = r;
    double theta = convertDegreeToRadian(modifier);
    r = r * cos(theta) + u * sin(theta);
    u = u * cos(theta) + (temp_vec * -1) * sin(theta);
}


void drawAxes() {
    if (drawaxes == 1) {
        glColor3f(1.0, 1.0, 1.0);
        glBegin(GL_LINES);
        {
            glColor3f(1.0, 0, 0);
            glVertex3f(1500, 0, 0);
            glVertex3f(-1500, 0, 0);

            glColor3f(0, 1.0, 0);
            glVertex3f(0, -1500, 0);
            glVertex3f(0, 1500, 0);

            glColor3f(0, 0, 1.0);
            glVertex3f(0, 0, 1500);
            glVertex3f(0, 0, -1500);
        }
        glEnd();
    }
}


void capture() {

    //cout << "number of objects: " << objects.size() << endl;

    bitmap_image image((int) image_width, (int) image_height);

    // setting up the image
    //cout << "calling capture" << endl;
    for (int i = 0; i < image_width; i++) {
        for (int j = 0; j < image_height; j++) {
            image.set_pixel(i, j, 0, 0, 0);
        }
    }


    double planeDistance = (windowHeight / 2.0) / tan(convertDegreeToRadian(view_angle) / 2.0);
    Vector3D topleft = pos + l * planeDistance - r * (windowWidth / 2) + u * (windowHeight / 2);
    double du = windowWidth / image_width;
    double dv = windowHeight / image_height;

    //cout << planeDistance << " " << du << " " << dv << endl;
    //printVector3D(topleft);

    // Choose middle of the grid cell
    topleft = topleft + r * (0.5 * du) - u * (0.5 * dv);
    //printVector3D(topleft);

    // the index of the nearest object dummy_color


    for (int i = 0; i < image_width; i++) {
        for (int j = 0; j < image_height; j++) {


            int nearest = -1;
            double t_min = 10000, t;

            // calculate current pixel vector/point
            Vector3D currentPixel = topleft + r * (i * du) - u * (j * dv);
            Vector3D R_d = {currentPixel.x - pos.x, currentPixel.y - pos.y, currentPixel.z - pos.z};

            Ray ray(pos, R_d);
            vector<double> dummy_color;
            dummy_color.push_back(0);
            dummy_color.push_back(0);
            dummy_color.push_back(0);
            //cout<<getVectorMagnitude(r.dir);

            // iterating over the objects
            for (int k = 0; k < objects.size(); k++) {

                //cout<<"for pixel i and j" << i << " " << j<<endl;
                t = objects[k]->intersect(ray, dummy_color, 0, k);
                if (t < t_min) {
                    nearest = k; //storing the index of the nearest object
                    t_min = t;
                    //cout<<"storing "<<k<<endl;
                }
            }

            // we again call intersect to set the color
            vector<double> color;
            color.push_back(0);
            color.push_back(0);
            color.push_back(0);

            if (nearest != -1) {
                double temp = objects[nearest]->intersect(ray, color, 1, nearest);

                // we set image pixel here
                // we clip colors here before we set the pixels here
                if (color[0] > 1) {
                    color[0] = 1;
                } else if (color[0] < 0) {
                    color[0] = 0;
                }

                if (color[1] > 1) {
                    color[1] = 1;
                } else if (color[1] < 0) {
                    color[1] = 0;
                }

                if (color[2] > 1) {
                    color[2] = 1;
                } else if (color[2] < 0) {
                    color[2] = 0;
                }


                image.set_pixel(i, j, color[0] * 255, color[1] * 255, color[2] * 255);
            }

            color.clear();
            dummy_color.clear();

        }
    }

    image.save_image(
            "C:\\Users\\ragnarok_79\\Documents\\Important Docs\\Academics\\4_1\\Lab\\Graphics\\Ray-Tracer\\ray_tracer\\1605079_out.bmp");
}

void keyboardListener(unsigned char key, int x, int y) {

    switch (key) {

        case '0':
            capture();
            break;
        case '1':
            camLookLr(1);
            break;
        case '2':
            camLookLr(-1);
            break;
        case '3':
            camLookUpDown(1);
            break;
        case '4':
            camLookUpDown(-1);
            break;
        case '5':
            camTilt(1);
            break;
        case '6':
            camTilt(-1);
            break;
        default:
            break;
    }
}

void specialKeyListener(int key, int x, int y) {
    switch (key) {
        case GLUT_KEY_DOWN:
            //down arrow key
            pos.x -= k * l.x;
            pos.y -= k * l.y;
            pos.z -= k * l.z;
            break;
        case GLUT_KEY_UP:        // up arrow key
            //cameraHeight += 3.0;
            pos.x += k * l.x;
            pos.y += k * l.y;
            pos.z += k * l.z;
            break;

        case GLUT_KEY_RIGHT:
            pos.x += k * r.x;
            pos.y += k * r.y;
            pos.z += k * r.z;
            break;
        case GLUT_KEY_LEFT:
            pos.x -= k * r.x;
            pos.y -= k * r.y;
            pos.z -= k * r.z;
            break;

        case GLUT_KEY_PAGE_UP:
            pos.x += k * u.x;
            pos.y += k * u.y;
            pos.z += k * u.z;
            break;
        case GLUT_KEY_PAGE_DOWN:
            pos.x -= k * u.x;
            pos.y -= k * u.y;
            pos.z -= k * u.z;

        default:
            break;
    }
}


void loadData() {
    //cout << "out" << endl;
    readFile(
            "C:\\Users\\ragnarok_79\\Documents\\Important Docs\\Academics\\4_1\\Lab\\Graphics\\Ray-Tracer\\ray_tracer\\scene.txt");

    // Reading in the data
    vector<string> lines;
    lines = tokenizeString(fileLines[0]);
    istringstream(lines[0]) >> recursion_level;

    //cout << recursion_level << endl;

    lines = tokenizeString(fileLines[1]);
    istringstream(lines[0]) >> image_width;
    istringstream(lines[0]) >> image_height;

    //cout << "EKHNE" << recursion_level << " " << image_width << " " << image_height << endl;

    for (int i = 4; i < fileLines.size(); i++) {
        Object *temp;
        if (fileLines[i] == "sphere") {
            //cout<<"In sphere"<<endl;
            //cout << fileLines[i] << i << endl;

            Vector3D center = parseLine(fileLines[i + 1]);
            double radius;
            istringstream(fileLines[i + 2]) >> radius;

            temp = new Sphere(center, radius);

            Vector3D color = parseLine(fileLines[i + 3]);
            temp->setColor(color.x, color.y, color.z);

            vector<string> lines;
            lines = tokenizeString(fileLines[i + 4]);
            double coeff1, coeff2, coeff3, coeff4, shine;
            istringstream(lines[0]) >> coeff1;
            istringstream(lines[1]) >> coeff2;
            istringstream(lines[2]) >> coeff3;
            istringstream(lines[3]) >> coeff4;

            temp->setCoEfficients(coeff1, coeff2, coeff3, coeff4);

            istringstream(fileLines[i + 5]) >> shine;
            temp->setShine(shine);

            objects.push_back(temp);
            i += 5;

        } else if (fileLines[i] == "triangle") {

            Vector3D a = parseLine(fileLines[i + 1]);
            Vector3D b = parseLine(fileLines[i + 2]);
            Vector3D c = parseLine(fileLines[i + 3]);

            temp = new Triangle(a, b, c);

            Vector3D color = parseLine(fileLines[i + 4]);
            temp->setColor(color.x, color.y, color.z);

            vector<string> lines;
            lines = tokenizeString(fileLines[i + 5]);
            double coeff1, coeff2, coeff3, coeff4, shine;
            istringstream(lines[0]) >> coeff1;
            istringstream(lines[1]) >> coeff2;
            istringstream(lines[2]) >> coeff3;
            istringstream(lines[3]) >> coeff4;

            temp->setCoEfficients(coeff1, coeff2, coeff3, coeff4);

            istringstream(fileLines[i + 6]) >> shine;
            temp->setShine(shine);

            objects.push_back(temp);

            i += 6;

        } else if (fileLines[i] == "general") {

            //cout << "In General" << endl;
            vector<string> lines;
            vector<double> coeffs;

            lines = tokenizeString(fileLines[i + 1]);

            // reading in the coefficients from the file
            for (int i = 0; i < lines.size(); i++) {
                double t;
                istringstream(lines[i]) >> t;
                coeffs.push_back(t);
                //cout<<i<<" ";
            }

//            for (int i = 0; i < coeffs.size(); i++) {
//
//                cout << coeffs[i] << endl;
//            }

            // we read in the cube data
            lines = tokenizeString(fileLines[i + 2]);
            Vector3D cube_reference_point = parseLine(lines[0] + " " + lines[1] + " " + lines[2]);
            double c_l, c_w, c_h;
            istringstream(lines[3]) >> c_l;
            istringstream(lines[4]) >> c_w;
            istringstream(lines[5]) >> c_h;

            temp = new GeneralQuadrates(coeffs, cube_reference_point, c_l, c_w, c_h);

            // Setting the color
            Vector3D color = parseLine(fileLines[i + 3]);
            temp->setColor(color.x, color.y, color.z);

            // setting the coefficients
            lines = tokenizeString(fileLines[i + 4]);
            double coeff1, coeff2, coeff3, coeff4, shine;
            istringstream(lines[0]) >> coeff1;
            istringstream(lines[1]) >> coeff2;
            istringstream(lines[2]) >> coeff3;
            istringstream(lines[3]) >> coeff4;

            temp->setCoEfficients(coeff1, coeff2, coeff3, coeff4);

            istringstream(fileLines[i + 5]) >> shine;
            temp->setShine(shine);

            objects.push_back(temp);

            i += 5;
        } else if (fileLines[i].empty()) {
            //cout<<"here?"<<endl;
            continue;
        } else {
            //cout<<"In light section"<<endl;

            double num_lights;
            istringstream(fileLines[i]) >> num_lights;
            //cout<<i+1+(int)num_lights*2<<endl;

            for (int j = i + 1; j < i + 1 + (int) num_lights * 2; j += 2) {
                //cout << fileLines[j] << endl;
                //cout << fileLines[j + 1] << endl;

                Vector3D light_pos = parseLine(fileLines[j]);
                //printVector3D(light_pos);
                Vector3D color = parseLine(fileLines[j + 1]);
                Light light_source(light_pos, color.x, color.y, color.z);
                lights.push_back(light_source);
            }

            i += 1 + (int) num_lights * 2;

        }
        //cout<<fileLines[i]<<endl;
    }

    // code for creating the floor object
    Object *temp;
    temp = new Floor(1000, 20); // you can change these values

    // Setting the color
    temp->setColor(1, 1, 1);
    temp->setCoEfficients(0.4, 0.2, 0.3, 0.3);
    temp->setShine(5);

    objects.push_back(temp);

    //cout << objects.size() << endl;


    /*for (int i = 0; i < objects.size(); i++) {
       cout << objects[i]->length << endl;
    }*/

    /*for(auto & light : lights){
        cout<<light.light_pos.x<<endl;
    }*/


}

void display() {

    //clear the display
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    glClearColor(0, 0, 0, 0);    //color black
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

    /********************
    / set-up camera here
    ********************/
    //load the correct matrix -- MODEL-VIEW matrix
    glMatrixMode(GL_MODELVIEW);

    //initialize the matrix
    glLoadIdentity();

    //now give three info
    //1. where is the camera (viewer)?
    //2. where is the camera looking?
    //3. Which direction is the camera's UP direction?

    //gluLookAt(100,100,100,	0,0,0,	0,0,1);
    //gluLookAt(200*cos(cameraAngle), 200*sin(cameraAngle), cameraHeight,		0,0,0,		0,0,1);
    gluLookAt(pos.x, pos.y, pos.z, pos.x + l.x, pos.y + l.y, pos.z + l.z, u.x, u.y, u.z);


    //again select MODEL-VIEW
    glMatrixMode(GL_MODELVIEW);


    /****************************
    / Add your objects from here
    ****************************/
    //add objects

    //glRotatef(-90, 0, 1, 0);

    drawAxes();

    for (int i = 0; i < objects.size(); i++) {
        objects[i]->draw();
    }

    for (int i = 0; i < lights.size(); i++) {
        lights[i].draw();
    }

    glutSwapBuffers();
}


void animate() {
    //codes for any changes in Models, Camera
    glutPostRedisplay();
}

void init() {
    //codes for initialization
    drawaxes = 1;
    angle = 0;

    //clear the screen
    glClearColor(0, 0, 0, 0);

    /************************
    / set-up projection here
    ************************/
    //load the PROJECTION matrix
    glMatrixMode(GL_PROJECTION);

    //initialize the matrix
    glLoadIdentity();

    gluPerspective(80, 1, 1, 1000.0);
    //give PERSPECTIVE parameters
    //field of view in the Y (vertically)
    //aspect ratio that determines the field of view in the X direction (horizontally)
    //near distance
    //far distance
}

int main(int argc, char **argv) {
    loadData();
    glutInit(&argc, argv);
    glutInitWindowSize(windowHeight, windowWidth);
    glutInitWindowPosition(0, 0);
    glutInitDisplayMode(GLUT_DEPTH | GLUT_DOUBLE | GLUT_RGB);    //Depth, Double buffer, RGB color

    glutCreateWindow("Ray Tracer");

    init();

    glEnable(GL_DEPTH_TEST);    //enable Depth Testing

    glutDisplayFunc(display);    //display callback function
    glutIdleFunc(animate);        //what you want to do in the idle time (when no drawing is occuring)

    glutKeyboardFunc(keyboardListener);
    glutSpecialFunc(specialKeyListener);

    glutMainLoop();        //The main loop of OpenGL
    return 0;
}
