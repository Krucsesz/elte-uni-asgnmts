#include <cmath>
#include <cstdlib>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

#ifdef __APPLE__             // Mac OS X uses different header
#  include <GLUT/glut.h>
#else                        // Unix and Windows
#  include <GL/gl.h>
#  include <GL/glu.h>
#  include <GL/glut.h>
#endif

using namespace std;

#include "vector.hpp"        // vectors with components of type double
#include "odeint.hpp"        // ODE integration routines, Runge-Kutta ...
using namespace cpl;

const double pi = 4 * atan(1.0);

const double g = 9.8;        // acceleration of gravity

// Paraméterek a kettős ingához
double m1 = 1.0, m2 = 1.0;
double L1 = 1.0, L2 = 0.5;

Vector f(const Vector& x) {  // extended derivative vector
    double t = x[0];
    double th1 = x[1], th2 = x[2];
    double w1  = x[3], w2  = x[4];
    Vector f(5);             // Vector with *5* components
    f[0] = 1;
    f[1] = w1;
    f[2] = w2;
    
    double delta = th1 - th2;
    double den1 = 2*m1 + m2 - m2*cos(2*delta);

    f[3] = (-g*(2*m1+m2)*sin(th1) - m2*g*sin(th1-2*th2) - 2*sin(delta)*m2*(w2*w2*L2 + w1*w1*L1*cos(delta))) / (L1 * den1);
    f[4] = (2*sin(delta)*(w1*w1*L1*(m1+m2) + g*(m1+m2)*cos(th1) + w2*w2*L2*m2*cos(delta))) / (L2 * den1);

    return f;
}

Vector x(5);
// Itt tároljuk a múltbeli pontokat
vector<double> hist_th1, hist_w1, hist_th2, hist_w2;

void getInput() {
    cout << " Double pendulum demo based on pendulum-gl.cpp\n"
        << " --------------------------------\n"
        << " Initial angles th1 and th2 (rad): ";
    double th1, th2, w1, w2;
    cin >> th1 >> th2;
    cout << " Initial angular velocities w1 and w2 (rad/s): ";
    cin >> w1 >> w2;

    x[0] = 0;
    x[1] = th1;
    x[2] = th2;
    x[3] = w1;
    x[4] = w2;
    hist_th1.clear(); hist_w1.clear(); hist_th2.clear(); hist_w2.clear();
}

double dt = 0.05;
double accuracy = 1e-6;

void step() {
    adaptiveRK4Step(x, dt, accuracy, f);
    
    // Szögek visszatekerése és elmentése a memóriába
    double th1 = x[1]; while(th1 >= pi) th1 -= 2*pi; while(th1 < -pi) th1 += 2*pi;
    double th2 = x[2]; while(th2 >= pi) th2 -= 2*pi; while(th2 < -pi) th2 += 2*pi;
    
    hist_th1.push_back(th1);
    hist_w1.push_back(x[3] / 3.0);
    hist_th2.push_back(th2);
    hist_w2.push_back(x[4] / 3.0);
}

double frames_per_second = 30;   // for animation in real time

void animation_step() {
    double start = x[0];
    clock_t start_time = clock();
    step();
    double tau = 1.0 / frames_per_second;
    while (x[0] - start < tau)
        step();
    while ((double(clock()) - start_time)/CLOCKS_PER_SEC < tau)
        ;
    glutPostRedisplay();
}

void drawText(const string& str, double x, double y) {
    glRasterPos2d(x, y);
    int len = str.find('\0');
    for (int i = 0; i < len; i++)
       glutBitmapCharacter(GLUT_BITMAP_HELVETICA_12, str[i]);
}

void drawVariables() {
    // write t, theta, omega values     
    glColor3ub(0, 0, 255);
    ostringstream os;
    os << "t = " << x[0] << ends;
    drawText(os.str(), -2.3, -2.3); // Módosított pozi
    os.seekp(0);
    os << "theta1 = " << x[1] << ", theta2 = " << x[2] << ends;
    drawText(os.str(), -2.3, 2.3); 
    os.seekp(0);
    os << "omega1 = " << x[3] << ", omega2 = " << x[4] << ends;
    drawText(os.str(), -0.5, 2.3);
}

void displayPendulum() {
    glClear(GL_COLOR_BUFFER_BIT);

    // draw the pendulum rod
    double th1 = x[1], th2 = x[2];
    double xp1 = L1 * sin(th1);
    double yp1 = -L1 * cos(th1);
    double xp2 = xp1 + L2 * sin(th2);
    double yp2 = yp1 - L2 * cos(th2);

    glColor3ub(0, 255, 0);
    glLineWidth(2.0);
    glBegin(GL_LINES);
        glVertex2d(0, 0);
        glVertex2d(xp1, yp1);
        glVertex2d(xp1, yp1);
        glVertex2d(xp2, yp2);
    glEnd();

    // draw the pendulum bob
    glPushMatrix();
    glTranslated(xp1, yp1, 0);
    glColor3ub(255, 0, 0);
    const double r = 0.1;
    glPolygonMode(GL_FRONT, GL_FILL);
    glBegin(GL_POLYGON);
        for (int i = 0; i < 12; i++) glVertex2d(r * cos(2*pi*i/12.0), r * sin(2*pi*i/12.0));
    glEnd();
    glPopMatrix();

    glPushMatrix();
    glTranslated(xp2, yp2, 0);
    glColor3ub(255, 0, 0);
    glBegin(GL_POLYGON);
        for (int i = 0; i < 12; i++) glVertex2d(r * cos(2*pi*i/12.0), r * sin(2*pi*i/12.0));
    glEnd();
    glPopMatrix();

    // write t, theta, and omega
    drawVariables();

    // we are using double buffering - write buffer to screen
    glutSwapBuffers();
}


void displayPhasePlot() {
    glClear(GL_COLOR_BUFFER_BIT); // Most már MINDIG törlünk villogás ellen!
    
    // Tengelyek
    glColor3ub(0, 255, 0);
    glBegin(GL_LINES);
        glVertex2d(0, -2.5); glVertex2d(0, 2.5);
        glVertex2d(-pi, 0); glVertex2d(pi, 0);
    glEnd();

    // Vonalak megrajzolása a memóriából
    if (hist_th1.size() > 1) {
        // Belső inga (Lila)
        glColor3ub(128, 0, 128);
        glBegin(GL_LINE_STRIP);
        for(size_t i = 0; i < hist_th1.size(); i++) {
            // Ha átfordult az inga (nagy ugrás), új vonalat kezdünk
            if (i > 0 && abs(hist_th1[i] - hist_th1[i-1]) > pi) {
                glEnd();
                glBegin(GL_LINE_STRIP);
            }
            glVertex2d(hist_th1[i], hist_w1[i]);
        }
        glEnd();

        // Külső inga (Piros)
        glColor3ub(255, 0, 0);
        glBegin(GL_LINE_STRIP);
        for(size_t i = 0; i < hist_th2.size(); i++) {
            if (i > 0 && abs(hist_th2[i] - hist_th2[i-1]) > pi) {
                glEnd();
                glBegin(GL_LINE_STRIP);
            }
            glVertex2d(hist_th2[i], hist_w2[i]);
        }
        glEnd();
    }

    // Szöveg kiírása
    drawVariables();
    glutSwapBuffers();
}

void reshape(int w, int h) {
    glViewport(0, 0, w, h);
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    double aspectRatio = w/double(h);
    double d = 2.5;
    if (aspectRatio > 1)
        glOrtho(-d*aspectRatio, d*aspectRatio, -d, d, -1.0, 1.0);
    else
        glOrtho(-d, d, -d/aspectRatio, d/aspectRatio, -1.0, 1.0);
    glMatrixMode(GL_MODELVIEW);
}

bool running;                    // is the animation running?
bool phasePlot;                  // switch to phase plot if true

void mouse(int button, int state, int mouse_x, int mouse_y) {
    switch (button) {
      case GLUT_LEFT_BUTTON:
        if (state == GLUT_DOWN) {
            if (running) {
                glutIdleFunc(NULL);
                running = false;
            } else {
                glutIdleFunc(animation_step);
                running = true;
            }
        }
        break;
      case GLUT_RIGHT_BUTTON:
        if (state == GLUT_DOWN) {
            if (phasePlot) {
                glutDisplayFunc(displayPendulum);
                phasePlot = false;
            } else {
                glutDisplayFunc(displayPhasePlot);
                phasePlot = true;
            }
            glutPostRedisplay();
        }
        break;
    }
}

int main(int argc, char *argv[]) {
    glutInit(&argc, argv);
    getInput();
    glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGB);
    glutInitWindowSize(600, 600);
    glutInitWindowPosition(100, 100);
    glutCreateWindow(argv[0]);
    glClearColor(1.0, 1.0, 1.0, 0.0);
    glShadeModel(GL_FLAT);
    glutDisplayFunc(displayPendulum);
    glutReshapeFunc(reshape);
    glutMouseFunc(mouse);
    glutIdleFunc(NULL);
    glutMainLoop();
}

