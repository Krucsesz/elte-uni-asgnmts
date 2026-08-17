#include <GL/glut.h>
#include <vector>
#include <cstdlib>
#include <ctime>
#include <cmath>
#include <iostream>

using namespace std;

enum BoundaryMode { OPEN, PERIODIC, ALIVE, RANDOM_STATIC };

int width = 100;          // rácssuélesség (cella)
int height = 100;         // rácsmagasság (cella)
int n_rule = 1;
BoundaryMode boundary_mode = PERIODIC;    // aktív peremfeltétel
double initial_density = 0.3;          // kezdeti sejtsűrűség

vector<vector<int>> grid;

// sejtinfók
int get_cell(int x, int y) {
    if (x >= 0 && x < width && y >= 0 && y < height) {
        return grid[y][x];
    }
    
    switch (boundary_mode) {
        case PERIODIC:
            return grid[(y + height) % height][(x + width) % width];
        case OPEN:
            return 0;
        case ALIVE:
            return 1;
        case RANDOM_STATIC:
            int hash = (x * 73856093 ^ y * 19349663);
            return (abs(hash) % 2);
    }
    return 0;
}

void init_random() {
    grid.assign(height, vector<int>(width, 0));
    srand(time(nullptr));
    for (int y = 0; y < height; ++y) {
        for (int x = 0; x < width; ++x) {
            if ((rand() / (double)RAND_MAX) < initial_density) {
                grid[y][x] = 1;
            }
        }
    }
}

void step() {
    vector<vector<int>> next_grid = grid;
    for (int y = 0; y < height; ++y) {
        for (int x = 0; x < width; ++x) {
            int live_neighbors = 0;
            
            // 8 szomszéd vizsgálata
            for (int dy = -1; dy <= 1; ++dy) {
                for (int dx = -1; dx <= 1; ++dx) {
                    if (dx == 0 && dy == 0) continue;
                    live_neighbors += get_cell(x + dx, y + dy);
                }
            }
            
            if (live_neighbors == n_rule) {
                next_grid[y][x] = grid[y][x]; // identitás
            } else if (live_neighbors == n_rule + 1) {
                next_grid[y][x] = 1;          // születés
            } else {
                next_grid[y][x] = 0;          // pusztulás
            }
        }
    }
    grid = next_grid;
}

void display() {
    glClear(GL_COLOR_BUFFER_BIT);

    glColor3f(0.0, 0.8, 0.2); 
    
    glBegin(GL_QUADS);
    for (int y = 0; y < height; ++y) {
        for (int x = 0; x < width; ++x) {
            if (grid[y][x]) {
                glVertex2i(x, y);
                glVertex2i(x + 1, y);
                glVertex2i(x + 1, y + 1);
                glVertex2i(x, y + 1);
            }
        }
    }
    glEnd();
    
    glutSwapBuffers();
}

void timer(int) {
    step();
    glutPostRedisplay();
    glutTimerFunc(100, timer, 0);
}

void reshape(int w, int h) {
    glViewport(0, 0, w, h);
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    gluOrtho2D(0, width, height, 0); 
    glMatrixMode(GL_MODELVIEW);
}

int main(int argc, char** argv) {
    init_random();

    glutInit(&argc, argv);
    glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGB);
    
    glutInitWindowSize(800, 800); 
    glutCreateWindow("Cellular Automaton");

    glClearColor(0.1, 0.1, 0.1, 1.0);

    glutDisplayFunc(display);
    glutReshapeFunc(reshape);
    glutTimerFunc(100, timer, 0);

    glutMainLoop();
    return 0;
}