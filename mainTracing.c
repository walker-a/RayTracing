/* 
* Made by Jack Wines and Alex Walker
* 01/39/17
* compile with the following
* clang mainTracing.c 000pixel.h 000pixel.o -lglfw -framework OpenGL
*/

#include <math.h>
#include <stdlib.h>
#include "vector.c"
#include "matrix.c"
#include "shapes.c"
#include "000pixel.h"
#include "camera.c"

#define WIDTH 512
#define HEIGHT 512

#define ORTHOGRAPHIC 0
#define PERSPECTIVE 0

int projectionType;
int worldHeight = 512;
int worldWidth = 512;

int numShapes = 1;
shape *shapes[1];

#define indexSPHERE1 0

camCamera cam;
double camToScreenDist;

int getScreenCoordX(double coord) {
    double ratio = WIDTH / worldWidth;
    return coord * ratio;
}

int getScreenCoordY(double coord) {
    double ratio = HEIGHT / worldHeight;
    return coord * ratio;
}

void launchRays()  {
    double isom[4][4];
    mat44InverseIsometry(cam.rotation, cam.translation, isom);
    
    for (int i = 0; i < numShapes; i++) {
        for (int j = 0; j < 1; j++) {
            double *worldCoords = shapes[i] -> unif;
            worldCoords[3] = 1;
            mat441Multiply(isom, worldCoords, shapes[i] -> unif);
        }
    }
    
    for (int i = 0; i < worldWidth; ++i)  {
        for (int j = 0; j < worldHeight; ++j)  {
            double *sTemp = cam.camPos;
            double s[3] = {s[0], s[1], s[2]};
            double d[3] = {0, 0, 1};
            
            if (projectionType == PERSPECTIVE) {
                d[0] = i;
                d[1] = j;
                d[2] = camToScreenDist;
                vecUnit(3, d, d);
            }
            
            double rgb[3] = {0, 0, 0};
            double loc[3] = {-1, -1, -1};
            int maxIndex = -1;
            int max = -1;

            // for (int k = 0; k < shapeNum; ++k)  {
            //     if(shapes[k] -> intersection(shapes[k], s, d, loc, rgb) != 0)  {
            //         if(loc[2] > max)  {
            //             max = loc[2];
            //             maxIndex = k;
            //         }
            //     }
            // }
            if(shapes[0] -> intersection(shapes[0], s, d, loc, rgb) == 0)  {
                printf("got here\n");
                int iScreen = getScreenCoordX(i);
                int jScreen = getScreenCoordY(j);
                pixSetRGB(iScreen, jScreen, rgb[0], rgb[1], rgb[2]);
            }
        }
    }

}

void sphereSetup(double radius, double center[], int shapeIndex)  {
    shapes[shapeIndex] = sphereMalloc();
    sphereInit(shapes[shapeIndex], center, radius);
}

void initializeScene() {
    double sphere1radius = 30;
    double sphere1Center[3] = {256, 256, 50};
    sphereSetup(sphere1radius, sphere1Center, indexSPHERE1);
}


// probably want to abstract the camera later
void initializeCamera() {
    double *target = shapes[indexSPHERE1] -> unif;
    double camPos[3] = {256, 256, -camToScreenDist};
    camLookAtAndFrom(&cam, camPos, target); 
}

int main(int argc, const char **argv)  {
    projectionType = PERSPECTIVE;
    pixInitialize(WIDTH, HEIGHT, "ray tracing");
    pixClearRGB(0.0, 0.0, 0.0);
    initializeScene();
    initializeCamera();
    camToScreenDist = 1000;
    launchRays();
    pixRun();
}