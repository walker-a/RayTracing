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

camCamera camera;
double camToScreenDist;

int getScreenCoordX(double coord) {
    double ratio = WIDTH / worldWidth;
    return coord * ratio;
}

int getScreenCoordY(double coord) {
    double ratio = HEIGHT / worldHeight;
    return coord * ratio;
}

void launchRays(camCamera *cam)  {
    printf("cam pos coords0: %f, %f, %f\n", cam->camPos[0], cam->camPos[1], cam->camPos[2]);
    printf("cam pos coords1: %f, %f, %f\n", cam->camPos[0], cam->camPos[1], cam->camPos[2]);
    double isom[4][4];
    mat44InverseIsometry(cam -> rotation, cam -> translation, isom);
    double camEyePos[4];
    double camWorldPos[4] = {cam -> camPos[0], cam -> camPos[1], cam -> camPos[2], 1};
    mat441Multiply(isom, camWorldPos, camEyePos);
    printf("cam eye pos coords: %f, %f, %f\n", camEyePos[0], camEyePos[1], camEyePos[2]);
    
    for (int i = 0; i < numShapes; i++) {
        for (int j = 0; j < 1; j++) {
            double *worldCoords = shapes[i] -> unif;
            worldCoords[3] = 0;
            mat441Multiply(isom, worldCoords, shapes[i] -> unif);
        }
    }
    printf("sphere center coords: %f, %f, %f\n", shapes[0] -> unif[0], shapes[0] -> unif[1], shapes[0] -> unif[2]);
    
    for (int i = -worldWidth / 2; i < worldWidth / 2; ++i)  {
        for (int j = -worldHeight / 2; j < worldHeight / 2; ++j)  {
            double s[3] = {i, j, camEyePos[2] + camToScreenDist};
            double d[3] = {0, 0, 1};
            
            if (projectionType == PERSPECTIVE) {
                vecSubtract(3, s, camEyePos, d);
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
// need to initialize scene first
void initializeCamera(camCamera *cam) {
    camToScreenDist = 500;
    double *target = shapes[indexSPHERE1] -> unif;
    double camPos[3] = {256, 256, -camToScreenDist};
    camLookAtAndFrom(cam, camPos, target);
}

int main(int argc, const char **argv)  {
    projectionType = PERSPECTIVE;
    pixInitialize(WIDTH, HEIGHT, "ray tracing");
    pixClearRGB(0.0, 0.0, 0.0);
    initializeScene();
    initializeCamera(&camera);
    launchRays(&camera);
    pixRun();
}