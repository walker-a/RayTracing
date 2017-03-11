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
#define PERSPECTIVE 1

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

// launches a ray, does neccessary reflection.
// writes answer to rgb
// returns 0 on intersection, 1 on no intersection
int shootRay(double s[3], double d[3], double rgb[3])  {
            double candidateLoc[3];
            double minIntersectLoc[3];
            double minIntersectDist = -1;
            double minNormal[3];
            int minIndex = -1;
            int max = -1;

            for (int k = 0; k < numShapes; ++k)  {
                double candidateNormal[3];
                if(shapes[k] -> intersection(shapes[k], s, d, candidateLoc, candidateNormal) == 0)  {
                    double stoCandidateVec[3];
                    vecSubtract(3, candidateLoc, s, stoCandidateVec);
                    double stoCandidateDist = vecLength(3, stoCandidateVec);
                    if(minIntersectDist == -1 || 
                        stoCandidateDist < minIntersectDist)  {
                        vecCopy(3, candidateLoc, minIntersectLoc);
                        vecCopy(3, candidateNormal, minNormal);
                        minIntersectDist = stoCandidateDist;
                    }
                }
            }

            // if we hit something
            if(minIndex != -1)  {
                shapes[minIndex] -> color(shapes[minIndex], minIntersectLoc, rgb);

                // if there's a reflection
                if(shapes[minIndex] -> reflectivity > 0)  {
                    // r = d − 2(d⋅n)n
                    double twodnn[3];
                    double r[3];
                    double reflectedRGB[3];
                    vecSubtract(3, d, twodnn, r);
                    vecUnit(3, r, r);
                    shootRay(minIntersectLoc, r, reflectedRGB);
                    vecScale(3, shapes[minIndex] -> reflectivity, reflectedRGB, reflectedRGB);
                    vecScale(3, 1 - shapes[minIndex] -> reflectivity, rgb, rgb);
                    vecAdd(3, reflectedRGB, rgb, rgb);
                }
                return 0;
            }
            return 1;
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
                double pixelCoord[3] = {i, j, camToScreenDist};
                vecSubtract(3, pixelCoord, s, d);
                vecUnit(3, d, d);
            }

            double rgb[3] = {.1,.1,.1};

            shootRay(s, d, rgb);
            pixSetRGB(getScreenCoordX(i), getScreenCoordY(j), rgb[0], rgb[1], rgb[2]);
        }
    }

}

void sphereSetup(double radius, double center[], int shapeIndex)  {
    shapes[shapeIndex] = sphereMalloc();
    sphereInit(shapes[shapeIndex], center, radius, .5);
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