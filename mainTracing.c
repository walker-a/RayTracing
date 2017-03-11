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

#define WIDTH 512
#define HEIGHT 512

#define ORTHOGRAPHIC 0
#define PERSPECTIVE 1

int projectionType;
int screenHeight = 512;
int screenWidth = 512;

int numShapes = 1;
shape *shapes[1];

#define indexSPHERE1 0

// variables necessary for tracing with the camera
double target[3];
double targetToScreen;
double screenToCam;
double viewDir[3];
double screenCenter[3];
double camPos[3];
// lrVec and udVec must be perpendicular to viewDir
double lrVec[3];
double udVec[3];

int getScreenCoordX(double coord) {
    double ratio = WIDTH / screenWidth;
    coord = coord + screenWidth / 2;
    return coord * ratio;
}

int getScreenCoordY(double coord) {
    double ratio = HEIGHT / screenHeight;
    coord = coord + screenHeight / 2;
    return coord * ratio;
}

// given a starting point, a direction, and a length to go in that direction, find finish point, store in finishPoint
// dir should be a unit vector
void getNewPoint(double startPoint[3], double dir[3], double dist, double finishPoint[3]) {
    double vecScalar = sqrt(dist * dist / (dir[0] * dir[0] + dir[1] * dir[1] + dir[2] * dir[2]));
    double largerVec[3];
    vecScale(3, vecScalar, dir, largerVec);
    vecAdd(3, startPoint, largerVec, finishPoint);
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
                        minIndex = k;
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
    double pixPos[3];
    for (int i = -screenWidth / 2; i < screenWidth / 2; i++) {
        for (int j = -screenHeight / 2; j < screenHeight / 2; j++) {
            getNewPoint(screenCenter, lrVec, i, pixPos);
            getNewPoint(pixPos, udVec, j, pixPos);
            
            double rayDir[3];
            vecSubtract(3, pixPos, camPos, rayDir);
            vecUnit(3, rayDir, rayDir);

            double rgb[3] = {.1,.1,.1};

            shootRay(pixPos, rayDir, rgb);
            pixSetRGB(getScreenCoordX(i), getScreenCoordY(j), rgb[0], rgb[1], rgb[2]);
        }
    }
}

// prepares the camera and viewing plane for ray tracing, given a target in the scene
// call after initializeShapes()
void sceneInitialize(double targetPos[3], double targetToScreenDist, double screenToCamDist) {
    // defining vector directions for screen plane
    viewDir[0] = 0; viewDir[1] = 0; viewDir[2] = -1;
    lrVec[0] = 0; lrVec[1] = 1; lrVec[2] = 0;
    udVec[0] = 1; udVec[1] = 0; udVec[2] = 0;
    vecUnit(3, viewDir, viewDir);
    vecUnit(3, lrVec, lrVec);
    vecUnit(3, udVec, udVec);
    // set targetPos and distances to screen and camera
    target[0] = targetPos[0]; target[1] = targetPos[1]; target[2] = targetPos[2];
    targetToScreen = targetToScreenDist;
    screenToCam = screenToCamDist;
    //calculating screen center point
    printf("targetPos: %f, %f, %f\n", targetPos[0], targetPos[1], targetPos[2]);
    printf("viewDir: %f, %f, %f\n", viewDir[0], viewDir[1], viewDir[2]);
    printf("targetToScreenDist: %f\n", targetToScreenDist);
    getNewPoint(targetPos, viewDir, targetToScreenDist, screenCenter);
    printf("screenCenter: %f, %f, %f\n", screenCenter[0], screenCenter[1], screenCenter[2]);
    //calculating cam pos
    getNewPoint(screenCenter, viewDir, screenToCamDist, camPos);
    printf("camPos: %f, %f, %f\n", camPos[0], camPos[1], camPos[2]);
}

// sets up a sphere, given a radius, center, and shapeIndex in our shape array
void sphereSetup(double radius, double center[], int shapeIndex)  {
    shapes[shapeIndex] = sphereMalloc();
    sphereInit(shapes[shapeIndex], center, radius, 0);
}

void initializeShapes() {
    double sphere1Radius = 30;
    double sphere1Center[3] = {256, 256, 50};
    sphereSetup(sphere1Radius, sphere1Center, indexSPHERE1);
    sceneInitialize(sphere1Center, sphere1Radius * 2, 500);
}

int main(int argc, const char **argv)  {
    projectionType = PERSPECTIVE;
    pixInitialize(WIDTH, HEIGHT, "ray tracing");
    pixClearRGB(0.0, 0.0, 0.0);
    initializeShapes();
    launchRays();
    pixRun();
}