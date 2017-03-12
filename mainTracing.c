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
shape *shapes[2];

#define indexSPHERE1 0
#define indexSPHERE2 1

// variables necessary for tracing with the camera
double target[3];
double targetToScreen;
double screenToCam;
double viewDir[3];
double viewDirSpherical[3];
double screenCenter[3];
double camPos[3];
// lrVec and udVec must be perpendicular to viewDir
double lrVec[3];
double udVec[3];
double lrVecSpherical[3];
double udVecSpherical[3];
double backgroundColor[3] = {0,0,0};
int maxDepth = 40;


// converts vector coords to spherical coords
// stores info in array: array[0] = rho, array[1] = phi, array[2] = theta
void vecToSpherical(double vec[3], double spherical[3]) {
    spherical[0] = vecLength(3, vec);
    spherical[1] = atan2(vec[1], vec[0]);
    spherical[2] = acos(vec[2] / spherical[0]);
}

// converts spherical coords to vector coords
// spherical is rho, phi, theta
void sphericalToVec(double spherical[3], double vec[3]) {
    vec[0] = spherical[0] * sin(spherical[1]) * cos(spherical[2]);
    vec[1] = spherical[0] * sin(spherical[1]) * sin(spherical[2]);
    vec[2] = spherical[0] * cos(spherical[1]);
}

// adds phi and theta to the current phi and theta values in viewDir
// modifies left right vector and up down vectors to account for this as well
void modViewDir(double phi, double theta) {
    viewDirSpherical[1] = viewDirSpherical[1] + phi;
    viewDirSpherical[2] = viewDirSpherical[2] + theta;
    lrVecSpherical[1] = lrVecSpherical[1] + phi;
    lrVecSpherical[2] = lrVecSpherical[2] + theta;
    udVecSpherical[1] = udVecSpherical[1] + phi;
    udVecSpherical[2] = udVecSpherical[2] + theta;
    sphericalToVec(viewDirSpherical, viewDir);
    sphericalToVec(lrVecSpherical, lrVec);
    sphericalToVec(udVecSpherical, udVec);
    vecUnit(3, viewDir, viewDir);
    vecUnit(3, lrVec, lrVec);
    vecUnit(3, udVec, udVec);
}

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
int shootRay(double s[3], double d[3], double rgb[3], int depth)  {
    if(depth >= maxDepth)  {
        return 2;
    }
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
        // if(shapes[minIndex] -> reflectivity > 0)  {
        //     // r = d − 2(d⋅n)n
        //     double twodnn[3];
        //     double r[3];
        //     double reflectedRGB[3] = {0,0,0};
        //     vecScale(3, 2 * vecDot(3, d, minNormal), minNormal, twodnn);
        //     vecSubtract(3, d, twodnn, r);
        //     vecUnit(3, r, r);
        //     double nudgedS[3];
        //     double tinyR[3];
        //     vecScale(3, .000001, r, tinyR);
        //     vecAdd(3, minIntersectLoc, tinyR, nudgedS);
        //     shootRay(nudgedS, r, reflectedRGB, depth + 1);
        //     vecScale(3, shapes[minIndex] -> reflectivity, reflectedRGB, reflectedRGB);
        //     vecScale(3, 1 - shapes[minIndex] -> reflectivity, rgb, rgb);
        //     vecAdd(3, reflectedRGB, rgb, rgb);
        // }
    }
    return minIndex == -1;
}

void launchRays()  {
    double pixPos[3];
    double pixPosFinal[3];
    for (int i = -screenWidth / 2; i < screenWidth / 2; i++) {
        for (int j = -screenHeight / 2; j < screenHeight / 2; j++) {
            getNewPoint(screenCenter, lrVec, i, pixPos);
            getNewPoint(pixPos, udVec, j, pixPosFinal);
            
            double rayDir[3];
            vecSubtract(3, pixPosFinal, camPos, rayDir);
            vecUnit(3, rayDir, rayDir);

            double rgb[3] = {1, 1, 1};

            shootRay(pixPosFinal, rayDir, rgb, 0);
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
    vecToSpherical(viewDir, viewDirSpherical);
    vecToSpherical(lrVec, lrVecSpherical);
    vecToSpherical(udVec, udVecSpherical);
    // set targetPos and distances to screen and camera
    target[0] = targetPos[0]; target[1] = targetPos[1]; target[2] = targetPos[2];
    targetToScreen = targetToScreenDist;
    screenToCam = screenToCamDist;
    //calculating screen center point
    getNewPoint(targetPos, viewDir, targetToScreenDist, screenCenter);
    //calculating cam pos
    getNewPoint(screenCenter, viewDir, screenToCamDist, camPos);
}

// sets up a sphere, given a radius, center, and shapeIndex in our shape array
void sphereSetup(double radius, double center[], int shapeIndex, double color[3], double reflectivity)  {
    shapes[shapeIndex] = sphereMalloc();
    sphereInit(shapes[shapeIndex], center, radius, reflectivity, color);
}

void initializeShapes() {
    double sphere1Radius = 30;
    double sphere1Center[3] = {256, 256, 50};
    double sphere1Color[3] = {.2, .8, .1};
    sphereSetup(sphere1Radius, sphere1Center, indexSPHERE1, sphere1Color, .4);

    // double sphere2Center[3] = {100, 100, 45};
    double sphere2Radius = 50;
    double sphere2Center[3] = {200, 220, 80};
    double sphere2Color[3] = {.5, .8, .8};
    // sphereSetup(sphere2Radius, sphere2Center, indexSPHERE2, sphere2Color, 0);
    sceneInitialize(sphere2Center, sphere1Radius * 2, 500);
}

// Move camera on arrow keys + shift + option.
void handleKeyDown(int key, int shiftIsDown, int controlIsDown,
        int altOptionIsDown, int superCommandIsDown) {
    switch(key) {
        double theta = .1;
        double phi = .1;
        case 257:
        break;
        case 262:  // right
        modViewDir(0, theta);
        break;
        case 263:  // left
        modViewDir(0, -theta);
        break;
        case 264:  // down
        if(shiftIsDown)  {
            targetToScreen--;
        }  else  {
            modViewDir(-phi, 0);
        }
        break;
        case 265:  // up
        if(shiftIsDown)  {
            targetToScreen++;
        }  else  {
            modViewDir(phi, 0);
        }
        break;
        default:
        return;
    }
    launchRays();
    printf("we done shot some rays\n");
}

int main(int argc, const char **argv)  {
    projectionType = PERSPECTIVE;
    pixInitialize(WIDTH, HEIGHT, "ray tracing");
    pixClearRGB(0.0, 0.0, 0.0);
    initializeShapes();
    launchRays();
    printf("we done shot some rays\n");
    pixSetKeyDownHandler(handleKeyDown);
    pixSetKeyRepeatHandler(handleKeyDown);
    pixRun();
    for (int i = 0; i < 2; ++i)  {
        sphereDestroy(shapes[i]);
    }
}