/* 
* Made by Jack Wines and Alex Walker
* 01/39/17
* Compile with make fast, make accurate, or make debug
* run with rayTracing
*
* Alternitively, compile with the following
* clang mainTracing.c 000pixel.h 000pixel.o -lglfw -framework OpenGL
*/

#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include "vector.c"
#include "matrix.c"
#include "shapes.c"
#include "000pixel.h"
#include "light.c"

#define WIDTH 512
#define HEIGHT 512

#define ORTHOGRAPHIC 0
#define PERSPECTIVE 1

#define ANTIALIASING_OFF 0
#define ANTIALIASING_ON 1

int projectionType;
int aliasing;
int numPasses;
double screenHeight = 512;
double screenWidth = 512;

// int numShapes = 6;
// shape *shapes[6];
int numShapes = 11;
shape *shapes[11];

int numLights = 1;
light *lights[2];

#define indexSPHERE1 0
#define indexSPHERE2 1
#define indexSPHERE3 2
#define indexSPHERE4 3
#define indexPLANE1 4
#define indexPLANE2 5
#define indexPLANE3 6
#define indexPLANE4 7
#define indexPLANE5 8
#define indexPLANE6 9
#define indexLIGHTSPHERE 10

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
double backgroundColor[3] = {0, 0, 0};
int maxDepth = 4;
double ambientLight = .2;

int once;

// Uses angle axis rotation to rotate the direction the camera is looking at, as well as the vectors
// describing the screen directions
void rotateView(double theta, double axis[3]) {
    double rot[3][3];
    double tempView[3];
    double tempLR[3];
    double tempUD[3];
    vecCopy(3, viewDir, tempView);
    vecCopy(3, lrVec, tempLR);
    vecCopy(3, udVec, tempUD);
    mat33AngleAxisRotation(theta, axis, rot);
    mat331Multiply(rot, tempView, viewDir);
    mat331Multiply(rot, tempLR, lrVec);
    mat331Multiply(rot, tempUD, udVec);
    vecUnit(3, viewDir, viewDir);
    vecUnit(3, lrVec, lrVec);
    vecUnit(3, udVec, udVec);
}

// converts world x coordinates into screen coordinates from 0 to WIDTH - 1
int getScreenCoordX(double coord) {
    double ratio = WIDTH / screenWidth;
    coord = round(coord * ratio);
    return coord + WIDTH / 2;
}

// converts world y coordinates into screen coordinates from 0 to HEIGHT - 1
int getScreenCoordY(double coord) {
    double ratio = HEIGHT / screenHeight;
    coord = round(coord * ratio);
    return coord + HEIGHT / 2;
}

// given a starting point, a direction, and a length to go in that direction, find finish point, store in finishPoint
// dir should be a unit vector
void getNewPoint(double startPoint[3], double dir[3], double dist, double finishPoint[3]) {
    double vecScalar;
    if (dist < 0) {
        vecScalar = -sqrt(dist * dist / (dir[0] * dir[0] + dir[1] * dir[1] + dir[2] * dir[2]));
    } else {
        vecScalar = sqrt(dist * dist / (dir[0] * dir[0] + dir[1] * dir[1] + dir[2] * dir[2]));
    }
    double largerVec[3];
    vecScale(3, vecScalar, dir, largerVec);
    vecAdd(3, startPoint, largerVec, finishPoint);
}

int contains(shape *itemToLookFor, shape **l, int lLen){
    for (int i = 0; i < lLen; ++i)  {
        if(itemToLookFor == l[i])  {
            return 1;
        }
    }
    return 0;
}

int rayIntersectAvoid(double s[3], double d[3], double intersection[3], double normal[3], shape **avoid, int avoidNum)  {
    int minIndex = -1;
    int max = -1;
    double minIntersectDist = -1;

    for (int k = 0; k < numShapes; ++k)  {
        if(contains(shapes[k], avoid, avoidNum))
            continue;

        double candidateNormal[3];
        double candidateLoc[3];
        if(shapes[k] -> intersection(shapes[k], s, d, candidateLoc, candidateNormal) == 0)  {
            double stoCandidateVec[3];
            vecSubtract(3, candidateLoc, s, stoCandidateVec);
            double stoCandidateDist = vecLength(3, stoCandidateVec);
            if(minIntersectDist == -1 || 
                stoCandidateDist < minIntersectDist)  {
                vecCopy(3, candidateLoc, intersection);
                vecCopy(3, candidateNormal, normal);
                minIntersectDist = stoCandidateDist;
                minIndex = k;
            }
        }
    }
    return minIndex;
}

// returns the index in shapes of the first object the ray hits
// returns -1 otherwise
int rayIntersect(double s[3], double d[3], double intersection[3], double normal[3])  {
    return rayIntersectAvoid(s, d, intersection, normal, NULL, 0);
}

// nudges orig in direction, writes to result
void vecNudge(double orig[3], double direction[3], double result[3])  {
    double tinyDir[3];
    vecScale(3, .000001, direction, tinyDir);
    vecAdd(3, orig, tinyDir, result);
}

// same as haskell zipWith (*)
void vecZipWithMultiply(int length, double *a, double *b, double *result)  {
    for (int i = 0; i < length; ++i)  {
        result[i] = a[i] * b[i];
    }
}

int shootRay(double s[3], double d[3], double rgbFinal[3], int depth);

// handles the reflection
// r = d − 2(d⋅n)n
void reflection(shape *contact, double s[3], double d[3], double normal[3], double out[3], int depth)  {
    double twodnn[3];
    double r[3];
    double reflectedRGB[3] = {0,0,0};
    vecScale(3, 2 * vecDot(3, d, normal), normal, twodnn);
    vecSubtract(3, d, twodnn, r);
    vecUnit(3, r, r);
    double nudgedS[3];
    vecNudge(s, r, nudgedS);
    if (shootRay(nudgedS, r, reflectedRGB, depth + 1) > 0) {
        if(reflectedRGB[0] > 0)  {
            // printf("%s\n", );
        }
        vecScale(3, contact -> reflectivity, reflectedRGB, reflectedRGB);
        double scaledRGB[3];
        vecAdd(3, reflectedRGB, out, out);
    }
}

int lighting(shape *contact, double s[3], double intersectLoc[3], double normal[3], double surfaceCol[3], double rgb[3])  {
    int numUsedLights = 0;
    int intersect = -1;
    for (int i = 0; i < numLights; ++i)  {
        double dirFromContactToLight[3];
        double rayDir[3];
        double nudgedIntersect[3];
        vecSubtract(3, lights[i] -> loc, intersectLoc, dirFromContactToLight);
        vecUnit(3, dirFromContactToLight, rayDir);
        vecNudge(intersectLoc, rayDir, nudgedIntersect);
        double _[3];
        double shadowIntersect[3];

        // if nothing's blocking our shadow ray
        intersect = rayIntersectAvoid(nudgedIntersect, rayDir, shadowIntersect, _, lights[i] -> avoid, lights[i] -> avoidNum);
        double rayFromIntersectToLight[3];
        vecSubtract(3, lights[i] -> loc, shadowIntersect, rayFromIntersectToLight);
        vecUnit(3, rayFromIntersectToLight, rayFromIntersectToLight);

        if(intersect < 0 || vecDot(3, rayFromIntersectToLight, rayDir) < 0.1 )  {
            double dirToLight[3];
            numUsedLights++;
            vecSubtract(3, lights[i] -> loc, intersectLoc, dirToLight);
            vecUnit(3, dirToLight, dirToLight);
            double diffInt = vecDot(3, dirToLight, normal);
            diffInt = diffInt < 0? 0: diffInt;
            double lightColor[3];
            vecZipWithMultiply(3, lights[i] -> color, surfaceCol, lightColor);
            vecScale(3, diffInt, lightColor, lightColor);
            vecAdd(3, lightColor, rgb, rgb);
        }
    }
    return intersect;
}

// launches a ray, does neccessary reflection.
// writes answer to rgb
// returns 0 on intersection, 1 on no intersection
int shootRay(double s[3], double d[3], double rgbFinal[3], int depth)  {

    double black[3] = {0, 0, 0};
    vecCopy(3, black, rgbFinal);

    if(depth >= maxDepth)
        return -2;

    double intersectLoc[3];
    double normal[3];

    int minIndex = rayIntersect(s, d, intersectLoc, normal);

    // exit if we didn't hit anything
    if(minIndex == -1)  {
        vecCopy(3, backgroundColor, rgbFinal);
        return 0;
    }

    shape *contact = shapes[minIndex];
    double rgb[3] = {0,0,0};
    contact -> color(contact, intersectLoc, rgb);
    double surfaceCol[3];
    vecCopy(3, surfaceCol, rgb);

    // ambient lighting calculations
    // double rgbAmbient[3];
    // vecScale(3, ambientLight, rgb, rgbAmbient);
    // vecAdd(3, rgbAmbient, rgbFinal, rgbFinal);

    double reflectionRGB[3] = {0, 0, 0};

    // reflection calculations
    if(contact -> reflectivity > 0)  {
            reflection(contact, intersectLoc, d, normal, reflectionRGB, depth);
            vecZipWithMultiply(3, rgb, reflectionRGB, reflectionRGB);
            vecAdd(3, reflectionRGB, rgbFinal, rgbFinal);
    }

    double lightingRGB[3] = {0,0,0};
    lighting(contact, s, intersectLoc, normal, rgb, lightingRGB);
    vecAdd(3, lightingRGB, rgbFinal, rgbFinal);
    // if(once && vecDot(3, lightingRGB, lightingRGB) != 0)  {
    //     once = 0;
    //     vecPrint(3, lightingRGB);
    //     printf("uhhh\n");
    // }
    if(vecLength(3, rgbFinal) == 0)  {
        // printf("uhhh\n");
    }
    return minIndex;
}

// sets up the camera in relation to its target and the plane representing the screen
void camInitialize(double targetPos[3], double targetToScreenDist, double screenToCamDist) {
    // set targetPos and distances to screen and camera
    vecCopy(3, targetPos, target);
    targetToScreen = targetToScreenDist;
    screenToCam = screenToCamDist;
    //calculating screen center point
    getNewPoint(targetPos, viewDir, targetToScreenDist, screenCenter);
    //calculating cam pos
    getNewPoint(screenCenter, viewDir, screenToCamDist, camPos);
}

// launches all the rays, one for each pixel
void launchRays()  {
    camInitialize(target, targetToScreen, screenToCam);
    double pixPos[3];
    double pixPosFinal[3];
    for (double i = -screenWidth / 2; i < screenWidth / 2; i = i + screenWidth / WIDTH) {
        for (double j = -screenHeight / 2; j < screenHeight / 2; j = j + screenHeight / HEIGHT) {
            double rgb[3];
            vecCopy(3, backgroundColor, rgb);
            
            if (aliasing == ANTIALIASING_ON) {
                double tempRGB[3];
                int counter = 0;
                for (int k = 0; k < numPasses; k++) {
                    for (int l = 0; l < numPasses; l++) {
                        counter++;
                        vecCopy(3, rgb, tempRGB);
                        getNewPoint(screenCenter, lrVec, i + (screenWidth / WIDTH / numPasses) * (k + 1), pixPos);
                        getNewPoint(pixPos, udVec, j + (screenHeight / WIDTH / numPasses) * (l + 1), pixPosFinal);
            
                        //calculates the ray direction
                        double rayDir[3];
                        if (projectionType == ORTHOGRAPHIC) {
                            double rayDirTemp[3] = {-viewDir[0], -viewDir[1], -viewDir[2]};
                            vecCopy(3, rayDirTemp, rayDir);
                            vecUnit(3, rayDir, rayDir);
                        } else {
                            vecSubtract(3, pixPosFinal, camPos, rayDir);
                            vecUnit(3, rayDir, rayDir);
                        }
                        shootRay(pixPosFinal, rayDir, tempRGB, 0);
                        vecAdd(3, rgb, tempRGB, rgb);
                    }
                }
                pixSetRGB(getScreenCoordX(i), getScreenCoordY(j), rgb[0]/counter, 
                          rgb[1]/counter, rgb[2]/counter);
            } else {
                getNewPoint(screenCenter, lrVec, i, pixPos);
                getNewPoint(pixPos, udVec, j, pixPosFinal);

                //calculates the ray direction
                double rayDir[3];
                if (projectionType == ORTHOGRAPHIC) {
                    double rayDirTemp[3] = {-viewDir[0], -viewDir[1], -viewDir[2]};
                    vecCopy(3, rayDirTemp, rayDir);
                    vecUnit(3, rayDir, rayDir);
                } else {
                    vecSubtract(3, pixPosFinal, camPos, rayDir);
                    vecUnit(3, rayDir, rayDir);
                }
                shootRay(pixPosFinal, rayDir, rgb, 0);
                pixSetRGB(getScreenCoordX(i), getScreenCoordY(j), rgb[0], rgb[1], rgb[2]);
            }
        }
    }
}

// prepares the camera and viewing plane for ray tracing, given a target in the scene
// call after initializeShapes()
void sceneInitialize(double targetPos[3], double targetToScreenDist, double screenToCamDist) {
    // defining vector directions for screen plane
    viewDir[0] = 0; viewDir[1] = 0; viewDir[2] = -1;
    lrVec[0] = 1; lrVec[1] = 0; lrVec[2] = 0;
    udVec[0] = 0; udVec[1] = 1; udVec[2] = 0;
    vecUnit(3, viewDir, viewDir);
    vecUnit(3, lrVec, lrVec);
    vecUnit(3, udVec, udVec);
    
    camInitialize(targetPos, targetToScreenDist, screenToCamDist);
}

// sets up a sphere, given a radius, center, and shapeIndex in our shape array
void sphereSetup(double radius, double center[], int shapeIndex, double color[3], double reflectivity, double ambient)  {
    shapes[shapeIndex] = sphereMalloc();
    sphereInit(shapes[shapeIndex], center, radius, reflectivity, color, ambient);
}

// sets up a plane, given a radius, center, and shapeIndex in our shape array
void planeSetup(double normal[3], double center[3], int shapeIndex, double color[3], double reflectivity, double ambient) {
    shapes[shapeIndex] = planeMalloc();
    planeInit(shapes[shapeIndex], center, normal, reflectivity, color, ambient);
}

// sets up our shapes, which are currently three circles
void initializeShapes() {    
    double sphere1Radius = 90;
    double sphere1Center[3] = {160, -256 + sphere1Radius, 0};
    double sphere1Color[3] = {.2, .8, .1};
    sphereSetup(sphere1Radius, sphere1Center, indexSPHERE1, sphere1Color, .3, ambientLight);

    double sphere2Radius = 110;
    double sphere2Center[3] = {-100, -256 + sphere2Radius, -100};
    double sphere2Color[3] = {.5, .8, .8};
    sphereSetup(sphere2Radius, sphere2Center, indexSPHERE2, sphere2Color, .3, ambientLight);

    double sphere3Radius = 70;
    double sphere3Center[3] = {0, -256 + sphere3Radius, 185};
    double sphere3Color[3] = {.2, .5, .6};
    sphereSetup(sphere3Radius, sphere3Center, indexSPHERE3, sphere3Color, .3, ambientLight);

    double sphere4Radius = 30;
    double sphere4Center[3] = {0, 50, 0};
    double sphere4Color[3] = {.3, .3, .3};
    sphereSetup(sphere4Radius, sphere4Center, indexSPHERE4, sphere4Color, 1, ambientLight);
    
    double plane1Normal[3] = {0, 1, 0};
    double plane1Center[3] = {0, -256, 0};
    double plane1Color[3] = {.5, .5, .5};
    planeSetup(plane1Normal, plane1Center, indexPLANE1, plane1Color, .5, ambientLight);

    double plane2Normal[3] = {1, 0, 0};
    double plane2Center[3] = {-500, 0, 0};
    double plane2Color[3] = {.5, .5, .5};
    planeSetup(plane2Normal, plane2Center, indexPLANE2, plane2Color, .5, ambientLight);

    double plane3Normal[3] = {0, 0, 1};
    double plane3Center[3] = {0, 0, -500};
    double plane3Color[3] = {.5, .5, .5};
    planeSetup(plane3Normal, plane3Center, indexPLANE3, plane3Color, .5, ambientLight);

    double plane4Normal[3] = {0, 0, -1};
    double plane4Center[3] = {0, 0, 500};
    double plane4Color[3] = {.5, .5, .5};
    planeSetup(plane4Normal, plane4Center, indexPLANE4, plane4Color, .5, ambientLight);

    double plane5Normal[3] = {0, -1, 0};
    double plane5Center[3] = {0, 500, 0};
    double plane5Color[3] = {.5, .5, .5};
    planeSetup(plane5Normal, plane5Center, indexPLANE5, plane5Color, .5, ambientLight);

    double plane6Normal[3] = {-1, 0, 0};
    double plane6Center[3] = {500, 0, 0};
    double plane6Color[3] = {.5, .5, .5};
    planeSetup(plane6Normal, plane6Center, indexPLANE6, plane6Color, .5, ambientLight);

    double lightSphereRadius = 10;
    double lightSphereCenter[3] = {0, -256 + 400, 0};
    double lightSphereColor[3] = {1, 1, 1};
    shapes[indexLIGHTSPHERE] = sphereMalloc();
    sphereInit(shapes[indexLIGHTSPHERE], lightSphereCenter, lightSphereRadius, 0, lightSphereColor, 1);

    double camTarget[3] = {0, 0, 0};
    sceneInitialize(camTarget, 256, 1000);
}

// sets up our lights
void initializeLights()  {
    lights[0] = malloc(sizeof(light));
    double color[3] = {.6,.6,.6};
    double pos[3] =   {0, -256 + 400, 0};
    lightInit(lights[0], color, pos);
    lightSetAvoid(lights[0], shapes + indexLIGHTSPHERE, 1);

    lights[1] = malloc(sizeof(light));
    double color2[3] = {.6,.6,.6};
    double pos2[3] =   {100, 256, 100};
    lightInit(lights[1], color2, pos2);
}

// Move camera on arrow keys, zoom through shift + up down, enter to change projection type.
void handleKeyDown(int key, int shiftIsDown, int controlIsDown,
        int altOptionIsDown, int superCommandIsDown) {
    double adjustTheta = 3.14 / 16;
    double udVector[3] = {0,1,0};
    switch(key) {
        case 257: // ENTER
        if (projectionType == ORTHOGRAPHIC) {
            projectionType = PERSPECTIVE;
        } else {
            projectionType = ORTHOGRAPHIC;
        }
        case 65: // if 'a' is pressed
        if (aliasing == ANTIALIASING_ON) {
            aliasing = ANTIALIASING_OFF;
        } else {
            aliasing = ANTIALIASING_ON;
        }
        break;
        case 262:  // right
        rotateView(-adjustTheta, udVector);
        break;
        case 263:  // left
        rotateView(adjustTheta, udVector);
        break;
        case 264:  // down
        if(shiftIsDown)  {
            screenWidth *= .9;
            screenHeight *= .9;
        }  else  {
            rotateView(-adjustTheta, lrVec);
        }
        break;
        case 265:  // up
        if(shiftIsDown)  {
            screenWidth *= 1.1;
            screenHeight *= 1.1;
        }  else  {
            rotateView(adjustTheta, lrVec);
        }
        break;
        case 49:  // 1
        numPasses = 1;
        break;
        case 50:  // 2
        numPasses = 2;
        break;
        case 51:  // 3
        numPasses = 3;
        break;
        case 52:  // 4
        numPasses = 4;
        break;
        case 53:  // 5
        numPasses = 5;
        break;
        case 54:  // 6
        numPasses = 6;
        break;
        case 55:  // 7
        numPasses = 7;
        break;
        case 56:  // 8
        numPasses = 8;
        break;
        case 57:  // 9
        numPasses = 9;
        break;
        default:
        return;
    }
    pixClearRGB(backgroundColor[0], backgroundColor[1], backgroundColor[2]);
    launchRays();
}

int main(int argc, const char **argv)  {
    once = 0;
    // if(argc > 0)  {
    //     ambientLight = atof(argv[0]);
    //     maxDepth = atoi(argv[1]);
    // }
    projectionType = PERSPECTIVE;
    aliasing = ANTIALIASING_OFF;
    numPasses = 6;
    pixInitialize(WIDTH, HEIGHT, "ray tracing");
    pixClearRGB(0.0, 0.0, 0.0);
    initializeShapes();
    initializeLights();
    launchRays();
    pixSetKeyDownHandler(handleKeyDown);
    pixSetKeyRepeatHandler(handleKeyDown);
    pixRun();
    for (int i = 0; i < numShapes; ++i)  {
        sphereDestroy(shapes[i]);
    }
}