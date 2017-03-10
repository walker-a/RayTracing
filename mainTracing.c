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

shape *shapes[1];
int shapeNum;

void launchRays()  {
    for (int i = 0; i < WIDTH; ++i)  {
        for (int j = 0; j < HEIGHT; ++j)  {
            double s[3] = {i, j, 0};
            double d[3] = {0, 0, 1};
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
                pixSetRGB(i, j, rgb[0], rgb[1], rgb[2]);
            }
        }
    }

}

void sphereSetup(int shapeIndex)  {
    double radius = 30;
    double center[3] = {256, 256, 50};
    shapes[shapeIndex] = sphereMalloc();
    sphereInit(shapes[shapeIndex], center, radius);
}

void test()  {
    double s[3] = {256, 256, 0};
    double d[3] = {0,0, 1};
    double loc[3];
    double rgb[3];
    int j = shapes[0] -> intersection(shapes[0], s, d, loc, rgb);
}

int main(int argc, const char **argv)  {
    pixInitialize(WIDTH, HEIGHT, "ray tracing");
    pixClearRGB(0.0, 0.0, 0.0);
    sphereSetup(0);
    // test();
    launchRays();
    pixRun();
}