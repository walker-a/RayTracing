/*
* Made by Jack Wines and Alex Walker for their final project for graphics.
* Compile with make
* run with rayTracing
*
* Alternatively, complile with:
* clang -O3 mainTracing.c 000pixel.h 000pixel.o -lglfw -framework OpenGL
*/

// shape struct
// unif keeps track of coordinate points defining the object as well as color
// each shape has its own intersection and find color functions
// stores reflectiviity, a value for ambient light (used in lighting calculations) and refraction info
typedef struct shape  {
    double* unif;
    int unifDim;
    int (*intersection)(struct shape*, double[3], double[3], double[3], double[3]);
    void (*color)(struct shape*, double[3], double[3]);
    // fraction of light that should be reflected
    double reflectivity;
    double ambientLight;
    double refrIndex;
    double transparency;
} shape;

// returns 0 on success (real results), 1 on faliure (comples results)
int quadraticFormula (double results[2], double a, double b, double c)  {
    if(a * c > b*b)  {
        return 1;
    }
    double sqrtAC = sqrt(b*b - 4*a*c);
    results[0] = (-b + sqrtAC) / (2*a);
    results[1] = (-b - sqrtAC) / (2*a);
    return 0;
}

// stores the sphere's color into rgb
void sphereColor(shape *inputSphere, double intersectLoc[3], double rgb[3])  {
    vecCopy(3, inputSphere -> unif + 4, rgb);
}

// returns whether intersection occurs, and if so, writes to a location
// retuns 0 on intersection, 1 on not
// d is a unit direction vector
// s is where the ray started
// taken from https://en.wikipedia.org/wiki/Ray_tracing_(graphics)#Example
int sphereIntersect(shape *sphere, double s[3], double d[3], double intersectLoc[3], double normal[3])  {
    double radius = sphere -> unif[3];
    double *center = sphere -> unif;
    double quadResults[2];
    double v[3];
    vecSubtract(3, s, center, v);
    double vdotd = vecDot(3, v, d);
    double sqrtResults = vdotd * vdotd - (vecDot(3, v, v) - radius * radius);
    if(sqrtResults < 0)  {
        return 1;
    }
    sqrtResults = sqrt(sqrtResults);
    quadResults[0] = -vdotd + sqrtResults;
    quadResults[1] = -vdotd - sqrtResults;
    double t = quadResults[1] >= 0? quadResults[1]: quadResults[0];
    if(t < 0)  {
        return 2;
    }
    double dtScaled[3];
    vecScale(3, t, d, dtScaled);
    vecAdd(3, s, dtScaled, intersectLoc);

    double notUnitNormal[3];
    vecSubtract(3, intersectLoc, center, notUnitNormal);
    vecUnit(3, notUnitNormal, normal);
    return 0;
}


// formula used to find intersection:
// d = ((p0 - l0) . n)/(l . n)
//      where 
//          p0 = point on plane
//          n   = normal
//          l0 = point on line
//          l   = line vector
int planeIntersect(shape *plane, double l0[3], double l[3], double intersectLoc[3], double normal[3])  {
    double p0Minusl0[3];
    double *p0 = plane -> unif;
    double *n  = plane -> unif + 3;
    vecSubtract(3, p0, l0, p0Minusl0);
    double d = vecDot(3, p0Minusl0, n) / vecDot(3, l, n);
    if(d < 0)
        return 1;
    vecCopy(3, n, normal);
    double froml0ToPlane[3];
    vecScale(3, d, l, froml0ToPlane);
    vecAdd(3, l0, froml0ToPlane, intersectLoc);
    return 0;
}

// stores the planes color in rgb
void planeColor(shape *plane, double intersectLoc[3], double rgb[3])  {
    vecCopy(3, plane -> unif + 6, rgb);
}

// mallocs the space for a plane and its innards.
shape *planeMalloc()  {
    shape *toReturnShape = malloc(sizeof(shape));
    toReturnShape -> unif = malloc(sizeof(double) * 9);
    return toReturnShape;
}

// initializes a plane
void planeInit(shape *plane, double point[3], double normal[3], double reflectivity, double rgb[3], double ambient,
               double refractionIndex, double transparence)  {
    plane -> ambientLight = ambient;
    vecCopy(3, point,  plane -> unif);
    vecCopy(3, normal, plane -> unif + 3);
    vecCopy(3, rgb,    plane -> unif + 6);
    plane -> intersection = planeIntersect;
    plane -> color = planeColor;
    plane -> reflectivity = reflectivity;
    plane -> refrIndex = refractionIndex;
    plane -> transparency = transparence;
}

// mallocs the space for a sphere and its innards.
shape *sphereMalloc()  {
     shape *toReturnShape = malloc(sizeof(shape));
     toReturnShape -> unif = malloc(sizeof(double)*7);
     return toReturnShape;
}

// destroys all objects
void objectDestroy(shape *toDestroy)  {
    free(toDestroy -> unif);
    free(toDestroy);
}

// sets the shape's unifs. [x, y, z, r]
void sphereInit(shape *toReturnShape, double center[3], double radius, double reflectivity, double colors[3], double ambient,
                double refractionIndex, double transparence)  {
    toReturnShape -> ambientLight = ambient;
    toReturnShape -> unifDim = 4;
    toReturnShape -> intersection = sphereIntersect;
    toReturnShape -> color = sphereColor;
    vecCopy(3, center, toReturnShape -> unif);
    toReturnShape -> unif[3] = radius;
    toReturnShape -> reflectivity = reflectivity;
    vecCopy(3, colors, toReturnShape -> unif + 4);
    toReturnShape -> refrIndex = refractionIndex;
    toReturnShape -> transparency = transparence;
}