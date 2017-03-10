typedef struct shape{
    double* unif;
    int unifDim;
    int (*intersection)(struct shape*, double[3], double[3], double[3], double[3]);
}shape;

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

// returns whether intersection occurs, and if so, writes to a location
// retuns 0 on intersection, 1 on not
// d is a unit direction vector
// s is where the ray started
// taken from https://en.wikipedia.org/wiki/Ray_tracing_(graphics)#Example
int sphereIntersection(shape *inputSphere, double s[3], double d[3], double intersectionLoc[3], double rgb[3])  {
    double radius = inputSphere -> unif[3];
    double center* = inputSphere -> unif;
    double[2] quadResults;
    double v[3];
    double vdotd = vecDot(v, d);
    double sqrtResults = vdotd * vdotd - (vecDot(v, v) - radius^2);
    if(sqrtResults < 0)  {
        return 1;
    }
    sqrtResults = sqrt(sqrtResults);
    quadResults[0] = -vdotd + sqrtResults;
    quadResults[1] = -vdotd - sqrtResults;
    t = quadResults[1];
    if(t < 0)  {
        return 0;
    }
    double dtScaled[3];
    vecScale(3, t, d, dtScaled);
    vecAdd(3, s, dtScaled, intersectionLoc);
    double sphereColor[3] = {.5, .7, .3};
    vecCopy(sphereColor, rgb);
}

// mallocs the space for a sphere and its innards.
shape *sphereMalloc()  {
     shape *toReturnShape = malloc(sizeof(shape));
     toReturnShape -> unif = malloc(sizeof(double)*4);
}

// sets the shape's unifs. [x, y, z, r]
void sphereInit(shape *toReturnShape, double center[3], double radius)  {
    toReturnShape -> unifDim = 4;
    toReturnShape -> intersection = sphereIntersection;
    vecCopy(3, )
    toReturnShape -> unif[3] = radius;
}