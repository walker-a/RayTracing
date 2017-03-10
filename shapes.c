typedef struct shape{
    double* unif;
    int unifDim;
    int (*intersection)(struct shape*, double[3], double[3], double[3]);
}shape;

// returns whether intersection occurs, and if so, writes to a location
int sphereIntersection(shape *inputSphere, double start[3], double directionVector[3], double intersectionLocation[3])  {

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

   for (int i = 0; i < 3; ++i)  {
       toReturnShape -> unif[i] = center;
   }
   
   toReturnShape -> unif[3] = radius;
}


int main() {
    return 0;
}