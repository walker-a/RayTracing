typedef struct light  {
    double color[3];
    double loc[3];
} light;


void lightInit(light *lightInput, double color[3], double loc[3])  {
    vecCopy(3, color, lightInput -> color);
    vecCopy(3, loc, lightInput -> color);
}