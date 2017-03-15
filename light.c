typedef struct light  {
    double color[3];
    double loc[3];
    shape **avoid;
    int avoidNum;
} light;

void lightSetAvoid(light *lightInput, shape **avoidList, int avoidListNum)  {
    lightInput -> avoid = avoidList;
    lightInput -> avoidNum = avoidListNum;
}

void lightInit(light *lightInput, double color[3], double loc[3])  {
    vecCopy(3, color, lightInput -> color);
    vecCopy(3, loc, lightInput -> loc);
    lightInput -> avoid = NULL;
    lightInput -> avoidNum = 0;
}