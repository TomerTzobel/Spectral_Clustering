#define ERR_MSG "An Error Has Occured"
#include <stdio.h>
#include <math.h>
#include <assert.h>
#include <stdlib.h>
#include <string.h>
#include "nsc.h"

int main(int argc, char **argv) {
    int k;
    char *filename, *goal;
    assert(argc == 4 && ERR_MSG);
    k = atoi(argv[1]);
    goal = argv[2];
    filename = argv[3];
    nsc(k, goal, filename, 0);
    return 0;
}
