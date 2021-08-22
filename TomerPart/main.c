#include <stdio.h>
#include <math.h>
#include <assert.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>

#define MAX(x, y) (((x) > (y)) ? (x) : (y))
#define MIN(x, y) (((x) < (y)) ? (x) : (y))
#define FLT_MAX 3.402823e+38


double norm (double* point1, double* point2, int dimension){
    int i;
    double result = 0;
    double norm;
    for (i = 0; i < dimension; i++){
        result += (double)pow((point1[i]-point2[i]),2);
    norm = (double)sqrt(result);
    return norm;
    }
}

double** wam(double** datapoints, double pointnumber, int dimension) {
    //build matrix
    double ** wamMatrix;
    int i, j;
    double pointNorm;
    double result;
    wamMatrix = malloc(pointnumber * sizeof(double *));
    assert(wamMatrix != NULL && "malloc failed");
    for (i = 0; i < pointnumber; i++) {
        wamMatrix[i] = malloc(pointnumber * sizeof(double));
        assert(wamMatrix[i] != NULL && "malloc failed");
    }
    for (i = 0; i < pointnumber; i++) {
        for (j = 0; i < pointnumber; i++) {
            if (i == j)
                wamMatrix[i][j] = 0;
            else {
                pointNorm = norm(datapoints[i], datapoints[j], dimension);
                result = (double) pow(M_E, ((-pointNorm) / 2)); //M_E is the mathematical constant e
                if (result > 0)
                    wamMatrix[i][j] = result;
                else
                    wamMatrix[i][j] = 0;
            }
        }
    }
    return wamMatrix;
}

double** lnorm(double** datapoints,double pointnumber, int dimension){
    double** wamMatrix = wam(datapoints,pointnumber,dimension);
    double** DDGMatrix; //לחשוב איך לבצעאת זה בצורה חכמה בהתבסס על חישובים קודמים
    double** IMatrix;
    double** lnorm;
    int i,j;
    for (i = 0; i < pointnumber; i++){
        for (j = 0; i < pointnumber; i++){
            if (i == j)
                DDGMatrix[i][j] = (double)pow(DDGMatrix[i][j],-1/2);
            }
        }
    //lnorm = Multiplicationmatrices(DDGMatrix,wamMatrix); //להשלים פונקציה מאדם
    //lnorm = Multiplicationmatrices(lnorm,DDGMatrix); //D^(-1/2) * W * D^(-1/2)
    for (i = 0; i < pointnumber; i++){ //I - D^(-1/2) * W * D^(-1/2)
        for (j = 0; i < pointnumber; i++){
            if (i == j)
                lnorm[i][j] = 1 - lnorm[i][j];
            else
                lnorm[i][j] = 0 - lnorm[i][j];
        }
    }
    return lnorm;
}

int isNumber(char *stringNum)
{
    int length, i;
    length = strlen(stringNum);
    for (i = 0; i < length; i++)
        if (!isdigit(stringNum[i]))
        {
            return 0;
        }
    return 1;
}

int main(int argc, char **argv) {
    double **points, **centroids, **wam;
    int *clusters;
    int max_itter = 200;
    int dimension = 1;
    int pointsNumber = 0;
    int maxlinelen = 0;
    int currlinelen = 0;
    int k, i, j, ch, ClusterNumber;
    FILE *fp;
    int changed = 1;
    char *cordinate;
    char *line;

    assert(argc == 4 && "invalid number of args");
    assert(isNumber(argv[1]) && "1st arg is not a number");
    k = atoi(argv[1]);

    fp = fopen("input1.txt" , "r");
    assert(fp != NULL && "failed to open file");
    while ((ch = fgetc(fp)) != 10) /*check the dimension of the vectors*/
    {
        if (ch == ',')
            dimension++;
        printf("%d",dimension);
    }
    fp = fopen( "input1.txt" , "r");;
    while ((ch = fgetc(fp)) != EOF) /*check the number of the vectors*/
    {
        currlinelen++;
        if (ch == 10)
        {
            pointsNumber++;
            maxlinelen = MAX(maxlinelen, currlinelen);
            currlinelen = 0;
        }
    }

    /*malloc points matrix and read points */
    points = malloc(pointsNumber * sizeof(double *));
    assert(points != NULL && "malloc failed");

    for (i = 0; i < pointsNumber; i++)
    {
        points[i] = malloc(dimension * sizeof(double));
        assert(points[i] != NULL && "malloc failed");
    }
    centroids = malloc(k * sizeof(double *));
    assert(centroids != NULL && "malloc failed");
    for (i = 0; i < k; i++)
    {
        centroids[i] = malloc(dimension * sizeof(double));
        assert(centroids[i] != NULL && "malloc failed");
    }
    fp = fopen( "input1.txt" , "r");;
    i = 0;
    line = malloc(maxlinelen * sizeof(char));
    while (fgets(line, maxlinelen + 1, fp) != NULL)
    {
        cordinate = strtok(line, ",");
        for (j = 0; j < dimension; j++)
        {
            points[i][j] = atof(cordinate);
            if (i < k)
                centroids[i][j] = atof(cordinate);
            cordinate = strtok(NULL, ",");
        }
        i++;
    }
    fp = fopen( "input1.txt" , "r");;
    /*INIT CLUSTERS*/
    assert(NULL != (clusters = calloc(pointsNumber, sizeof(int))) && "calloc failed");
    for (i = 0; i < pointsNumber; i++)
    {
        if (i < k)
            clusters[i] = i;
        else
            clusters[i] = -1;
    }
    fclose(fp);


    for (i = 0; i < k; i++)
    {
        for (j = 0; j < dimension; j++)
        {
            printf("%.4f", points[i][j]);
            if (j < dimension - 1)
                printf(",");
        }
        printf("\n");
    }


    wam = (points,pointsNumber,dimension);
    for (i = 0; i < pointsNumber; i++)
    {
        for (j = 0; j < pointsNumber; j++)
        {
            printf("%.4f", wam[i][j]);
            if (j < pointsNumber - 1)
                printf(",");
        }
        printf("\n");
    }






    //האלגוריתם עצמו
//    while (max_itter > 0 && changed)
//    {
//        changed = 0;
//        max_itter--;
//        for (i = 0; i < pointsNumber; i++)
//        {
//            ClusterNumber = findMinCent(points[i], centroids, k, dimension);
//            clusters[i] = ClusterNumber;
//        }
//        changed = UpdateAllAvg(centroids, clusters, points, k, pointsNumber, dimension);
//    }

    /*FINAL PRINT*/
//    for (i = 0; i < k; i++)
//    {
//        for (j = 0; j < dimension; j++)
//        {
//            printf("%.4f", centroids[i][j]);
//            if (j < dimension - 1)
//                printf(",");
//        }
//        printf("\n");
//    }

    for (i = 0; i < k; i++)
        free(centroids[i]);
    free(centroids);
    for (i = 0; i < pointsNumber; i++)
        free(points[i]);
    free(points);
    free(clusters);
    free(line);



    printf("Hello, World!\n");
    return 0;
}
