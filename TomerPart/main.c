#include <stdio.h>
#include <math.h>
#include <assert.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>

#define MAX(x, y) (((x) > (y)) ? (x) : (y))
#define MIN(x, y) (((x) < (y)) ? (x) : (y))
#define FLT_MAX 3.402823e+38
#define ERR_MSG "An Error Has Occured"

///* init matrix of zeros */
double **init_matrix(int rows, int cols) {
    double **matrix = malloc(rows * sizeof(double *));
    assert(matrix != NULL && ERR_MSG);
    for (int i = 0; i < rows; i++) {
        matrix[i] = calloc(cols, sizeof(double));
        assert(matrix[i] != NULL && ERR_MSG);
    }
    return matrix;
}

void print_matrix(double **matrix, int rows, int cols) {
    int i, j;
    for (i = 0; i < rows; i++) {
        for (j = 0; j < cols; j++) {
            printf("%.4f", matrix[i][j]);
            if (j < cols - 1)
                printf(",");
        }
        printf("\n");
    }
}

double norm (double* point1, double* point2, int dimension){
    int i;
    double result = 0;
    double norm;
    for (i = 0; i < dimension; i++) {
        result += (double) pow((point1[i] - point2[i]), 2);
    }
    norm = (double)sqrt(result);
    return norm;
}

double** wam(double** datapoints, int pointnumber, int dimension) {
    //build matrix
    double** wamMatrix;
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
        for (j = 0; j < pointnumber; j++) {
            if (i == j) {
                wamMatrix[i][j] = 0;
            }
            else {
                pointNorm = norm(datapoints[i], datapoints[j], dimension);
                result = (double) pow(M_E, ((-pointNorm) / 2)); //M_E is the mathematical constant e
                wamMatrix[i][j] = result;
            }
        }
    }
    return wamMatrix;
}

double** eval_ddg(double **weighted_matrix, int n) {
    int i,j;
    double sum;
    double **diagonal_degree_matrix = init_matrix(n,n);
    for (i = 0; i < n; i++) {
        sum = 0;
        for (j = 0; j < n; j++) { // sum row
            sum += weighted_matrix[i][j];
        }
        diagonal_degree_matrix[i][i] = sum;
    }
    return diagonal_degree_matrix;
}

double** ddg(double **datapoints, int pointnumber, int dimension) {
    double **weighted_matrix = wam(datapoints, pointnumber, dimension);
    double **diag = eval_ddg(weighted_matrix, pointnumber);
    return diag;
}

double** multiply_matrices(double **A, int rows_A, int cols_A, double **B, int rows_B, int cols_B) {
    assert(cols_A == rows_B && "invalid multiplication");
    double **res = init_matrix(rows_A, cols_B);
    for (int i = 0; i < rows_A; ++i) {
        for (int j = 0; j < cols_B; ++j) {
            for (int k = 0; k < cols_A; ++k) {
                res[i][j] += A[i][k] * B[k][j];
            }
        }
    }
    return res;
}

double** lnorm(double** datapoints,double pointnumber, int dimension){
    double** wamMatrix = wam(datapoints,pointnumber,dimension);
    double** DDGMatrix = ddg(datapoints,pointnumber,dimension);  //לחשוב איך לבצע את זה בצורה חכמה בהתבסס על חישובים קודמים
    double** tmp_multiply;
    double** lnorm;
    int i,j;
    for (i = 0; i < pointnumber; i++){
        for (j = 0; j < pointnumber; j++){
            if (i == j)
                DDGMatrix[i][j] = (double)(1/(double)sqrt(DDGMatrix[i][j])); //D^(-1/2)
            }
        }
    tmp_multiply = multiply_matrices(DDGMatrix,pointnumber,pointnumber,wamMatrix,pointnumber,pointnumber);
    lnorm = multiply_matrices(tmp_multiply,pointnumber,pointnumber,DDGMatrix,pointnumber,pointnumber); //D^(-1/2) * W * D^(-1/2)
    for (i = 0; i < pointnumber; i++){ //I - D^(-1/2) * W * D^(-1/2)
        for (j = 0; j < pointnumber; j++){
            if (i == j)
                lnorm[i][j] = 1 - lnorm[i][j];
            else
                lnorm[i][j] = 0 - lnorm[i][j];
        }
    }
    return lnorm;
}

int isNumber(char *stringNum){
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
    double **points, **centroids, **wamMatrix;
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

    fp = fopen("C:\\Users\\user\\Desktop\\input_1.txt" , "r");
    assert(fp != NULL && "failed to open file");
    while ((ch = fgetc(fp)) != 10) /*check the dimension of the vectors*/
    {
        if (ch == ',')
            dimension++;
    }
    fp = fopen("C:\\Users\\user\\Desktop\\input_1.txt" , "r");
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
    fp = fopen("C:\\Users\\user\\Desktop\\input_1.txt" , "r");
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
    fp = fopen("C:\\Users\\user\\Desktop\\input_1.txt" , "r");
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


//    for (i = 0; i < pointsNumber; i++)
//    {
//        for (j = 0; j < dimension; j++)
//        {
//            printf("%.4f", points[i][j]);
//            if (j < dimension - 1)
//                printf(",");
//        }
//        printf("\n");
//    }

    double** ddgmatirx;
    double** lnorm1;
    double** test;
    wamMatrix = wam(points,pointsNumber,dimension);
    ddgmatirx = ddg(points,pointsNumber,dimension);
    lnorm1= lnorm(points,pointsNumber,dimension);
//    double A[2][2] = {{2.0,2.0},{2.0,2.0}};
//    double B[2][2] ={{3.0,1.0},{3.0,1.0}};
//    int c = 2;
//    test = multiply_matrices(&A,2, 2,&B, 2,2);
//    print_matrix(test,2,2);

    print_matrix(lnorm1,pointsNumber,pointsNumber);
//    for (i = 0; i < pointsNumber; i++)
//    {
//        for (j = 0; j < pointsNumber; j++)
//        {
//            printf("%.4f", wamMatrix[i][j]);
//            if (j < pointsNumber - 1)
//                printf(",");
//        }
//        printf("\n");
//    }






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



    return 0;
}
