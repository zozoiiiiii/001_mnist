#include "matrix.h"
#include "utils.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <math.h>

void free_matrix(matrix m)
{
    int i;
    for(i = 0; i < m.rows; ++i) free(m.vals[i]);
    free(m.vals);
}

matrix make_matrix(int rows, int cols)
{
    matrix m;
    m.rows = rows;
    m.cols = cols;
    m.vals = calloc(m.rows, sizeof(double *));
    int i;
    for(i = 0; i < m.rows; ++i) m.vals[i] = calloc(m.cols, sizeof(double));
    return m;
}

matrix hold_out_matrix(matrix *m, int n)
{
    int i;
    matrix h;
    h.rows = n;
    h.cols = m->cols;
    h.vals = calloc(h.rows, sizeof(double *));
    for(i = 0; i < n; ++i){
        int index = rand()%m->rows;
        h.vals[i] = m->vals[index];
        m->vals[index] = m->vals[--(m->rows)];
    }
    return h;
}

double *pop_column(matrix *m, int c)
{
    double *col = calloc(m->rows, sizeof(double));
    int i, j;

    for(i = 0; i < m->rows; ++i)
	{
        col[i] = m->vals[i][c];
        for(j = c; j < m->cols-1; ++j)
		{
            m->vals[i][j] = m->vals[i][j+1];
        }
    }
    --m->cols;
    return col;
}

matrix csv_to_matrix(char *filename)
{
	FILE *fp = fopen(filename, "r");
	if(!fp) file_error(filename);

    matrix m;
    m.cols = -1;

	char *line;

	int n = 0;
	int size = 1024;
	m.vals = calloc(size, sizeof(double*));
	while((line = fgetl(fp)))
	{
        if(m.cols == -1)
			m.cols = count_fields(line);

		if(n == size)
		{
			size *= 2;
			m.vals = realloc(m.vals, size*sizeof(double*));
		}
		m.vals[n] = parse_fields(line, m.cols);
		free(line);
		++n;
	}
	m.vals = realloc(m.vals, n*sizeof(double*));
    m.rows = n;
	return m;
}

void print_matrix(matrix m)
{
    int i, j;
    printf("%d X %d Matrix:\n",m.rows, m.cols);
    printf(" __");
    for(j = 0; j < 16*m.cols-1; ++j) printf(" ");
    printf("__ \n");

    printf("|  ");
    for(j = 0; j < 16*m.cols-1; ++j) printf(" ");
    printf("  |\n");

    for(i = 0; i < m.rows; ++i){
        printf("|  ");
        for(j = 0; j < m.cols; ++j){
            printf("%15.7f ", m.vals[i][j]);
        }
        printf(" |\n");
    }
    printf("|__");
    for(j = 0; j < 16*m.cols-1; ++j) printf(" ");
    printf("__|\n");
}
