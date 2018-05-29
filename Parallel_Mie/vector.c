//The code is taken from https://gist.github.com/EmilHernvall/953968/0fef1b1f826a8c3d8cfb74b2915f17d2944ec1d0
//Was modified for processing the lists of double 05/29/2018 M. Solyanik


#include <stdio.h>
#include <stdlib.h>
#include "vector.h"

void vector_init(Vector *vector) {
    // initialize size and capacity
    vector->size = 0;
    vector->capacity = VECTOR_INITIAL_CAPACITY;
    
    // allocate memory for vector->data
    vector->data = malloc(sizeof(double) * vector->capacity);
}

void vector_append(Vector *vector, double value) {
    // make sure there's room to expand into
    vector_double_capacity_if_full(vector);
    
    // append the value and increment vector->size
    vector->data[vector->size++] = value;
}

double vector_get(Vector *vector, int index) {
    if (index >= vector->size || index < 0) {
        printf("Index %d out of bounds for vector of size %d\n", index, vector->size);
        exit(1);
    }
    return vector->data[index];
}

void vector_set(Vector *vector, int index, double value) {
    // zero fill the vector up to the desired index
    while (index >= vector->size) {
        vector_append(vector, 0);
    }
    
    // set the value at the desired index
    vector->data[index] = value;
}

void vector_double_capacity_if_full(Vector *vector) {
    if (vector->size >= vector->capacity) {
        // double vector->capacity and resize the allocated memory accordingly
        vector->capacity *= 2;
        vector->data = realloc(vector->data, sizeof(double) * vector->capacity);
    }
}

void vector_free(Vector *vector) {
    free(vector->data);
}
