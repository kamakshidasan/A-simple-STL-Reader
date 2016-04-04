// Author : Adhitya Kamakshidasan

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <map>
#include <vector>
#include <cmath>
#include <limits>

// Error Codes
#define FAILURE 0
#define SUCCESS 1
#define INCORRECT_FILE_FORMAT 2
#define FILE_OPENING_ERROR 3
#define IMPROPER_SYNTAX 4
#define NULL_POINTER 5

// Different states while parsing
#define UNKNOWN 0
#define SOLID 1
#define FACET 2
#define NORMAL 3
#define OUTER 4
#define LOOP 5
#define VERTEX 6
#define ENDLOOP 7
#define ENDFACET 8
#define BRANCH 9
#define ENDSOLID 10

#define EPSILON 2.2204460492503131e-016

typedef int ERROR_TYPE;

struct TriangleMesh {
    int vertex_count;
    double* vertex_coords;

    int triangle_count;
    int* triangle_vertices;
    double* normal_coords;
};

// Structure for storing a vertex
struct Vertex {
    double x;
    double y;
    double z;
};

ERROR_TYPE init_mesh(TriangleMesh*);
ERROR_TYPE free_mesh(TriangleMesh*);
int compareToken(const char*, char*, int);
int checkFormat(const char*);
ERROR_TYPE test_mesh(TriangleMesh*);
ERROR_TYPE store_STL(TriangleMesh*, int, std::vector<Vertex>& vertex_coords, std::vector<int>& triangle_vertices, std::vector<Vertex>& normal_coords);
ERROR_TYPE read_STL(const char* filename, TriangleMesh* mesh);

// Initialise elements of a mesh
ERROR_TYPE init_mesh(TriangleMesh* mesh)
{
    mesh->vertex_count = 0;
    mesh->vertex_coords = NULL;

    mesh->triangle_count = 0;
    mesh->triangle_vertices = NULL;
    mesh->normal_coords = NULL;

    return 0;
}

// Free the elements of a mesh
ERROR_TYPE free_mesh(TriangleMesh* mesh)
{

    if (mesh->vertex_coords != NULL && mesh->normal_coords != NULL && mesh->triangle_vertices != NULL) {
        mesh->vertex_count = 0;
        mesh->triangle_count = 0;

        free(mesh->vertex_coords);
        mesh->vertex_coords = NULL;

        free(mesh->triangle_vertices);
        mesh->triangle_vertices = NULL;

        free(mesh->normal_coords);
        mesh->normal_coords = NULL;
        return SUCCESS;
    }
    else {
        return NULL_POINTER;
    }
}

// There are various token that are present,
// For example: solid, facet are all tokens
// If the expected token and the input does not match, this function returns an unknown state
int compareToken(const char* token, char* input, int state)
{
    if (strcmp(token, input) == 0) {
        return state;
    }
    else {
        return UNKNOWN;
    }
}

// Find if the file has .stl extension
int checkFormat(const char* filename)
{
    int filename_len = strlen(filename);
    //Take the last four characters of the filename and check!
    if (filename_len >= 4 && strcmp(filename + filename_len - 4, ".stl") == 0) {
        return 1;
    }
    else {
        return 0;
    }
}

// Print the elements of the mesh
ERROR_TYPE test_mesh(TriangleMesh* mesh)
{
    // Check whether none of the elements are null
    if (mesh->vertex_coords != NULL && mesh->normal_coords != NULL && mesh->triangle_vertices != NULL) {
        printf("\nNormal Coordinates:\n");
        for (int i = 0; i < mesh->triangle_count; i++) {
            int j = 3 * i;
            // Print only the first three decimal places
            printf("%.3lf %.3lf %.3lf\n", mesh->normal_coords[j + 0], mesh->normal_coords[j + 1], mesh->normal_coords[j + 2]);
        }

        printf("\nTriangle Vertices:\n");

        for (int i = 0; i < mesh->triangle_count; i++) {
            int j = 3 * i;
            printf("%d %d %d\n", mesh->triangle_vertices[j + 0], mesh->triangle_vertices[j + 1], mesh->triangle_vertices[j + 2]);
        }

        printf("\nVertex Coordinates:\n");

        for (int i = 0; i < mesh->vertex_count; i++) {
            int j = 3 * i;
            // Print only the first decimal place
            printf("%.1lf %.1lf %.1lf\n", mesh->vertex_coords[j + 0], mesh->vertex_coords[j + 1], mesh->vertex_coords[j + 2]);
        }

        printf("\nTriangle Count: %d\n", mesh->triangle_count);
        printf("\nVertex Count: %d\n", mesh->vertex_count);
        return SUCCESS;
    }
    else {
        return NULL_POINTER;
    }
}

// The function read_STL only *reads* and the checks the file for errors
// Values are stored in local variables in read_STL.
// In this function, those local variables are passed and are stored in the TriangleMesh structure
ERROR_TYPE store_STL(TriangleMesh* mesh, int triangle_count, std::vector<Vertex>& vertex_coords,
    std::vector<int>& triangle_vertices, std::vector<Vertex>& normal_coords)
{

    // Allocate exact memory that is needed - The number of normal coordinates = number of triangles
    // In each normal coordinate, there are three more values - Therefore we multiply by 3
    mesh->normal_coords = (double*)malloc(3 * triangle_count * sizeof(double));
    if (mesh->normal_coords != NULL) {
        for (int i = 0; i < normal_coords.size(); i++) {
            int j = 3 * i;
            mesh->normal_coords[j + 0] = normal_coords[i].x;
            mesh->normal_coords[j + 1] = normal_coords[i].y;
            mesh->normal_coords[j + 2] = normal_coords[i].z;
        }
    }
    else {
        return NULL_POINTER;
    }

    mesh->triangle_count = triangle_count;
    mesh->vertex_count = vertex_coords.size();

    mesh->triangle_vertices = (int*)malloc(3 * triangle_count * sizeof(int));
    if (mesh->normal_coords != NULL) {
        for (int i = 0; i < triangle_vertices.size(); i++) {
            mesh->triangle_vertices[i] = triangle_vertices[i];
        }
    }
    else {
        return NULL_POINTER;
    }

    mesh->vertex_coords = (double*)malloc(3 * vertex_coords.size() * sizeof(double));
    if (mesh->vertex_coords != NULL) {
        for (int i = 0; i < vertex_coords.size(); i++) {
            int j = 3 * i;
            mesh->vertex_coords[j + 0] = vertex_coords[i].x;
            mesh->vertex_coords[j + 1] = vertex_coords[i].y;
            mesh->vertex_coords[j + 2] = vertex_coords[i].z;
        }
    }
    else {
        return NULL_POINTER;
    }
    return SUCCESS;
}

int findVertex(std::vector<Vertex>& vertex_coords, Vertex vertex1)
{
    int j;
    double dx, dy, dz, answer;
    for (int i = 0; i < vertex_coords.size(); i++) {
        j = 3 * i;
        Vertex vertex2 = vertex_coords[i];
        dx = (vertex2.x - vertex1.x);
        dy = (vertex2.y - vertex1.y);
        dz = (vertex2.z - vertex1.z);

        answer = sqrt((dx * dx) + (dy * dy) + (dz * dz));
        if (answer >= 0 && answer < EPSILON) {
            return i;
        }
    }
    return -1;
}

// The parsing of .stl file can be thought of as a state diagram
// There are different states that a token can exist in
// Here, each token is seen one by one and the next expecting state is found
// If the expecting state and the token doesn't match, then it exits giving out an error!
ERROR_TYPE read_STL(const char* filename, TriangleMesh* mesh)
{
    if (checkFormat(filename)) {
        FILE* fp;
        char c;
        char buffer[50]; // Buffer to store the token
        Vertex vertex;
        int vertex_count = 0;
        int triangle_count = 0;
        int position = -1;
        std::vector<Vertex> vertex_coords; // unordered_map each vertex to the index where it was seen
        std::vector<int> triangle_vertices; // Store the indices of each the vertex
        std::vector<Vertex> normal_coords; // For each triangle, store the normal vector
        int state = SOLID; // Start the state from solid
        if ((fp = fopen(filename, "r"))) { // Open the file
            while (!feof(fp) && state != 0 && state != ENDSOLID) { //Stop parsing if an unknown state is reached or 'endsolid' tag or end of file
                fscanf(fp, "%s", buffer); // Get each token from the file
                //printf("%s\n",buffer);
                //printf("State: %d\n",state);
                switch (state) { // Check the state to which the token belongs to
                case SOLID:
                    state = compareToken("solid", buffer, state); // Check whether the expected token and the input matches
                    while ((c = fgetc(fp)) != '\n') { // Go to the end of the line - In .stl files there are comments present in the first line
                    } // Skip all of that!
                    state = FACET; // If it has reached here, then it means that an Unknown state was not reached and can proceed to next expected state
                    break;
                case FACET:
                    state = compareToken("facet", buffer, state); // If 'facet' token is obtained then move on to the 'normal'
                    state = NORMAL;
                    break;
                case NORMAL:
                    state = compareToken("normal", buffer, state);
                    fscanf(fp, "%lf %lf %lf", &vertex.x, &vertex.y, &vertex.z); // Get the normal vector
                    //printf("%.1lf %.1lf %.1lf\n",vertex.x, vertex.y, vertex.z);
                    normal_coords.push_back(vertex); //Push this vector into the Vector
                    state = OUTER;
                    break;
                case OUTER:
                    state = compareToken("outer", buffer, state);
                    state = LOOP;
                    break;
                case LOOP:
                    state = compareToken("loop", buffer, state);
                    state = VERTEX;
                    break;
                case VERTEX:
                    state = compareToken("vertex", buffer, state);
                    fscanf(fp, "%lf %lf %lf", &vertex.x, &vertex.y, &vertex.z);
                    position = findVertex(vertex_coords, vertex);
                    if (position == -1) { // Check if vertex already exists in the unordered_map, if it does not
                        vertex_coords.push_back(vertex); // then insert this vertex and assign the position where it has been found as its value
                        position = vertex_coords.size();
                    }
                    //printf("%.1lf %.1lf %.1lf\n",vertex.x, vertex.y, vertex.z);
                    triangle_vertices.push_back(position); //Add current parsed vertex into the set of triangle vertices
                    vertex_count++; // Once three vertices has been parsed in this triangle, move to the next state
                    if (vertex_count == 3) {
                        vertex_count = 0;
                        state = ENDLOOP;
                    }
                    break;
                case ENDLOOP:
                    state = compareToken("endloop", buffer, state);
                    state = ENDFACET;
                    break;
                case ENDFACET:
                    // If 'endfacet' is reached then there are two options, either next can be 'facet' or 'endsolid'
                    state = compareToken("endfacet", buffer, state);
                    triangle_count++;
                    state = BRANCH;
                    break;
                case BRANCH:
                    if (compareToken("facet", buffer, state)) {
                        state = NORMAL; // if the token is 'facet' then start expecting 'normal'
                    }
                    else if (compareToken("endsolid", buffer, state)) {
                        state = ENDSOLID; // if the token is 'endsolid' we should stop parsing - while loop will take care of this
                    }
                    else {
                        state = UNKNOWN; // Any other token, it's definitely an Unknown state!
                    }
                    break;
                default:
                    state = UNKNOWN;
                    break;
                }
            }
            fclose(fp);
            if (state == UNKNOWN) {
                printf("%s", "Please check the syntax of your file\n");
                return IMPROPER_SYNTAX;
            }
            else {
                // If the state was not unknown, then start storing the STL file
                return store_STL(mesh, triangle_count, vertex_coords, triangle_vertices, normal_coords);
            }
        }
        else {
            printf("%s", "There was an error opening the file\n");
            return FILE_OPENING_ERROR;
        }
    }
    else {
        printf("%s", "Please check if file format is correct\n");
        return INCORRECT_FILE_FORMAT;
    }
}

int main(int argc, char* argv[])
{
    if (argc == 2) {
        TriangleMesh mesh;
        init_mesh(&mesh);
        read_STL(argv[1], &mesh);
        test_mesh(&mesh);
        free_mesh(&mesh);
    }
    else {
        printf("%s", "Please check the number of arguments\n");
        return FAILURE;
    }
}
