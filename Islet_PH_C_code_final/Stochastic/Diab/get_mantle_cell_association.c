//#define PY_SSIZE_T_CLEAN
//#include <Python.h>

#include <stdio.h>
#include <stdlib.h>
#include <omp.h>
#include <math.h>
#include <unistd.h>
#include <sys/types.h>
#include <time.h>
#include <string.h>
#include <unistd.h>
#include <assert.h>
#include <time.h> 

#include <pthread.h>

typedef unsigned int EDGE_ID;
typedef unsigned int VERT_ID;

typedef double  PAR;

/*
 * file1: loop_pts (locations)
 * file2: undead loops
 * file3: cluster_pts
 */

int main(int argc, char* argv[]){

    

     int fname_len1 = strlen(argv[1]);
     int fname_len2 = strlen(argv[1]);
     int fname_len3 = strlen(argv[2]);
     int fname_len4 = strlen(argv[3]);

     char* loop_pts_file = malloc((fname_len1+100)*sizeof(char));
     char* loop_edges_file = malloc((fname_len2+100)*sizeof(char));
     char* undead_file = malloc((fname_len3+100)*sizeof(char));
     char* cluster_pts_file = malloc((fname_len4+100)*sizeof(char));

     strcpy(loop_pts_file, argv[1]);
     strcpy(loop_edges_file, argv[2]);
     strcpy(undead_file, argv[3]);
     strcpy(cluster_pts_file, argv[4]);


     //////////////////////////////////////////////////////////////
     // Read points
     EDGE_ID locs_max_len = 10;
     EDGE_ID locs_len = 0;
     PAR* locs = malloc(locs_max_len*sizeof(PAR));
     
     FILE *fp = fopen(loop_pts_file, "r");  

     if (fp == NULL){
          perror("Unable to open loop pts file!");
          exit(1);
     }

     while(getline(&line, &len, fp) != -1) {

            if (locs_len == locs_max_len){

                locs_max_len += 100;
                locs = realloc(locs, locs_max_len*sizeof(PAR));

            }

            dist = strtok(line, ",");
            while(dist != NULL){

              dist_d = atod(dist);
              locs[locs_len++] = dist_d;

              dist = strtok(NULL, ",");

            }

     }

     fclose(fp);

     locs = realloc(locs, locs_len*sizeof(PAR));
     locs_max_len = locs_len;


     //////////////////////////////////////////////////////////////
     // Read edges
     EDGE_ID loop_edges_max_len = 10;
     EDGE_ID loop_edges_len = 0;
     EDGE_ID** loop_edges = malloc(loop_max_len*sizeof(EDGE_ID**));
     
     FILE *fp = fopen(loop_edges_file, "r");  

     if (fp == NULL){
          perror("Unable to open loop edge file!");
          exit(1);
     }

     EDGE_ID pos;

     while(getline(&line, &len, fp) != -1) {

            if (locs_edges_len == locs_edges_max_len){

                locs_edges_max_len += 100;
                loop_edges = realloc(loop_edges, locs_edges_max_len*sizeof(PAR));

            }

            loop_edges[loop_edges_len] = malloc(2*sizeof(EDGE_ID));

            dist = strtok(line, ",");
            pos = 0;

            while(dist != NULL){

              dist_d = atoi(dist);
              loop_edges[loop_edges_len][pos++] = dist_d;

              dist = strtok(NULL, ",");

            }
            
            loop_edges_len++;

     }

     fclose(fp);

     loop_edges = realloc(loop_edges, loop_edges_len*sizeof(PAR));
     loop_edges_max_len = loop_edges_len;


     //////////////////////////////////////////////////////////////
     // Read undead loops











}
