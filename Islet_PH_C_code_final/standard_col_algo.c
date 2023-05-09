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

typedef struct{

    EDGE_ID g_cluster_len;
    EDGE_ID g_cluster_max_len;
    PAR** g_cluster_pts;

    EDGE_ID g_loop_len;
    EDGE_ID g_loop_max_len;
    PAR** g_loop_pts;

    EDGE_ID** g_loop_edges;
    EDGE_ID g_loop_edges_len;
    EDGE_ID g_loop_edges_max_len;

    EDGE_ID** g_RR;
    EDGE_ID* g_R_col_lens;
    EDGE_ID g_R_len;
    EDGE_ID g_R_max_len;

    EDGE_ID g_n_edges;

    EDGE_ID* g_pivots_H0;

    EDGE_ID* g_undead_edges;
    EDGE_ID g_undead_edges_len;
    EDGE_ID g_undead_edges_max_len;

    EDGE_ID* g_pivots;

    EDGE_ID** g_VV;

    EDGE_ID* g_V_col_lens;

    EDGE_ID* g_V_col_max_lens;

    EDGE_ID** g_undead_cycles;

    EDGE_ID g_undead_cycles_len;
    EDGE_ID g_undead_cycles_max_len;

    EDGE_ID* g_undead_cycle_len;
    EDGE_ID* g_undead_cycle_max_len;

    int g_suppress_output;

}filtration;


int cmpfunc (const void * a, const void * b) {
   return ( *(int*)a - *(int*)b );
}

int ray_algorithm(filtration* \
                  , EDGE_ID \
                  , EDGE_ID \
                  );

PAR min(PAR, PAR);
PAR max(PAR, PAR);

/*
 file1: loop-edges
 file2: triangles_file
 file3: undead file to save
 file4: cluster_pts files
 */

int main(int argc, char* argv[]){

      // Hard coding dimension for now
      int dim = 2;

     int fname_len1 = strlen(argv[1]);
     int fname_len2 = strlen(argv[2]);
     int fname_len3 = strlen(argv[3]);
     int fname_len4 = strlen(argv[4]);
     int fname_len5 = strlen(argv[5]);

     char* edge_file = malloc((fname_len1+100)*sizeof(char));
     char* triangles_file = malloc((fname_len2+100)*sizeof(char));
     char* undead_file = malloc((fname_len3+100)*sizeof(char));
     char* cluster_pts_file = malloc((fname_len4+100)*sizeof(char));
     char* loop_pts_file = malloc((fname_len5+100)*sizeof(char));

     strcpy(edge_file, argv[1]);
     strcpy(triangles_file, argv[2]);
     strcpy(undead_file, argv[3]);
     strcpy(cluster_pts_file, argv[4]);
     strcpy(loop_pts_file, argv[5]);

     filtration* self;

     self = malloc(sizeof(filtration));

     
     self->g_R_max_len = 10;
     self->g_RR = malloc(self->g_R_max_len*sizeof(EDGE_ID*));
     self->g_R_col_lens = malloc(self->g_R_max_len*sizeof(EDGE_ID));

     self->g_loop_edges_max_len = 10;
     self->g_loop_edges_len = 0;
     self->g_loop_edges = malloc(self->g_loop_edges_max_len*sizeof(EDGE_ID*));

     int pos;

     self->g_R_len = 0;

     EDGE_ID max_vert = 0;

     char* line = NULL;
     size_t len = 0;
     char* end;
     char* dist;
     EDGE_ID dist_d;

     FILE *fp = fopen(edge_file, "r");  

     if (fp == NULL){
          perror("Unable to open edge file!");
          exit(1);
     }

     while(getline(&line, &len, fp) != -1) {

            if (self->g_R_len == self->g_R_max_len){
                self->g_R_max_len += 100;
                self->g_RR = (EDGE_ID**)realloc(self->g_RR, self->g_R_max_len*sizeof(EDGE_ID*));
                self->g_R_col_lens = (EDGE_ID*)realloc(self->g_R_col_lens, self->g_R_max_len*sizeof(EDGE_ID));

            }

            if (self->g_loop_edges_len == self->g_loop_edges_max_len){

                self->g_loop_edges_max_len += 100;
                self->g_loop_edges = realloc(self->g_loop_edges\
                                            , self->g_loop_edges_max_len*sizeof(EDGE_ID*));

            }

            self->g_RR[self->g_R_len] = malloc(2*sizeof(EDGE_ID));
            self->g_R_col_lens[self->g_R_len] = 2;

            self->g_loop_edges[self->g_loop_edges_len] = malloc(2*sizeof(EDGE_ID));

            dist = strtok(line, ",");
            pos = 0;

            while(dist != NULL){

              dist_d = atoi(dist);

              self->g_RR[self->g_R_len][pos] = dist_d;
              self->g_loop_edges[self->g_loop_edges_len][pos] = dist_d;

              pos++;

              if (dist_d > max_vert){
                  max_vert = dist_d;
              }

              dist = strtok(NULL, ",");

            }

            self->g_R_len++;
            self->g_loop_edges_len++;

     }

     fclose(fp);


     self->g_RR = realloc(self->g_RR, self->g_R_len*sizeof(EDGE_ID*));
     self->g_R_col_lens = realloc(self->g_R_col_lens, self->g_R_len*sizeof(EDGE_ID));
     self->g_R_max_len = self->g_R_len;

     self->g_n_edges = self->g_R_len;

     self->g_loop_edges = realloc(self->g_loop_edges\
                                            , self->g_loop_edges_len*sizeof(EDGE_ID*));
     self->g_loop_edges_max_len = self->g_loop_edges_len;


     self->g_undead_edges_len = 0;
     self->g_undead_edges_max_len = 10;
     self->g_undead_edges = malloc(self->g_undead_edges_max_len*sizeof(EDGE_ID));


     // Because of 0-indexing
     max_vert = max_vert + 1;

     self->g_pivots_H0 = malloc(max_vert*sizeof(EDGE_ID));

     for (EDGE_ID ii = 0; ii < max_vert ; ii++){
        self->g_pivots_H0[ii] = self->g_n_edges;
     }

     // Define V
     self->g_VV = malloc(self->g_n_edges*sizeof(EDGE_ID*));
     self->g_V_col_lens = calloc(self->g_n_edges, sizeof(EDGE_ID));
     self->g_V_col_max_lens = calloc(self->g_n_edges, sizeof(EDGE_ID));

     // Standard column reduction
     EDGE_ID* temp;
     EDGE_ID* reduce_with;
     EDGE_ID reduce_with_len;
     EDGE_ID match_pivot;

     EDGE_ID this_pivot, this_len;

     for (EDGE_ID mm = 0; mm < self->g_n_edges; mm++){


        this_len = self->g_R_col_lens[mm];

        this_pivot = self->g_RR[mm][this_len-1];

        match_pivot = self->g_pivots_H0[this_pivot];

        while (match_pivot != self->g_n_edges){
              
          reduce_with = self->g_RR[match_pivot];
          reduce_with_len = self->g_R_col_lens[match_pivot];
            
          // Merge reduction operations

          // First check space
          if (self->g_V_col_lens[mm] + self->g_V_col_lens[match_pivot] + 1 > self->g_V_col_max_lens[mm]){
             
             self->g_V_col_max_lens[mm] += 5 + self->g_V_col_lens[match_pivot];
             if (!self->g_V_col_lens[mm]){
                 self->g_VV[mm] = malloc(self->g_V_col_max_lens[mm]*sizeof(EDGE_ID)); 
             }
             else{
                 self->g_VV[mm] = realloc(self->g_VV[mm], self->g_V_col_max_lens[mm]*sizeof(EDGE_ID)); 
             }

          }

          // Merge which column is added (because we do not store diagonal entries)
          self->g_VV[mm][self->g_V_col_lens[mm]++] = match_pivot;
          
          // Merge reduction operations of match_pivot
          for (int jj = 0; jj < self->g_V_col_lens[match_pivot] ; jj++){
             self->g_VV[mm][self->g_V_col_lens[mm]++] = self->g_VV[match_pivot][jj];
          }

          // Reduce
          temp = (EDGE_ID*)malloc((this_len+reduce_with_len)*sizeof(EDGE_ID));

          EDGE_ID ii = 0;
          EDGE_ID jj = 0;

          EDGE_ID ptr = 0;

          while((ii < this_len) && (jj < reduce_with_len)){

             
             if (self->g_RR[mm][ii] < reduce_with[jj]){
                 
                 temp[ptr++] = self->g_RR[mm][ii++];

             }
             else if (self->g_RR[mm][ii] > reduce_with[jj]){

                 temp[ptr++] = reduce_with[jj++];

             }
             else{
                 ii++;
                 jj++;
             }
             
             
          }

          while (ii < this_len){
                 temp[ptr++] = self->g_RR[mm][ii++];
          }

          while (jj < reduce_with_len){
                 temp[ptr++] = reduce_with[jj++];
          }

          this_len = ptr;

          if (this_len){
             
             this_pivot = temp[this_len-1];
             match_pivot = self->g_pivots_H0[this_pivot];

             self->g_RR[mm] = (EDGE_ID*)realloc(self->g_RR[mm], this_len*sizeof(EDGE_ID));
             self->g_R_col_lens[mm] = this_len;

             for (EDGE_ID kk = 0; kk < this_len; kk++){
                 
                 self->g_RR[mm][kk] = temp[kk];

             }

             free(temp);

          }
          else{
            
             // Reduced to 0
              
             if (self->g_undead_edges_len == self->g_undead_edges_max_len){
                
                self->g_undead_edges_max_len += 100;
                self->g_undead_edges = (EDGE_ID*)realloc(self->g_undead_edges\
                                                        , self->g_undead_edges_max_len*sizeof(EDGE_ID));

             }

             self->g_undead_edges[self->g_undead_edges_len++] = mm;
             
             free(temp);
             self->g_R_col_lens[mm] = 0;
             free(self->g_RR[mm]);
             break;

          }

        }

        if (self->g_R_col_lens[mm]){

           this_pivot = self->g_RR[mm][self->g_R_col_lens[mm]-1];
           
           self->g_pivots_H0[this_pivot] = mm;

        }

     }

     self->g_undead_edges = (EDGE_ID*)realloc(self->g_undead_edges\
                                              , self->g_undead_edges_len*sizeof(EDGE_ID));


     // Reduce each V column

     for (EDGE_ID mm = 0; mm < self->g_n_edges; mm++){
        
          if (self->g_V_col_lens[mm] < 2){
            continue;
          }

          // Sort
          qsort(self->g_VV[mm], self->g_V_col_lens[mm], sizeof(EDGE_ID), cmpfunc);

          // Reduce in place
          EDGE_ID cur_ptr = 0;
          int coeff = 1;
          for (EDGE_ID nn = 1; nn < self->g_V_col_lens[mm]; nn++){
 
            if (self->g_VV[mm][nn] == self->g_VV[mm][cur_ptr]){
                coeff = 1 - coeff;
            }
            else{
                if (coeff){
                    cur_ptr++;
                }
                self->g_VV[mm][cur_ptr] = self->g_VV[mm][nn];
                coeff = 1;

            }

          }

          self->g_V_col_lens[mm] = cur_ptr+1;

          self->g_VV[mm] = (EDGE_ID*)realloc(self->g_VV[mm], self->g_V_col_lens[mm]*sizeof(EDGE_ID));

     }


     //// Print V-cycles
     //for (EDGE_ID mm = 0; mm < self->g_R_len; mm++){

     //  if (self->g_V_col_lens[mm]>0){

     //     printf("\nReduction operations %d of length %d: ", mm, self->g_V_col_lens[mm]);

     //     for (EDGE_ID jj = 0; jj < self->g_V_col_lens[mm]; jj++){
     //        
     //        printf("%d,", self->g_VV[mm][jj]);

     //     }

     //  }
     //   
     //}


     // Free RR and R_col_lens
     for (EDGE_ID ii = 0; ii < self->g_n_edges; ii++){
          if (self->g_R_col_lens[ii]){
              free(self->g_RR[ii]);
          }
     }

     free(self->g_RR);
     free(self->g_R_col_lens);




 ////////////////////////////////////////////////////////
 ////////////////////////////////////////////////////////
 //
 //   COMPUTING H1
 //
 ////////////////////////////////////////////////////////
 ////////////////////////////////////////////////////////
 
     //printf("\nComputing H1...");
     //getchar();

     fp = fopen(triangles_file, "r");  

     if (fp == NULL){
          perror("Unable to open triangles file!");
          exit(1);
     }


     self->g_R_max_len = 10;
     self->g_RR = (EDGE_ID**)malloc(self->g_R_max_len*sizeof(EDGE_ID*));

     self->g_R_col_lens = (EDGE_ID*)malloc(self->g_R_max_len*sizeof(EDGE_ID));

     self->g_R_len = 0;

     pos = 0;

     while(getline(&line, &len, fp) != -1) {


            if (self->g_R_len == self->g_R_max_len){
                self->g_R_max_len += 100;
                self->g_RR = (EDGE_ID**)realloc(self->g_RR, self->g_R_max_len*sizeof(EDGE_ID*));
                self->g_R_col_lens = (EDGE_ID*)realloc(self->g_R_col_lens, self->g_R_max_len*sizeof(EDGE_ID));

            }

            self->g_RR[self->g_R_len] = (EDGE_ID*)malloc(3*sizeof(EDGE_ID));
            self->g_R_col_lens[self->g_R_len] = 3;

            dist = strtok(line, ",");

            pos = 0;

            while(dist != NULL){

              dist_d = atoi(dist);

              self->g_RR[self->g_R_len][pos++] = dist_d;

              dist = strtok(NULL, ",");
              
            }

            self->g_R_len++;

     }

     fclose(fp);


     self->g_RR = (EDGE_ID**)realloc(self->g_RR, self->g_R_len*sizeof(EDGE_ID*));
     self->g_R_col_lens = (EDGE_ID*)realloc(self->g_R_col_lens, self->g_R_len*sizeof(EDGE_ID));
     self->g_R_max_len = self->g_R_len;

     self->g_pivots = (EDGE_ID*)malloc((self->g_n_edges)*sizeof(EDGE_ID));

     for (EDGE_ID mm = 0; mm < self->g_n_edges; mm++){
        self->g_pivots[mm] = self->g_R_len;
     }



     // Standard column reduction
     for (EDGE_ID mm = 0; mm < self->g_R_len; mm++){

       this_len = self->g_R_col_lens[mm];

       this_pivot = self->g_RR[mm][this_len-1];

       match_pivot = self->g_pivots[this_pivot];

       while (match_pivot != self->g_R_len){

         reduce_with = self->g_RR[match_pivot];
         reduce_with_len = self->g_R_col_lens[match_pivot];

         // Reduce
         temp = (EDGE_ID*)malloc((this_len+reduce_with_len)*sizeof(EDGE_ID));

         EDGE_ID ii = 0;
         EDGE_ID jj = 0;

         EDGE_ID ptr = 0;

         while((ii < this_len) && (jj < reduce_with_len)){

            
            if (self->g_RR[mm][ii] < reduce_with[jj]){
                
                temp[ptr++] = self->g_RR[mm][ii++];

            }
            else if (self->g_RR[mm][ii] > reduce_with[jj]){

                temp[ptr++] = reduce_with[jj++];

            }
            else{
                ii++;
                jj++;
            }
            
            
         }

         while (ii < this_len){
                temp[ptr++] = self->g_RR[mm][ii++];
         }

         while (jj < reduce_with_len){
                temp[ptr++] = reduce_with[jj++];
         }

         this_len = ptr;

         if (this_len){
            
            this_pivot = temp[this_len-1];
            match_pivot = self->g_pivots[this_pivot];

            self->g_RR[mm] = (EDGE_ID*)realloc(self->g_RR[mm], this_len*sizeof(EDGE_ID));
            self->g_R_col_lens[mm] = this_len;

            for (EDGE_ID kk = 0; kk < this_len; kk++){
                
                self->g_RR[mm][kk] = temp[kk];

            }

            free(temp);

         }
         else{
            
            free(temp);
            self->g_R_col_lens[mm] = 0;
            free(self->g_RR[mm]);
            break;

         }
          
       }

       if (self->g_R_col_lens[mm]){

          this_pivot = self->g_RR[mm][self->g_R_col_lens[mm]-1];
          
          self->g_pivots[this_pivot] = mm;

       }

     }


     // Undead cycles -- if R(i) is zero and edge (i, j) is not pivot

     self->g_undead_cycles_max_len = 10;
     self->g_undead_cycles = (EDGE_ID**)malloc(self->g_undead_edges_max_len*sizeof(EDGE_ID*));
     self->g_undead_cycles_len = 0;
     self->g_undead_cycle_len = (EDGE_ID*)malloc(self->g_undead_edges_max_len*sizeof(EDGE_ID));
     
     // Get V-cycles for cycles that do not die
     
     for (EDGE_ID ii = 0; ii < self->g_undead_edges_len; ii++){

        //printf("\n%d", ii);

        EDGE_ID this_edge = self->g_undead_edges[ii];

        if (self->g_pivots[this_edge]==self->g_R_len){


            if (self->g_undead_cycles_len == self->g_undead_cycles_max_len){

                self->g_undead_cycles_max_len += 100;
                self->g_undead_cycles = (EDGE_ID**)realloc(self->g_undead_cycles\
                                                          , self->g_undead_cycles_max_len*sizeof(EDGE_ID*));
                self->g_undead_cycle_len = (EDGE_ID*)realloc(self->g_undead_cycle_len\
                                                          , self->g_undead_cycles_max_len*sizeof(EDGE_ID));

            }
            
            self->g_undead_cycle_len[self->g_undead_cycles_len] = self->g_V_col_lens[this_edge]+1;

            self->g_undead_cycles[self->g_undead_cycles_len] =\
                                            (EDGE_ID*)malloc(self->g_undead_cycle_len[self->g_undead_cycles_len]*sizeof(EDGE_ID));

            //printf("\nPrinting V col of:%d,", this_edge);

            for (EDGE_ID jj = 0; jj < self->g_V_col_lens[this_edge]; jj++){
               
               //printf("%d,", self->g_VV[this_edge][jj]);
               self->g_undead_cycles[self->g_undead_cycles_len][jj] =\
                                                        self->g_VV[this_edge][jj];

            }

            self->g_undead_cycles[self->g_undead_cycles_len][self->g_undead_cycle_len[self->g_undead_cycles_len]-1] =\
                                                                                                          this_edge;

            self->g_undead_cycles_len++;

        }

     }

     

     EDGE_ID max_c1, max_c2, max_diff;

     //printf("\ngot undead %d..Press key to continue...", self->g_undead_cycles_len);
     //getchar();

     if (self->g_undead_cycles_len > 1){

        while(1){

            max_diff = 0;

            for (EDGE_ID ii = 0; ii < self->g_undead_cycles_len; ii++){

                for (EDGE_ID jj = ii+1; jj < self->g_undead_cycles_len; jj++){

                     EDGE_ID ptr = 0;
                     EDGE_ID iii = 0;
                     EDGE_ID jjj = 0;

                     EDGE_ID diff1, diff2, diff;

                     diff1 = 0;
                     diff2 = 0;

                     while((iii < self->g_undead_cycle_len[ii]) && (jjj < self->g_undead_cycle_len[jj])){

                         if (self->g_undead_cycles[ii][iii] < self->g_undead_cycles[jj][jjj]){
                             
                             ptr++;
                             iii++;

                         }
                         else if (self->g_undead_cycles[ii][iii] > self->g_undead_cycles[jj][jjj]){

                             ptr++;
                             jjj++;

                         }
                         else{
                             iii++;
                             jjj++;
                         }

                     }

                     while (iii < self->g_undead_cycle_len[ii]){
                            ptr++;
                            iii++;
                     }

                     while (jjj < self->g_undead_cycle_len[jj]){
                            ptr++;
                            jjj++;
                     }

                     if (ptr < self->g_undead_cycle_len[ii]){
                           diff1 = self->g_undead_cycle_len[ii] - ptr;
                     }

                     if (ptr < self->g_undead_cycle_len[jj]){
                           diff2 = self->g_undead_cycle_len[jj] - ptr;
                     }

                     if ((diff1 > diff2) || (diff1 == diff2)){
                         if (diff1 > max_diff){
                             max_diff = diff1;
                             max_c1 = ii;
                             max_c2 = jj;
                         }
                     }
                     else if (diff2 > diff1){
                         if (diff2 > max_diff){
                             max_diff = diff2;
                             max_c1 = jj;
                             max_c2 = ii;
                         }
                     }

                }
                 
            }

            if (!max_diff) break;


            //printf("\nmax diff %d if %d - %d", max_diff, max_c1, max_c2);

            EDGE_ID len1 = self->g_undead_cycle_len[max_c1];
            EDGE_ID len2 = self->g_undead_cycle_len[max_c2];

            EDGE_ID* cyc1 = self->g_undead_cycles[max_c1];
            EDGE_ID* cyc2 = self->g_undead_cycles[max_c2];

            //printf("\nReducing:%d:%d", max_c1, max_c2);

            //printf(":");
            //for (EDGE_ID ii = 0; ii < len1; ii++){
            //   printf("%d,", cyc1[ii]);
            //}
            //printf(":");
            //for (EDGE_ID ii = 0; ii < len2; ii++){
            //   printf("%d,", cyc2[ii]);
            //}

            temp = (EDGE_ID*)malloc((len1+len2)*sizeof(EDGE_ID));

            EDGE_ID ii = 0;
            EDGE_ID jj = 0;
            EDGE_ID ptr = 0;

            while((ii < len1) && (jj < len2)){

                if (cyc1[ii] < cyc2[jj]){
                    temp[ptr++] = cyc1[ii++];
                }
                else if (cyc1[ii] > cyc2[jj]){
                    temp[ptr++] = cyc2[jj++];
                }
                else{
                    ii++;
                    jj++;
                }

            }

            while (ii < len1){
                   temp[ptr++] = cyc1[ii++];
            }

            while (jj < len2){
                   temp[ptr++] = cyc2[jj++];
            }

            self->g_undead_cycle_len[max_c1] = ptr;
            self->g_undead_cycles[max_c1] = (EDGE_ID*)realloc(self->g_undead_cycles[max_c1]\
                                                              , ptr*sizeof(EDGE_ID));

            for (EDGE_ID kk = 0; kk < ptr; kk++){
                   self->g_undead_cycles[max_c1][kk] = temp[kk];
                   //printf("%d,", self->g_undead_cycles[max_c1][kk]);
            }

            //printf(":");
            //for (EDGE_ID ii = 0; ii < ptr; ii++){
            //   printf("%d,", temp[ii]);
            //}

            free(temp);


        }

     }

     //printf("\nshortened undead ..Press key to continue...");
     //getchar();
     //


     // Now to get loop:cluster association and save it to file

     // First read the cluster pts
     self->g_cluster_max_len = 100;
     self->g_cluster_len = 0;
     self->g_cluster_pts = malloc(self->g_cluster_max_len*sizeof(PAR*));

     fp = fopen(cluster_pts_file, "r");
     if (fp == NULL){
          perror("Unable to open cluster locs file!");
          exit(1);
     }

     char* dummy;

     PAR val;

     while(getline(&line, &len, fp) != -1) {

          if (self->g_cluster_len == self->g_cluster_max_len){

                self->g_cluster_max_len += 100;
                self->g_cluster_pts = realloc(self->g_cluster_pts\
                                          , self->g_cluster_max_len*sizeof(PAR*));

          }

          self->g_cluster_pts[self->g_cluster_len] = malloc(dim*sizeof(PAR));

          dist = strtok(line, ",");
          pos = 0;

          while(dist != NULL){

              val = strtod(dist, &dummy);

              self->g_cluster_pts[self->g_cluster_len][pos++] = val;
                
              dist = strtok(NULL, ",");

          }

          self->g_cluster_len++;


     }

     self->g_cluster_pts = realloc(self->g_cluster_pts\
                                 , self->g_cluster_len*sizeof(PAR*));

     fclose(fp);


     // Read loop points
     
     //printf("\nREADING LOOP POINTS");

     self->g_loop_max_len = 100;
     self->g_loop_len = 0;
     self->g_loop_pts = malloc(self->g_loop_max_len*sizeof(PAR*));

     fp = fopen(loop_pts_file, "r");
     if (fp == NULL){
          perror("Unable to open loop locs file!");
          exit(1);
     }


     while(getline(&line, &len, fp) != -1) {

          if (self->g_loop_len == self->g_loop_max_len){

                self->g_loop_max_len += 100;
                self->g_loop_pts = realloc(self->g_loop_pts\
                                          , self->g_loop_max_len*sizeof(PAR*));

          }

          self->g_loop_pts[self->g_loop_len] = malloc(dim*sizeof(PAR));


          dist = strtok(line, ",");
          pos = 0;

          while(dist != NULL){

              val = strtod(dist, &dummy);
              //printf("\nadding %lf", val);

              self->g_loop_pts[self->g_loop_len][pos++] = val;
                
              dist = strtok(NULL, ",");

          }

          self->g_loop_len++;


     }


     self->g_loop_pts = realloc(self->g_loop_pts\
                                 , self->g_loop_len*sizeof(PAR*));

     fclose(fp);

     //printf("\ncluster pts %d, loop pts %d", self->g_cluster_len, self->g_loop_len);

     //// Check points loaded
     //for (EDGE_ID kk = 0; kk < self->g_cluster_len; kk++){

     //     printf("\n %lf, %lf", self->g_cluster_pts[kk][0]\
     //                       , self->g_cluster_pts[kk][1]);

     //}

     //for (EDGE_ID kk = 0; kk < self->g_loop_len; kk++){

     //     printf("\n %lf, %lf", self->g_loop_pts[kk][0]\
     //                       , self->g_loop_pts[kk][1]);

     //}
     //getchar();


     //printf("\nCluster len %d", self->g_cluster_len);
     //printf("\nloop len %d", self->g_loop_len);

     fp = fopen(undead_file, "w");

     //printf("\nWriting to %s", undead_file);

     for (EDGE_ID ii = 0; ii < self->g_undead_cycles_len; ii++){

          if (!self->g_undead_cycle_len[ii]) continue;

          //printf("\nLoop ");

          for (EDGE_ID jj = 0; jj < self->g_undead_cycle_len[ii]; jj++){
                
                fprintf(fp, "%d,", self->g_undead_cycles[ii][jj]);

                //printf("%d,", self->g_undead_cycles[ii][jj]);
          }

          fprintf(fp, ":");
          // CHECK CONTAINMENT
          for (EDGE_ID kk = 0; kk < self->g_cluster_len; kk++){

                //printf("\nChecking with cell %d", kk);

                int contains = ray_algorithm(self\
                                            , ii\
                                            , kk\
                                            );
                if (contains){
                    fprintf(fp, "%d,", kk);
                    //printf("\nContains cell %d,", kk);
                }

          }


          fprintf(fp, "\n");

          //getchar();
        
     }

     fclose(fp);


     // Print undead that wrap around cluster cells



}

int ray_algorithm(filtration* self\
                  , EDGE_ID loop_id\
                  , EDGE_ID cluster_cell_id\
                  ){


    PAR x = self->g_cluster_pts[cluster_cell_id][0];
    PAR y = self->g_cluster_pts[cluster_cell_id][1];

    //printf("\nCell is %lf, %lf", x, y);

    int intersect = 0;

    for (EDGE_ID ptr = 0; ptr < self->g_undead_cycle_len[loop_id]; ptr++){

        EDGE_ID edge = self->g_undead_cycles[loop_id][ptr];

        EDGE_ID p1 = self->g_loop_edges[edge][0];
        EDGE_ID p2 = self->g_loop_edges[edge][1];

        PAR x1 = self->g_loop_pts[p1][0];
        PAR y1 = self->g_loop_pts[p1][1];
          
        PAR x2 = self->g_loop_pts[p2][0];
        PAR y2 = self->g_loop_pts[p2][1];

        //printf("\nEdge (%lf, %lf) - (%lf, %lf)"\
        //                                    , x1\
        //                                    , y1\
        //                                    , x2\
        //                                    , y2\
        //                                    );

        // If horizontal line, check if point is on the edge
        if (y1 == y2){
            PAR xmin = min(x1, x2);
            PAR xmax = max(x1, x2);
            if ((x > xmin && x < xmax) || (x == xmin) || (x == xmax)){
                //printf("\nReturning 1 because on horizontal edge");
                return 1;
            }

        }
        else{

            // Otherwise check intersection (horizontal ray algorithm)
            PAR t = (y-y1)/(y2-y1);
            if ((t > 0) && (t < 1)){
                PAR x_intersect = (1-t)*x1 + t*x2;
                if (x_intersect >= x){
                    //printf("\nAdding 1 because intersects inside");
                    intersect += 1;
                }
            }
            else if (t == 0){
                //printf("\nAdding 0.5 because intersects at end-point");
                intersect += 0.5;
            }
            else if (t == 1){
                //printf("\nAdding 0.5 because intersects at end-point");
                intersect += 0.5;
            }

        }


    }

    //printf("\nIntersection is %d", intersect);

    if (intersect%2==0)
        return 0;
    else
        return 1;


}

PAR min(PAR a, PAR b){

    if (a > b){
      return b;
    }
    else if (a < b){
      return a;
    }
    else{
      return a;
    }

}

PAR max(PAR a, PAR b){

    if (a > b){
      return a;
    }
    else if (a < b){
      return b;
    }
    else{
      return a;
    }

}

//double isLeft( double P0x,double P0y, double P1x,double P1y, double P2x, double P2y )
// {
//   //isLeft(): tests if a point is Left|On|Right of an infinite line.
//   //    Input:  three points P0, P1, and P2
//   //    Return: >0 for P2 left of the line through P0 and P1
//   //            =0 for P2  on the line
//   //            <0 for P2  right of the line
//   //    See: Algorithm 1 "Area of Triangles and Polygons"
//   //algorithm from http://geomalgorithms.com/a03-_inclusion.html
//   return ( (P1x - P0x) * (P2y - P0y)
//            - (P2x -  P0x) * (P1y - P0y) );
// }
//
//
////int WindingNumber(Graph* Tree, vector<vertex_t>* Loop,double Px,double Py ){
//int WindingNumber(filtration* self\
//                  , EDGE_ID loop_id\
//                  , EDGE_ID cluster_cell_id\
//                  ){
//
//  PAR x = self->g_cluster_pts[cluster_cell_id][0];
//  PAR y = self->g_cluster_pts[cluster_cell_id][1];
//
//  //Gp = Tree;
//  // wn_PnPoly(): winding number test for a point in a polygon
//  //      Input:   P = a point,
//  //               V[] = vertex points of a polygon V[n+1] with V[n]=V[0]
//  //      Return:  wn = the winding number (=0 only when P is outside)
//  //algorithm from http://geomalgorithms.com/a03-_inclusion.html
//  int wn = 0;    // the  winding number counter
//  
//  // loop through all edges of the polygon
//  //cout << "Inside winding number " << endl;
//  for (int i=0; i<(*Loop).size()-1; i++) {   // edge from V[i] to  V[i+1]
//    if ((*Tree)[(*Loop)[i]].y <= Py) {          // start y <= P.y
//	    if ((*Tree)[(*Loop)[i+1]].y  > Py)      // an upward crossing
//	      if (isLeft((*Tree)[(*Loop)[i]].x\
//                          ,(*Tree)[(*Loop)[i]].y\
//                          ,(*Tree)[(*Loop)[i+1]].x\
//                          ,(*Tree)[(*Loop)[i+1]].y\
//                          ,Px\
//                          ,Py) > 0)  // P left of  edge
//	          ++wn;            // have  a valid up intersect
//    }
//    else {                        // start y > P.y (no test needed)
//	    if ((*Tree)[(*Loop)[i+1]].y  <= Py)     // a downward crossing
//	      if (isLeft((*Tree)[(*Loop)[i]].x\
//              ,(*Tree)[(*Loop)[i]].y\
//              ,(*Tree)[(*Loop)[i+1]].x\
//              ,(*Tree)[(*Loop)[i+1]].y\
//              ,Px\
//              ,Py ) < 0)  // P right of  edge
//	        --wn;            // have  a valid down intersect
//    }
//  }
//  //cout << "Finished winding number " << endl;
//  return wn;
//}
//
//
