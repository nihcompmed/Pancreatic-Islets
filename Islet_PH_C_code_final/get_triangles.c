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


PAR distxy(PAR, PAR, PAR, PAR);

PAR det2by2(PAR M[2][2]);

int check_containment(PAR p0[2], PAR p1[2], PAR p2[2], PAR p[2]);

// MERGE SORT ALGORITHM
void mergeSort(PAR* , EDGE_ID* , EDGE_ID , EDGE_ID ) ;
void merge(PAR* , EDGE_ID* , EDGE_ID , EDGE_ID , EDGE_ID ) ;

// MERGE SORT Triangles
void mergeSort_tri(EDGE_ID** , EDGE_ID , EDGE_ID ) ;
void merge_tri(EDGE_ID** , EDGE_ID , EDGE_ID , EDGE_ID ) ;

/*
 * Input: file1, file2, file3
 * file1: Points of loop graph
 * file2: Edges of loop graph (HAVE TO BE PRE_SORTED BY LENGTHS)
 * file3: points to be looped around/cluster points
 * file4: filename to save triangles
 */

int main(int argc, char* argv[]){

     int fname_len1 = strlen(argv[1]);
     int fname_len2 = strlen(argv[2]);
     int fname_len3 = strlen(argv[3]);
     int fname_len4 = strlen(argv[4]);

     char* loop_pts_file = malloc((fname_len1+100)*sizeof(char));
     char* loop_edges_file = malloc((fname_len2+100)*sizeof(char));
     char* cluster_pts_file = malloc((fname_len3+100)*sizeof(char));
     char* loop_triangles_file = malloc((fname_len3+100)*sizeof(char));

     strcpy(loop_pts_file, argv[1]);
     strcpy(loop_edges_file, argv[2]);
     strcpy(cluster_pts_file, argv[3]);
     strcpy(loop_triangles_file, argv[4]);


     PAR** loop_pts;
     PAR** cluster_pts;

     VERT_ID cluster_len, loop_len;
     VERT_ID cluster_max, loop_max;
     
     loop_max = 10;
     cluster_max = 10;

     loop_pts = (PAR**)malloc(loop_max*sizeof(PAR*));
     cluster_pts = (PAR**)malloc(cluster_max*sizeof(PAR*));

     cluster_len = 0;
     loop_len = 0;

     PAR loc_x, loc_y;

     // READ loop_pts
     char* line = NULL;
     size_t len = 0;
     char* end;
     char* dist;
     FILE *fp = fopen(loop_pts_file, "r");  
     if (fp == NULL){
          printf("\n ERROR %s", loop_pts_file);
          perror("Unable to open loop_pts_file!");
          getchar();
          exit(1);
     }

     int col = 0;


     while(getline(&line, &len, fp) != -1) {

          if (loop_len == loop_max){

              loop_max += 100;
              loop_pts = (PAR**)realloc(loop_pts, loop_max*sizeof(PAR*));

          }

          loop_pts[loop_len] = (PAR*)malloc(2*sizeof(PAR));

          // Read first field
          dist = strtok(line, ",");

          col = 0;

          while(dist != NULL){

            loop_pts[loop_len][col++] = atof(dist);

            // Read next field
            dist = strtok(NULL, ",");
              
          }


          loop_len++;

     }

     loop_pts = (PAR**)realloc(loop_pts, loop_len*sizeof(PAR*));

     fclose(fp);


     // Read cluster-pts
     fp = fopen(cluster_pts_file, "r");  
     if (fp == NULL){
          printf("\n ERROR %s", cluster_pts_file);
          perror("Unable to open cluster_pts_file!");
          getchar();
          exit(1);
     }

     col = 0;


     while(getline(&line, &len, fp) != -1) {

          if (cluster_len == cluster_max){

              cluster_max += 100;
              cluster_pts = (PAR**)realloc(cluster_pts, cluster_max*sizeof(PAR*));

          }
          cluster_pts[cluster_len] = (PAR*)malloc(2*sizeof(PAR));

          // Read first field
          dist = strtok(line, ",");


          col = 0;
          while(dist != NULL){


            cluster_pts[cluster_len][col++] = atof(dist);

            // Read next field
            dist = strtok(NULL, ",");
              
          }

          cluster_len++;

     }

     fclose(fp);

     cluster_pts = (PAR**)realloc(cluster_pts, cluster_len*sizeof(PAR*));


     // Read loop-edges
     EDGE_ID loop_edge_len = 0;
     EDGE_ID loop_edge_max_len = 100;

     EDGE_ID* loop_edge_list = (EDGE_ID*)malloc(loop_edge_max_len*sizeof(EDGE_ID));
     

     fp = fopen(loop_edges_file, "r");  
     if (fp == NULL){
          printf("\n ERROR %s", loop_edges_file);
          perror("Unable to open loop_edges_file!");
          getchar();
          exit(1);
     }

     col = 0;

     while(getline(&line, &len, fp) != -1) {

          col = 0;
          if (loop_edge_len == loop_edge_max_len){

              loop_edge_max_len += 100;
              loop_edge_list = (EDGE_ID*)realloc(loop_edge_list, loop_edge_max_len*sizeof(EDGE_ID));

          }

          // Read first field
          dist = strtok(line, ",");


          while(dist != NULL){


            loop_edge_list[loop_edge_len++] = atoi(dist);

            // Read next field
            dist = strtok(NULL, ",");
              
          }

     }

     fclose(fp);

     loop_edge_list = (EDGE_ID*)realloc(loop_edge_list, loop_edge_len*sizeof(EDGE_ID));


     // Create loop-edge length list
     EDGE_ID loop_par_len = 0;
     PAR* loop_par_list = (PAR*)malloc(((EDGE_ID)(loop_edge_len/2))*sizeof(PAR));

     EDGE_ID ptr = 0;

     while (ptr < loop_edge_len){

        EDGE_ID c1 = loop_edge_list[ptr++];
        EDGE_ID c2 = loop_edge_list[ptr++];

        PAR x1 = loop_pts[c1][0];
        PAR y1 = loop_pts[c1][1];

        PAR x2 = loop_pts[c2][0];
        PAR y2 = loop_pts[c2][1];

        PAR diff1 = (x1 - x2);
        PAR diff2 = (y1 - y2);

        diff1 = diff1*diff1;
        diff2 = diff2*diff2;

        PAR sum = diff1 + diff2;

        loop_par_list[loop_par_len++] = sum;

     }
     

     ////////////////////////
     // DO NOT SORT EDGES!!!!!!!!!!
     // THE EDGES ARE ALREADY SORTED BY LENGTHS
     //// Sort by edge lengths
     //printf("\nSorting...");
     //mergeSort(ad_par_list, ad_edge_list, 0, ad_par_len-1);
     //mergeSort(b_par_list, b_edge_list, 0, b_par_len-1);
     //printf("\nSorted");
     ////////////////////////

     // Initiate loop_vertex to edge association
     EDGE_ID** loop_vert_to_edge;
     loop_vert_to_edge = (EDGE_ID**)malloc(loop_len*sizeof(EDGE_ID*));

     for (EDGE_ID ii = 0 ; ii < loop_len; ii++){

        loop_vert_to_edge[ii] = (EDGE_ID*)malloc(loop_len*sizeof(EDGE_ID));

        for (EDGE_ID jj = 0 ; jj < loop_len; jj++){
            loop_vert_to_edge[ii][jj] = loop_edge_len;
        }

     }

     // Go through ad_edge_list
     EDGE_ID zz = 0;
     EDGE_ID order = 0;
     while (zz < loop_edge_len){

        EDGE_ID c1 = loop_edge_list[zz++];
        EDGE_ID c2 = loop_edge_list[zz++];

        loop_vert_to_edge[c1][c2] = order;
        loop_vert_to_edge[c2][c1] = order;
        order++;

     }
     
     // loop-triangles
     EDGE_ID n_tri = 0;

     //printf("\nad cells %d", loop_len);

     PAR cell1[2], cell2[2], cell3[2];


     PAR cluster_pt[2];

     EDGE_ID count_has_beta=0, count_not_have_beta=0;

     EDGE_ID** loop_triangles;
     EDGE_ID loop_triangles_len = 0;
     EDGE_ID loop_triangles_max_len = 10;

     loop_triangles = (EDGE_ID**)malloc(loop_triangles_max_len*sizeof(EDGE_ID*));

     // Go through triangles that can be formed with edges
     EDGE_ID pp = 0;
     
     for (VERT_ID c1 = 0; c1 < loop_len; c1++){


        cell1[0] = loop_pts[c1][0];
        cell1[1] = loop_pts[c1][1];

        for (VERT_ID c2 = c1+1; c2 < loop_len; c2++){

            if (loop_vert_to_edge[c1][c2] == loop_edge_len) continue;

            cell2[0] = loop_pts[c2][0];
            cell2[1] = loop_pts[c2][1];

            for (VERT_ID c3 = c2+1; c3 < loop_len; c3++){

                if (loop_vert_to_edge[c1][c3] == loop_edge_len) continue;
                if (loop_vert_to_edge[c2][c3] == loop_edge_len) continue;

                cell3[0] = loop_pts[c3][0];
                cell3[1] = loop_pts[c3][1];

                // Now check for containing b-cell
                int contain = 0;

                for (VERT_ID b1 = 0; b1 < cluster_len; b1++){

                    cluster_pt[0] = cluster_pts[b1][0];
                    cluster_pt[1] = cluster_pts[b1][1];

                    if (check_containment(cell1, cell2, cell3, cluster_pt)){
                      contain = 1;
                      break;
                    }

                }

                if (!contain){
                    EDGE_ID e1 = loop_vert_to_edge[c1][c2];
                    EDGE_ID e2 = loop_vert_to_edge[c1][c3];
                    EDGE_ID e3 = loop_vert_to_edge[c2][c3];

                    EDGE_ID boundary[3];

                    if (e1 < e2){
                        if (e1 < e3){
                            boundary[0] = e1;
                            if (e2 < e3){//e1 < e2 < e3
                                boundary[1] = e2;
                                boundary[2] = e3;
                            }
                            else{// e1 < e3 < e2
                                boundary[1] = e3;
                                boundary[2] = e2;
                            }
                        }
                        else{// e3 < e1 < e2

                            boundary[0] = e3;
                            boundary[1] = e1;
                            boundary[2] = e2;

                        }
                    }
                    else{//e2 < e1
                        if (e1 < e3){//e2 < e1 < e3
                            boundary[0] = e2;
                            boundary[1] = e1;
                            boundary[2] = e3;
                        }
                        else{//e3 < e1
                            if(e2 < e3){//e2 < e3 < e1
                                boundary[0] = e2;
                                boundary[1] = e3;
                                boundary[2] = e1;
                            }
                            else{// e3 < e2 < e1

                                boundary[0] = e3;
                                boundary[1] = e2;
                                boundary[2] = e1;
                            }

                        }

                    }

                    if (loop_triangles_len == loop_triangles_max_len){

                        loop_triangles_max_len += 50;
                        loop_triangles =\
                                  realloc(loop_triangles, loop_triangles_max_len*sizeof(EDGE_ID*));

                    }

                    loop_triangles[loop_triangles_len] = malloc(3*sizeof(EDGE_ID));

                    loop_triangles[loop_triangles_len][0] = boundary[0];
                    loop_triangles[loop_triangles_len][1] = boundary[1];
                    loop_triangles[loop_triangles_len][2] = boundary[2];
                    loop_triangles_len++;

                }

            }

        }

     }


     //printf("\nNumber of triangles %d", loop_triangles_len);
     //getchar();
     loop_triangles = (EDGE_ID**)realloc(loop_triangles, loop_triangles_len*sizeof(EDGE_ID*));

     if (loop_triangles_len > 1){
        mergeSort_tri(loop_triangles, 0, loop_triangles_len-1);
     }


     fp = fopen(loop_triangles_file, "w");  

     for (EDGE_ID ii = 0; ii < loop_triangles_len; ii++){

          fprintf(fp, "%d,%d,%d\n"\
                            , loop_triangles[ii][0]\
                            , loop_triangles[ii][1]\
                            , loop_triangles[ii][2]\
                            );

     }

     fclose(fp);


}


PAR det2by2(PAR M[2][2]){
    return M[0][0]*M[1][1] - M[0][1]*M[1][0];
}

int check_containment(PAR p0[2], PAR p1[2], PAR p2[2], PAR p[2]){

    //v1 = p1 - p0
    //v2 = p2 - p0
    
    PAR v1[2], v2[2];

    v1[0] =  p1[0] - p0[0];
    v1[1] =  p1[1] - p0[1];

    v2[0] =  p2[0] - p0[0];
    v2[1] =  p2[1] - p0[1];


    PAR M[2][2];

    M[0][0] = v1[0];
    M[1][0] = v1[1];

    M[0][1] = v2[0];
    M[1][1] = v2[1];
    
    PAR bottom = det2by2(M);

    if (bottom == 0){
        //#print('degenerate triangle')
        //#exit()
        return 0;
    }

    PAR v3[2];

    //c = p - p0
    v3[0] =  p[0] - p0[0];
    v3[1] =  p[1] - p0[1];



    //s_det = det2by2(np.vstack((c, v2)).T)

    M[0][0] = v3[0];
    M[1][0] = v3[1];

    M[0][1] = v2[0];
    M[1][1] = v2[1];
    PAR s_det = det2by2(M);

    //t_det = det2by2(np.vstack((v1, c)).T)

    M[0][0] = v1[0];
    M[1][0] = v1[1];

    M[0][1] = v3[0];
    M[1][1] = v3[1];
    PAR t_det = det2by2(M);

    PAR s = s_det/bottom;

    PAR t = t_det/bottom;

    if ((s >= 0) && (t >= 0)){
        if (s+t <= 1){
            return 1;
        }
    }

    return 0;
}

PAR distxy(PAR x1, PAR y1, PAR x2, PAR y2){
    PAR num1 = x1-x2;
    PAR num2 = y1-y2;

    return num1*num1 + num2*num2;
}

//////////////////////////////////////////////////////////
// MERGE SORT ALGORITHMS
//////////////////////////////////////////////////////////

// Merges two subarrays of arr[]. 
// First subarray is arr[l..m] 
// Second subarray is arr[m+1..r] 
void merge(PAR* arr, EDGE_ID* aux, EDGE_ID l, EDGE_ID m, EDGE_ID r) 
{ 
    EDGE_ID i, j, k; 
    EDGE_ID n1 = m - l + 1; 
    EDGE_ID n2 =  r - m; 
    //printf("\nn1, n2: %u, %u", n1, n2);
  
    /* create temp arrays */
    PAR *L, *R;
    L = (PAR*)malloc(n1*sizeof(PAR));
    R = (PAR*)malloc(n2*sizeof(PAR));

    /* create temp arrays */
    EDGE_ID* L_aux;
    EDGE_ID* R_aux;
    L_aux = (EDGE_ID*)malloc(2*n1*sizeof(EDGE_ID));
    R_aux = (EDGE_ID*)malloc(2*n2*sizeof(EDGE_ID));
    
    //int L[n1], R[n2]; 
  
    /* Copy data to temp arrays L[] and R[] */
    for (i = 0; i < n1; i++){ 
        L[i] = arr[l + i]; 
        L_aux[2*i] = aux[2*(l + i)]; 
        L_aux[2*i+1] = aux[2*(l + i) + 1]; 
    }

    for (j = 0; j < n2; j++) {
        R[j] = arr[m + 1+ j]; 
        R_aux[2*j] = aux[2*(m + 1+ j)]; 
        R_aux[2*j+1] = aux[2*(m + 1+ j)+1]; 
    }
  
    /* Merge the temp arrays back into arr[l..r]*/
    i = 0; // Initial index of first subarray 
    j = 0; // Initial index of second subarray 
    k = l; // Initial index of merged subarray 


    while (i < n1 && j < n2) 
    { 

          if (L[i] <= R[j]) 
          { 
              arr[k] = L[i]; 
	            aux[2*k] = L_aux[2*i];
	            aux[2*k+1] = L_aux[2*i+1];
              i++; 
          } 
          else
          { 
              arr[k] = R[j]; 
              aux[2*k] = R_aux[2*j]; 
              aux[2*k+1] = R_aux[2*j+1]; 
              j++; 
          } 

          k++;

    }

  
    /* Copy the remaining elements of L[], if there 
       are any */
    while (i < n1) 
    { 
        arr[k] = L[i]; 
        aux[2*k] = L_aux[2*i]; 
        aux[2*k+1] = L_aux[2*i+1]; 
        i++; 
        k++; 
    } 
  
    /* Copy the remaining elements of R[], if there 
       are any */
    while (j < n2) 
    { 
        arr[k] = R[j]; 
        aux[2*k] = R_aux[2*j]; 
        aux[2*k+1] = R_aux[2*j+1]; 
        j++; 
        k++; 
    } 

    free(L);
    free(R);
    free(L_aux);
    free(R_aux);
} 
  
/* l is for left index and r is right index of the 
   sub-array of arr to be sorted */
void mergeSort_tri(EDGE_ID** arr, EDGE_ID l, EDGE_ID r) 
{ 
    if (l < r) 
    { 
        // Same as (l+r)/2, but avoids overflow for 
        // large l and h 
        EDGE_ID m = l+(r-l)/2; 
  
        // Sort first and second halves 
        mergeSort_tri(arr, l, m); 
        mergeSort_tri(arr, m+1, r); 
  
        merge_tri(arr, l, m, r); 
    } 
} 

// Merges two subarrays of arr[]. 
// First subarray is arr[l..m] 
// Second subarray is arr[m+1..r] 
void merge_tri(EDGE_ID** arr, EDGE_ID l, EDGE_ID m, EDGE_ID r) 
{ 
    EDGE_ID i, j, k; 
    EDGE_ID n1 = m - l + 1; 
    EDGE_ID n2 =  r - m; 
    //printf("\nn1, n2: %u, %u", n1, n2);
  
    /* create temp arrays */
    EDGE_ID **L, **R;
    L = (EDGE_ID**)malloc(n1*sizeof(EDGE_ID*));
    R = (EDGE_ID**)malloc(n2*sizeof(EDGE_ID*));

    ///* create temp arrays */
    //EDGE_ID* L_aux;
    //EDGE_ID* R_aux;
    //L_aux = (EDGE_ID*)malloc(2*n1*sizeof(EDGE_ID));
    //R_aux = (EDGE_ID*)malloc(2*n2*sizeof(EDGE_ID));

    
    //int L[n1], R[n2]; 
  
    /* Copy data to temp arrays L[] and R[] */
    for (i = 0; i < n1; i++){ 
        L[i] = arr[l + i]; 
        //L_aux[2*i] = aux[2*(l + i)]; 
        //L_aux[2*i+1] = aux[2*(l + i) + 1]; 
    }

    for (j = 0; j < n2; j++) {
        R[j] = arr[m + 1+ j]; 
        //R_aux[2*j] = aux[2*(m + 1+ j)]; 
        //R_aux[2*j+1] = aux[2*(m + 1+ j)+1]; 
    }
  
    /* Merge the temp arrays back into arr[l..r]*/
    i = 0; // Initial index of first subarray 
    j = 0; // Initial index of second subarray 
    k = l; // Initial index of merged subarray 


    while (i < n1 && j < n2) 
    { 

          if (L[i][2] <= R[j][2]) 
          { 
              arr[k] = L[i]; 
	            //aux[2*k] = L_aux[2*i];
	            //aux[2*k+1] = L_aux[2*i+1];
              i++; 
          } 
          else
          { 
              arr[k] = R[j]; 
              //aux[2*k] = R_aux[2*j]; 
              //aux[2*k+1] = R_aux[2*j+1]; 
              j++; 
          } 

          k++;

    }

  
    /* Copy the remaining elements of L[], if there 
       are any */
    while (i < n1) 
    { 
        arr[k] = L[i]; 
        //aux[2*k] = L_aux[2*i]; 
        //aux[2*k+1] = L_aux[2*i+1]; 
        i++; 
        k++; 
    } 
  
    /* Copy the remaining elements of R[], if there 
       are any */
    while (j < n2) 
    { 
        arr[k] = R[j]; 
        //aux[2*k] = R_aux[2*j]; 
        //aux[2*k+1] = R_aux[2*j+1]; 
        j++; 
        k++; 
    } 

    free(L);
    free(R);
    //free(L_aux);
    //free(R_aux);
} 
  
/* l is for left index and r is right index of the 
   sub-array of arr to be sorted */
void mergeSort(PAR* arr, EDGE_ID* aux, EDGE_ID l, EDGE_ID r) 
{ 
    if (l < r) 
    { 
        // Same as (l+r)/2, but avoids overflow for 
        // large l and h 
        EDGE_ID m = l+(r-l)/2; 
  
        // Sort first and second halves 
        mergeSort(arr, aux, l, m); 
        mergeSort(arr, aux, m+1, r); 
  
        merge(arr, aux, l, m, r); 
    } 
} 


