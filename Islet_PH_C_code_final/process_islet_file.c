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

int main(int argc, char* argv[]){

    //char str1[20] = "C programming";
    //char* str2 = (char*)malloc(40*sizeof(char));

    //// copying str1 to str2
    //strcpy(str2, str1);
    //strcat(str2, "hellp");

    //printf("\n%s", str2);

    //exit(1);

     int file_len = 10000;

     char* fprefix = (char*)malloc(file_len*sizeof(char));

     strcpy(fprefix, argv[1]);

     char ad_cell_suffix[] = "advertices_PH.csv";
     char b_cell_suffix[] = "bvertices_PH.csv";

     char ad_edge_suffix[] = "adedges_PH.csv";
     char b_edge_suffix[] = "bedges_PH.csv";

     char ad_triangles_suffix[] = "adtriangles_PH.csv";
     char b_triangles_suffix[] = "btriangles_PH.csv";


     int b_cell_fname_len = strlen(fprefix) + strlen(b_cell_suffix) + 100;
     int ad_cell_fname_len = strlen(fprefix) + strlen(ad_cell_suffix) + 100;

     int b_edge_fname_len = strlen(fprefix) + strlen(b_edge_suffix) + 100;
     int ad_edge_fname_len = strlen(fprefix) + strlen(ad_edge_suffix) + 100;

     int b_triangles_fname_len = strlen(fprefix) + strlen(b_triangles_suffix) + 100;
     int ad_triangles_fname_len = strlen(fprefix) + strlen(ad_triangles_suffix) + 100;


     char* b_cell_fname;
     char* ad_cell_fname;

     char* b_edge_fname;
     char* ad_edge_fname;

     char* b_triangles_fname;
     char* ad_triangles_fname;


     b_cell_fname = (char*)malloc(b_cell_fname_len*sizeof(char));
     strcpy(b_cell_fname, fprefix);
     strcat(b_cell_fname, b_cell_suffix);

     ad_cell_fname = (char*)malloc(ad_cell_fname_len*sizeof(char));
     strcpy(ad_cell_fname, fprefix);
     strcat(ad_cell_fname, ad_cell_suffix);

     b_edge_fname = (char*)malloc(b_edge_fname_len*sizeof(char));
     strcpy(b_edge_fname, fprefix);
     strcat(b_edge_fname, b_edge_suffix);

     ad_edge_fname = (char*)malloc(ad_edge_fname_len*sizeof(char));
     strcpy(ad_edge_fname, fprefix);
     strcat(ad_edge_fname, ad_edge_suffix);

     ad_triangles_fname = (char*)malloc(ad_triangles_fname_len*sizeof(char));
     strcpy(ad_triangles_fname, fprefix);
     strcat(ad_triangles_fname, ad_triangles_suffix);

     b_triangles_fname = (char*)malloc(b_triangles_fname_len*sizeof(char));
     strcpy(b_triangles_fname, fprefix);
     strcat(b_triangles_fname, b_triangles_suffix);




     //printf("\n%s", b_cell_fname);
     //printf("\n%s", ad_cell_fname);

     //printf("\n%s", b_edge_fname);
     //printf("\n%s", ad_edge_fname);

     //printf("\nPress key to continue\n\n");
     //getchar();


     PAR** b_cells;
     PAR** ad_cells;
     VERT_ID b_len, ad_len;
     VERT_ID b_max_len, ad_max_len;
     
     b_max_len = 10;
     ad_max_len = 10;

     b_cells = (PAR**)malloc(b_max_len*sizeof(PAR*));
     ad_cells = (PAR**)malloc(ad_max_len*sizeof(PAR*));

     b_len = 0;
     ad_len = 0;

     PAR loc_x, loc_y;

     // READ ad_cells
     char* line = NULL;
     size_t len = 0;
     char* end;
     char* dist;
     FILE *fp = fopen(ad_cell_fname, "r");  
     if (fp == NULL){
          printf("\n ERROR %s", ad_cell_fname);
          perror("Unable to open file!");
          getchar();
          exit(1);
     }

     int col = 0;


     while(getline(&line, &len, fp) != -1) {

          if (ad_len == ad_max_len){

              ad_max_len += 100;
              ad_cells = (PAR**)realloc(ad_cells, ad_max_len*sizeof(PAR*));

          }

          ad_cells[ad_len] = (PAR*)malloc(2*sizeof(PAR));

          // Read first field
          dist = strtok(line, ",");

          col = 0;

          while(dist != NULL){

            ad_cells[ad_len][col++] = atof(dist);

            // Read next field
            dist = strtok(NULL, ",");
              
          }


          ad_len++;

     }

     ad_cells = (PAR**)realloc(ad_cells, ad_len*sizeof(PAR*));

     fclose(fp);


     // Read b-cells
     fp = fopen(b_cell_fname, "r");  
     if (fp == NULL){
          printf("\n ERROR %s", b_cell_fname);
          perror("Unable to open file!");
          getchar();
          exit(1);
     }

     col = 0;


     while(getline(&line, &len, fp) != -1) {

          if (b_len == b_max_len){

              b_max_len += 100;
              b_cells = (PAR**)realloc(b_cells, b_max_len*sizeof(PAR*));

          }
          b_cells[b_len] = (PAR*)malloc(2*sizeof(PAR));

          // Read first field
          dist = strtok(line, ",");


          col = 0;
          while(dist != NULL){


            b_cells[b_len][col++] = atof(dist);

            // Read next field
            dist = strtok(NULL, ",");
              
          }

          b_len++;

     }

     fclose(fp);

     b_cells = (PAR**)realloc(b_cells, b_len*sizeof(PAR*));


     // Read ad-edges
     EDGE_ID ad_edge_len = 0;
     EDGE_ID ad_edge_max_len = 100;

     EDGE_ID* ad_edge_list = (EDGE_ID*)malloc(ad_edge_max_len*sizeof(EDGE_ID));
     

     fp = fopen(ad_edge_fname, "r");  
     if (fp == NULL){
          printf("\n ERROR %s", ad_edge_fname);
          perror("Unable to open file!");
          getchar();
          exit(1);
     }

     col = 0;

     while(getline(&line, &len, fp) != -1) {

          col = 0;
          if (ad_edge_len == ad_edge_max_len){

              ad_edge_max_len += 100;
              ad_edge_list = (EDGE_ID*)realloc(ad_edge_list, ad_edge_max_len*sizeof(EDGE_ID));

          }

          // Read first field
          dist = strtok(line, ",");


          while(dist != NULL){


            ad_edge_list[ad_edge_len++] = atoi(dist);

            // Read next field
            dist = strtok(NULL, ",");
              
          }

     }

     fclose(fp);

     ad_edge_list = (EDGE_ID*)realloc(ad_edge_list, ad_edge_len*sizeof(EDGE_ID));

     // Read b-edges
     EDGE_ID b_edge_len = 0;
     EDGE_ID b_edge_max_len = 100;

     EDGE_ID* b_edge_list = (EDGE_ID*)malloc(b_edge_max_len*sizeof(EDGE_ID));
     fp = fopen(b_edge_fname, "r");  
     if (fp == NULL){
          printf("\n ERROR %s", b_edge_fname);
          perror("Unable to open file!");
          getchar();
          exit(1);
     }

     col = 0;

     while(getline(&line, &len, fp) != -1) {

          col = 0;
          if (b_edge_len == b_edge_max_len){

              b_edge_max_len += 100;
              b_edge_list = (EDGE_ID*)realloc(b_edge_list, b_edge_max_len*sizeof(EDGE_ID));

          }

          // Read first field
          dist = strtok(line, ",");


          while(dist != NULL){


            b_edge_list[b_edge_len++] = atoi(dist);

            // Read next field
            dist = strtok(NULL, ",");
              
          }

     }

     b_edge_list = (EDGE_ID*)realloc(b_edge_list, b_edge_len*sizeof(EDGE_ID));

     fclose(fp);

     // Create ad-edge length list
     EDGE_ID ad_par_len = 0;
     PAR* ad_par_list = (PAR*)malloc(((EDGE_ID)(ad_edge_len/2))*sizeof(PAR));

     EDGE_ID ptr = 0;

     while (ptr < ad_edge_len){

        EDGE_ID c1 = ad_edge_list[ptr++];
        EDGE_ID c2 = ad_edge_list[ptr++];

        PAR x1 = ad_cells[c1][0];
        PAR y1 = ad_cells[c1][1];

        PAR x2 = ad_cells[c2][0];
        PAR y2 = ad_cells[c2][1];

        PAR diff1 = (x1 - x2);
        PAR diff2 = (y1 - y2);

        diff1 = diff1*diff1;
        diff2 = diff2*diff2;

        PAR sum = diff1 + diff2;

        ad_par_list[ad_par_len++] = sum;

     }
     
     // Create b-edge length list
     EDGE_ID b_par_len = 0;
     PAR* b_par_list = (PAR*)malloc(((EDGE_ID)(b_edge_len/2))*sizeof(PAR));

     ptr = 0;


     while (ptr < b_edge_len){

        EDGE_ID c1 = b_edge_list[ptr++];
        EDGE_ID c2 = b_edge_list[ptr++];

        PAR x1 = b_cells[c1][0];
        PAR y1 = b_cells[c1][1];

        PAR x2 = b_cells[c2][0];
        PAR y2 = b_cells[c2][1];

        PAR diff1 = (x1 - x2);
        PAR diff2 = (y1 - y2);

        diff1 = diff1*diff1;
        diff2 = diff2*diff2;

        PAR sum = diff1 + diff2;

        b_par_list[b_par_len++] = sum;

     }


     //// Check cells
     //for (EDGE_ID jj = 0 ; jj < ad_len; jj++){

     //     printf("\nad: (%lf, %lf)", ad_cells[jj][0], ad_cells[jj][1]);

     //}
     //for (EDGE_ID jj = 0 ; jj < b_len; jj++){

     //     printf("\nb: (%lf, %lf)", b_cells[jj][0], b_cells[jj][1]);

     //}


     ////////////////////////
     // DO NOT SORT!!!!!!!!!!
     //// Sort by edge lengths
     //printf("\nSorting...");
     //mergeSort(ad_par_list, ad_edge_list, 0, ad_par_len-1);
     //mergeSort(b_par_list, b_edge_list, 0, b_par_len-1);
     //printf("\nSorted");
     ////////////////////////

     // Make ad_vertex to edge association
     EDGE_ID** ad_vert_to_edge;
     ad_vert_to_edge = (EDGE_ID**)malloc(ad_len*sizeof(EDGE_ID*));
     for (EDGE_ID ii = 0 ; ii < ad_len; ii++){
        ad_vert_to_edge[ii] = (EDGE_ID*)malloc(ad_len*sizeof(EDGE_ID));
        for (EDGE_ID jj = 0 ; jj < ad_len; jj++){
            ad_vert_to_edge[ii][jj] = ad_edge_len;
        }

     }

     // Go through ad_edge_list
     EDGE_ID zz = 0;
     EDGE_ID order = 0;
     while (zz < ad_edge_len){

        EDGE_ID c1 = ad_edge_list[zz++];
        EDGE_ID c2 = ad_edge_list[zz++];

        ad_vert_to_edge[c1][c2] = order;
        ad_vert_to_edge[c2][c1] = order;
        order++;

     }
     
     // Make b_vertex to edge association
     EDGE_ID** b_vert_to_edge;
     b_vert_to_edge = (EDGE_ID**)malloc(b_len*sizeof(EDGE_ID*));
     for (EDGE_ID ii = 0 ; ii < b_len; ii++){
        b_vert_to_edge[ii] = (EDGE_ID*)malloc(b_len*sizeof(EDGE_ID));
        for (EDGE_ID jj = 0 ; jj < b_len; jj++){
            b_vert_to_edge[ii][jj] = b_edge_len;
        }

     }

     // Go through b_edge_list
     zz = 0;
     order = 0;
     while (zz < b_edge_len){

        EDGE_ID c1 = b_edge_list[zz++];
        EDGE_ID c2 = b_edge_list[zz++];

        b_vert_to_edge[c1][c2] = order;
        b_vert_to_edge[c2][c1] = order;
        order++;

     }


     // ad-triangles
     EDGE_ID n_tri = 0;

     //printf("\nad cells %d", ad_len);

     PAR cell1[2], cell2[2], cell3[2];


     PAR bcell[2];

     EDGE_ID count_has_beta=0, count_not_have_beta=0;

     EDGE_ID** ad_triangles;
     EDGE_ID ad_triangles_len = 0;
     EDGE_ID ad_triangles_max_len = 10;

     ad_triangles = (EDGE_ID**)malloc(ad_triangles_max_len*sizeof(EDGE_ID*));

     // Go through triangles that can be formed with edges
     EDGE_ID pp = 0;

     while (pp < ad_edge_len){

        VERT_ID c1 = ad_edge_list[pp++];
        VERT_ID c2 = ad_edge_list[pp++];

        cell1[0] = ad_cells[c1][0];
        cell1[1] = ad_cells[c1][1];

        cell2[0] = ad_cells[c2][0];
        cell2[1] = ad_cells[c2][1];

        for (VERT_ID c3 = 0; c3 < ad_len; c3++){

            if (c3 == c1) continue;
            if (c3 == c2) continue;

            if (ad_vert_to_edge[c1][c3] == ad_edge_len) continue;
            if (ad_vert_to_edge[c2][c3] == ad_edge_len) continue;

            cell3[0] = ad_cells[c3][0];
            cell3[1] = ad_cells[c3][1];

            // Now check for containing b-cell
            int contain = 0;

            for (VERT_ID b1 = 0; b1 < b_len; b1++){

                bcell[0] = b_cells[b1][0];
                bcell[1] = b_cells[b1][1];

                if (check_containment(cell1, cell2, cell3, bcell)){
                  contain = 1;
                  break;
                }

            }

            if (!contain){
                EDGE_ID e1 = ad_vert_to_edge[c1][c2];
                EDGE_ID e2 = ad_vert_to_edge[c1][c3];
                EDGE_ID e3 = ad_vert_to_edge[c2][c3];

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

                if (ad_triangles_len == ad_triangles_max_len){

                    ad_triangles_max_len += 50;
                    ad_triangles = (EDGE_ID**)realloc(ad_triangles, ad_triangles_max_len*sizeof(EDGE_ID*));

                }

                ad_triangles[ad_triangles_len] = (EDGE_ID*)malloc(3*sizeof(EDGE_ID));

                ad_triangles[ad_triangles_len][0] = boundary[0];
                ad_triangles[ad_triangles_len][1] = boundary[1];
                ad_triangles[ad_triangles_len][2] = boundary[2];
                ad_triangles_len++;



                //printf("\nTriangle contains bcell");
                //getchar();
            }

        }

     }

     //printf("\nNumber of triangles %d", ad_triangles_len);
     //getchar();
     ad_triangles = (EDGE_ID**)realloc(ad_triangles, ad_triangles_len*sizeof(EDGE_ID*));

     if (ad_triangles_len > 1){
        mergeSort_tri(ad_triangles, 0, ad_triangles_len-1);
     }


     fp = fopen(ad_triangles_fname, "w");  

     for (EDGE_ID ii = 0; ii < ad_triangles_len; ii++){

          fprintf(fp, "%d,%d,%d\n"\
                            , ad_triangles[ii][0]\
                            , ad_triangles[ii][1]\
                            , ad_triangles[ii][2]\
                            );

     }

     fclose(fp);


     // b-triangles
     PAR adcell[2];

     EDGE_ID** b_triangles;
     EDGE_ID b_triangles_len = 0;
     EDGE_ID b_triangles_max_len = 10;

     b_triangles = (EDGE_ID**)malloc(b_triangles_max_len*sizeof(EDGE_ID*));

     // Go through triangles that can be formed with edges
     pp = 0;

     while (pp < b_edge_len){

        VERT_ID c1 = b_edge_list[pp++];
        VERT_ID c2 = b_edge_list[pp++];

        cell1[0] = b_cells[c1][0];
        cell1[1] = b_cells[c1][1];

        cell2[0] = b_cells[c2][0];
        cell2[1] = b_cells[c2][1];

        for (VERT_ID c3 = 0; c3 < b_len; c3++){

            if (c3 == c1) continue;
            if (c3 == c2) continue;

            if (b_vert_to_edge[c1][c3] == b_edge_len) continue;
            if (b_vert_to_edge[c2][c3] == b_edge_len) continue;

            cell3[0] = b_cells[c3][0];
            cell3[1] = b_cells[c3][1];

            // Now check for containing ad-cell
            int contain = 0;

            for (VERT_ID b1 = 0; b1 < ad_len; b1++){

                adcell[0] = ad_cells[b1][0];
                adcell[1] = ad_cells[b1][1];

                if (check_containment(cell1, cell2, cell3, adcell)){
                  contain = 1;
                  break;
                }

            }

            if (!contain){
                EDGE_ID e1 = b_vert_to_edge[c1][c2];
                EDGE_ID e2 = b_vert_to_edge[c1][c3];
                EDGE_ID e3 = b_vert_to_edge[c2][c3];

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

                if (b_triangles_len == b_triangles_max_len){

                    b_triangles_max_len += 50;
                    b_triangles = (EDGE_ID**)realloc(b_triangles, b_triangles_max_len*sizeof(EDGE_ID*));

                }

                b_triangles[b_triangles_len] = (EDGE_ID*)malloc(3*sizeof(EDGE_ID));

                b_triangles[b_triangles_len][0] = boundary[0];
                b_triangles[b_triangles_len][1] = boundary[1];
                b_triangles[b_triangles_len][2] = boundary[2];
                b_triangles_len++;

                //printf("\nTriangle contains adcell");
                //getchar();
            }

        }

     }

     b_triangles = (EDGE_ID**)realloc(b_triangles, b_triangles_len*sizeof(EDGE_ID*));

     if (b_triangles_len > 1){
        mergeSort_tri(b_triangles, 0, b_triangles_len-1);
     }

     fp = fopen(b_triangles_fname, "w");  

     for (EDGE_ID ii = 0; ii < b_triangles_len; ii++){

          fprintf(fp, "%d,%d,%d\n"\
                            , b_triangles[ii][0]\
                            , b_triangles[ii][1]\
                            , b_triangles[ii][2]\
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

