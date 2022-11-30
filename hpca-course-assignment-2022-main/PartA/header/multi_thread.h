#include <pthread.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <pthread.h>
#include <smmintrin.h>
#include <immintrin.h>
#include <stdint.h>
#include <string.h>


struct arg_struct {
    int N;
    int *matA;
    int *matB;
    int *output;
    int i;
    int *colsadded;
    int j;
};

void* rowaddition(void* arguments)
{
    //cout<<"entered row addition function"<<endl;
    struct arg_struct *args = (struct arg_struct *)arguments;
    //cout<<"entered row addition function"<<args->i<<endl;
    //int i = step_i++; //i denotes row number of resultant matC
    
    //int rowcounter = args->i;
   for (int ii = 0; ii < args->N; ii+=8){
      //  cout<<"in the for loop"<<args->i<<endl;
    
        __m256i firstrow = _mm256_loadu_si256((__m256i*)&args->matA[(args->i) * (args->N) + ii]);
        __m256i secondrow = _mm256_loadu_si256((__m256i*)&args->matA[((args->i) + 1) * (args->N) + ii]);
       // union{
        //	__m256i addedrows;
        //	int s[8];
        //};

       __m256i addedrows = _mm256_add_epi32(firstrow,secondrow);
        //cout<<0;
        //uint16_t val[8];
    	//memcpy(val, &addedrows, sizeof(val));
    	//printf("Numerical: %i %i %i %i %i %i %i %i \n", 
          // val[0], val[1], val[2], val[3], val[4], val[5], 
           //val[6], val[7]);
        
        //for(int i=0; i<8; i++)
        //	cout << s[i] << " ";
	//cout << "\n";
        _mm256_storeu_si256((__m256i*) &args->output[(args->i)/2 * (args->N) + ii], addedrows);
        //cout << "Thread " << args->i << "\n";
        //if (args->i == 1) {
        //cout << "Thread 0:\n";
        //for(int i=0; i<8; i++){
    	//	cout<<args->output[args->i*args->N + ii+i]<<" ";
    	//}
    	}
    //cout<<endl;
        //cout<<1;
   }
    //rowcounter++;
   
    
      



void* coladdition(void* arguments)
{
    struct arg_struct *args = (struct arg_struct *)arguments;
    
    //cout<<endl;
    for (int j = 0; j < args->N; j+=8)
    {
        __m256i firstcol = _mm256_setr_epi32(args->matB[j*args->N + args->i], args->matB[(j+1)*args->N + args->i], args->matB[(j+2)*args->N + args->i], args->matB[(j+3)*args->N + args->i], args->matB[(j+4)*args->N + args->i], args->matB[(j+5)*args->N + args->i], args->matB[(j+6)*args->N + args->i], args->matB[(j+7)*args->N + args->i]);
          __m256i secondcol = _mm256_setr_epi32(args->matB[j*args->N + args->i+1], args->matB[(j+1)*args->N + args->i+1], args->matB[(j+2)*args->N + args->i+1], args->matB[(j+3)*args->N + args->i+1], args->matB[(j+4)*args->N + args->i+1], args->matB[(j+5)*args->N + args->i+1], args->matB[(j+6)*args->N + args->i+1], args->matB[(j+7)*args->N + args->i+1]);
          
          __m256i addedcols = _mm256_add_epi32(firstcol,secondcol);
          _mm256_storeu_si256((__m256i*) &args->output[(args->i/2)*args->N + j], addedcols);

          //for (int i = 0; i < 8; i++)
          //{
            //cout<<args->output[(args->i/2)*args->N + i]<<endl;
            
            
          //}
          
    }
    
    
}


void* multiplication(void* arguments)
{
    struct arg_struct *args = (struct arg_struct *)arguments;
    for (int k = 0; k < args->N/2; k=k+1)
    {
        for(int jj=0; jj<args->N; jj+=8){
      __m256i firstmultiplier = _mm256_loadu_si256((__m256i*)&args->matA[args->i*args->N+jj]);
      __m256i secondmultiplier = _mm256_loadu_si256((__m256i*)&args->colsadded[k*args->N+jj]);
      
      __m256i results = _mm256_mullo_epi32(firstmultiplier, secondmultiplier);
      
        int finalres = _mm256_extract_epi32(results, 0) + _mm256_extract_epi32(results, 1) + _mm256_extract_epi32(results, 2) + _mm256_extract_epi32(results, 3) + _mm256_extract_epi32(results, 4) + _mm256_extract_epi32(results, 5) + _mm256_extract_epi32(results, 6) + _mm256_extract_epi32(results, 7);
     // args->output[(args->i)>>2*(k>>2)] += finalres;
     args->output[args->i*args->N/2 + k]  += finalres;  
      }
    }
    
}

// Fill in this function
void multiThread(int N, int *matA, int *matB, int *output)
{
	int *coladdedarr = new int[N*(N>>1)];
	int *rowaddedarr = new int[N*(N>>1)];
    pthread_t rowadderthreads[N/2];
    struct arg_struct args_ptrs[N/2];
    for (int iii = 0; iii < N; iii+=2) {
        struct arg_struct args;
        args.matA = matA;
        args.matB = matB;
        args.N = N;
        args.output = rowaddedarr;
        args.i = iii;
        args_ptrs[iii/2] = args;
    }
     for (int iii = 0; iii < N; iii+=2) {
    
   	 pthread_create(&rowadderthreads[iii/2], NULL, rowaddition, (void*) &args_ptrs[iii/2]);
    
    }
    
    for (int i = 0; i < N/2; i++)
        pthread_join(rowadderthreads[i], NULL);
    
  //  for(int i=0; i<N*N/2; i++){
   // 	cout<<rowaddedarr[i]<<" ";
   // }
    //cout<<"join done!\n";    
        
    //for(int i=0; i<N; i++){
    	//cout<<matA[i]<<" ";
    //}
    //cout<<endl;
    //cout<<"rows added"<<endl;
    
    pthread_t coladderthreads[N>>1];
    //int count=0;
    for (int ii = 0; ii < N; ii+=2)
    {
        struct arg_struct args;
        args.matA = matA;
        args.matB = matB;
        args.N = N;
        args.output = coladdedarr;
        args.i = ii;
        args_ptrs[ii/2] = args;
        //count++;
    }
    //cout<<"colls added"<<endl;
    
    for (int ii = 0; ii < N; ii+=2)
    {
    	pthread_create(&coladderthreads[ii/2], NULL, coladdition, (void*) &args_ptrs[ii/2]);
	}
    for (int i = 0; i < N/2; i++)
        pthread_join(coladderthreads[i], NULL);      
    // for(int i=0; i<N*N/2; i++){
    //	cout<<coladdedarr[i]<<" ";
    //}
    
    //for(int i=0; i<N; i++){
    	//cout<<coladdedarr[i]<<" ";
    //}
    //cout<<endl;
    //int count = 0;
    pthread_t multiplier[(N>>1)];
    for (int ii = 0; ii < (N>>1); ii++)
    {
 	
        struct arg_struct args;
        args.matA = rowaddedarr;
        args.matB = matB;
        args.N = N;
        args.output = output;
        args.colsadded = coladdedarr;
        args.i = ii;
       // args.j = jj;
       	args_ptrs[ii] = args;
        //count++;
        
    }
    for (int ii = 0; ii < N>>1; ii++)
    {
    	
     pthread_create(&multiplier[ii], NULL, multiplication, (void*) &args_ptrs[ii]);
     
     }
    for (int i = 0; i < (N>>1); i++)
        pthread_join(multiplier[i], NULL);   
        
    //for(int i=0; i<N/2*N/2; i++){
    	//cout<<output[i]<<" ";
    //}
}



/*
int BARRIER = 0;



void *workerThreadFunc(void *arg)
{
    while (BARRIER == 0){

    }
    ThreadArgs * arg_p = (ThreadArgs *) arg;
    long *tid = (long *) arg_p->tid;
    int  easyID = (int ) arg_p->easyID;
    static int shared_var = 0;
    int local_var = 0;
    while (shared_var < 10){

      shared_var++;
      printf("ThreadID: %d, shared_var: %d\n",easyID, shared_var);
      //remove sleep to see behavior where single thread increments up to ten
      sleep(1);


    }
    printf("BARRIER Var Value: %d, Static Var Value: %d, Local Var Value: %d, Thread ID: %ld, easyID: %d \n",
    BARRIER, shared_var, ++local_var, *tid, easyID);
}

int main()
{

    // pthread_t tid0;
    // pthread_t tid1;
    // pthread_t tid2;
    // pthread_t *pthreads[] = {&tid0,&tid1,&tid2};


    for (int i = 0; i < N; i++) {
       ThreadArgs *args = (ThreadArgs *) malloc(sizeof(ThreadArgs));
       args->tid = pthreads[i];
       args->easyID = i;
	     pthread_create(args->tid, NULL, workerThreadFunc, (void *)args);

        };
        //Remove sleep to see thread behavior act sequentially
        sleep(1);
        BARRIER = 1;

    pthread_exit(NULL);
    return 0;
}
*/
