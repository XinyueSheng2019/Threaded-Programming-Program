// B156924

#include <stdio.h>
#include <math.h>
#include <omp.h> 
#include "loops.h"


int main(int argc, char *argv[]) { 

  double start1,start2,end1,end2;
  int r;

  init1(); 

  start1 = omp_get_wtime(); 

  for (r=0; r<reps; r++){ 
    runloop(1);
  } 

  end1  = omp_get_wtime();  

  valid1(); 

  printf("Total time for %d reps of loop 1 = %f\n",reps, (float)(end1-start1)); 


  init2(); 

  start2 = omp_get_wtime(); 

  for (r=0; r<reps; r++){ 
    runloop(2);
  } 

  end2  = omp_get_wtime(); 

  valid2(); 

  printf("Total time for %d reps of loop 2 = %f\n",reps, (float)(end2-start2)); 

} 

void init1(void)
{
  int i,j; 

  for (i=0; i<N; i++){ 
    for (j=0; j<N; j++){ 
      a[i][j] = 0.0; 
      b[i][j] = 3.142*(i+j); 
    }
  }

}

void init2(void)
{ 
  int i,j, expr; 

  for (i=0; i<N; i++){ 
    expr =  i%( 3*(i/30) + 1); 
    if ( expr == 0) { 
      jmax[i] = N;
    }
    else {
      jmax[i] = 1; 
    }
    c[i] = 0.0;
  }

  for (i=0; i<N; i++){ 
    for (j=0; j<N; j++){ 
      b[i][j] = (double) (i*j+1) / (double) (N*N); 
    }
  }
 
} 


void runloop(int loopid)  
{
  int p = omp_get_max_threads();

  int remain_array[p]; // create an array to store the number of remaining iterations
  int chunk_start[p];  //the start iteration of every chunk
  int chunk_finish[p]; //the finish iteration of every chunk
  int thread_start[p]; //the start iteration of every thread
  int thread_finish[p];//the finish iteration of every thread
  int chunk[p];
  int each_sum = (int) ceil((double)N/(double)p);  //the all number of iterations for each threads
  
  #pragma omp parallel default(none) shared(loopid, p, each_sum, chunk_start, chunk_finish, thread_start, thread_finish, remain_array, chunk) 
  {

    int myid = omp_get_thread_num();
    int complete = 0; //this is a flag, when complete==0, it means there are still some threads are executing their work, when complete==1, it means this thread has finished all its work.

    init_work(thread_start, thread_finish, remain_array, chunk_start, chunk_finish, myid, p, each_sum, chunk);
     
    #pragma omp barrier
    while(complete == 0)
    {
      switch (loopid) 
      { 
        case 1: loop1chunk(chunk_start[myid],chunk_finish[myid]); break;
        case 2: loop2chunk(chunk_start[myid],chunk_finish[myid]); break;
      } 

      #pragma omp critical
      {
        assign_work(thread_start, thread_finish, remain_array, chunk_start, chunk_finish, chunk, myid, p, &complete);
      }
    }


  }
  
}

void init_work(int* thread_start, int* thread_finish, int* remain_array, int* chunk_start, int* chunk_finish, int myid, int p, int each_sum, int* chunk)
{
    thread_start[myid] = myid * each_sum;
    thread_finish[myid] = ((myid+1) * each_sum) > N ? N:(myid+1) * each_sum;
    remain_array[myid] = thread_finish[myid] - thread_start[myid];
    chunk[myid]  = (int) ceil(remain_array[myid]/(double)p);  //initial the first chunk
    remain_array[myid] =  remain_array[myid] - chunk[myid];   //initial the remaining iterations in each thread
    chunk_start[myid] = thread_start[myid];                   //initial each chunk's start iteration and finish iteration
    chunk_finish[myid] = chunk_start[myid] + chunk[myid];
    
}
void assign_work(int* thread_start, int* thread_finish, int* remain_array, int* chunk_start, int* chunk_finish, int* chunk, int myid, int p, int* complete)
{
      chunk_start[myid] += chunk[myid];                       //move to the next chunk
      chunk[myid] = (int) ceil(remain_array[myid]/(double)p); //calculate new chunk
      chunk_finish[myid] = (chunk_start[myid] + chunk[myid]) > thread_finish[myid]? thread_finish[myid]:(chunk_start[myid] + chunk[myid]); //judge whether chunk finish is bigger than thread finish iteration.
     
      if (chunk_finish[myid] - chunk_start[myid] > 0)         //there is remaining iterations in this thread
      {
        chunk[myid] = chunk_finish[myid] - chunk_start[myid]; //refresh the chunk
        remain_array[myid] = remain_array[myid] - chunk[myid];
      } 
      else //there is no remaining iteration in this thread. Now we can consider how to steal work from other threads!
      {  
        *complete = steal_work(myid, remain_array, p, chunk_start, chunk_finish, thread_start, thread_finish, chunk);  
      } 
}

int steal_work(int myid,int *remain_array,int p, int* chunk_start, int* chunk_finish, int* thread_start, int* thread_finish, int* chunk)
{
  int id = 0; 
  int max_remain = 0; //to find the maxinum remaining iterations in all threads
  int find = myid;    //the thread id with the maxinum remaining iterations 
  int complete = 1;   //a flag to judge whether there are still some thread with remaining iterations.
                      // Firstly we assume that there is no available iterations to be stolen, so the flag sets 1.
  int steal_chunk;    //the chunk stolen from others.

  for(id = 0;id < p; id++) //look for the thread with the maxinum remaining iterations 
  {
    if (remain_array[id] > max_remain) 
    {
      max_remain = remain_array[id];
      find = id;
    }
  }

  if (find != myid) //have found the thread with the maximal remaining iterations 
  {
    complete = 0;
    steal_chunk = (int) ceil((double)max_remain/(double)p); //the work stolen from the thread with the maxinum remaining iterations 

    thread_start[myid] = chunk_start[find] + chunk[find];   //the thread which obtains the work should change its work boundaries.
    thread_finish[myid] = thread_start[myid] + steal_chunk;
    remain_array[find] -= steal_chunk;                      //the thread which is stolen should reduce the stolen work
    chunk[find] += steal_chunk; //add stolen work to the chunk which is being executed by the stolen thread, therefore, next time its chunk will start at the end iteration of the stolen work
    
    chunk[myid] =  (int) ceil((double)steal_chunk/(double)p); //create the new chunk for the thread which obtains the work
    chunk_start[myid] = thread_start[myid];                   //the start iteration of the chunk should be that of the thread. In other words, the thread is initialized again.
    chunk_finish[myid] = (chunk_start[myid] + chunk[myid])>thread_finish[myid]? thread_finish[myid]:(chunk_start[myid] + chunk[myid]); ////judge whether chunk finish is bigger than thread finish iteration.
    remain_array[myid] += steal_chunk-chunk[myid]; 
      
  }
  return complete; 
}


void loop1chunk(int lo, int hi) 
{ 
  int i,j; 
  for (i=lo; i<hi; i++)
  { 
    for (j=N-1; j>i; j--)
    {
      a[i][j] += cos(b[i][j]);
    } 
  }
  

} 



void loop2chunk(int lo, int hi) 
{
  int i,j,k; 
  double rN2; 
  rN2 = 1.0 / (double) (N*N);  

  for (i=lo; i<hi; i++)
  { 
    for (j=0; j < jmax[i]; j++)
    {
      for (k=0; k<j; k++)
      { 
       c[i] += (k+1) * log (b[i][j]) * rN2;
      } 
    }
  }

}

void valid1(void) 
{ 
  int i,j; 
  double suma; 
  
  suma= 0.0; 
  for (i=0; i<N; i++){ 
    for (j=0; j<N; j++){ 
      suma += a[i][j];
    }
  }
  printf("Loop 1 check: Sum of a is %lf\n", suma);

} 


void valid2(void) 
{ 
  int i; 
  double sumc; 
  
  sumc= 0.0; 
  for (i=0; i<N; i++){ 
    sumc += c[i];
  }
  printf("Loop 2 check: Sum of c is %f\n", sumc);
} 
 
