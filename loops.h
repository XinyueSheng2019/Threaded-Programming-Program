//B156924

#define N 729
#define reps 1000

double a[N][N], b[N][N], c[N];
int jmax[N];  


void init1(void);
void init2(void);
void runloop(int); 
void loop1chunk(int, int);
void loop2chunk(int, int);
void valid1(void);
void valid2(void);
void init_work(int* thread_start, int* thread_finish, int* remain_array, int* chunk_start, int* chunk_finish, int myid, int p, int each_sum, int* chunk);
void assign_work(int* thread_start, int* thread_finish, int* remain_array, int* chunk_start, int* chunk_finish, int* chunk, int myid, int p, int* complete);
int  steal_work(int myid,int *remain_array,int p, int* chunk_start, int* chunk_finish, int* thread_start, int* thread_finish, int* chunk);

