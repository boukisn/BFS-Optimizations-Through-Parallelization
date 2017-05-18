// gcc -O1 -o BFS_Serial_Opt_Half BFS_Serial_Opt_Half.c -lrt -lm

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <pthread.h>
 
#define N 16384
//#define N 32768
int visited[N];
int in[N];
int out[N];
int parents[N];
int where_in[N];
int in_place = 0;
int out_place = 0;
int visit_count = 0;

#define INFTY 1000
#define ALPHA 4
#define BETA  10
#define OPTIONS 5
#define GIG 1000000000
#define CPG 2.9 
#define NUM_THREADS 4
#define BLOCK_SIZE 2048

#define ITERS 25

const char *option_desc[OPTIONS];


int OPTION;
struct timespec diff(struct timespec start, struct timespec end);
struct timespec time1, time2;
struct timespec time_stamp[OPTIONS][ITERS+1];;

struct thread_data{
  int thread_id;
  int v;
};

pthread_mutex_t lock;
pthread_barrier_t barr;

 
char * graph;
 
// graph as adjacency matrix
//0 1 2 3 4 5 6
/*
int graph[N * N] = {  0, 0, 1, 1, 0, 0, 0, 0, 0, 0,//0
        0, 0, 0, 0, 1, 1, 0, 0, 0, 0, //1
        1, 0, 0, 1, 1, 0, 0, 1, 0, 0, //2
        1, 0, 1, 0, 0, 0, 0, 0, 0, 0, //3
        0, 1, 1, 0, 0, 1, 0, 0, 1, 0, //4
        0, 1, 0, 0, 1, 0, 1, 0, 0, 1, //5
        0, 0, 0, 0, 0, 1, 0, 0, 0, 1,
        0, 0, 1, 0, 0, 0, 0, 0, 1, 0,
        0, 0, 0, 0, 1, 0, 0, 1, 0, 1,
        0, 0, 0, 0, 0, 1, 1, 0, 1, 0,
        };
 */
// queue implementation
int front;
int rear;
int q[N];
int check_count;



void bfs(int v);
void synch_bfs(int v);
void initialize_graph(int w);
int check_graph();
int check_visited();
void print_graph();
void print_row(int row);
void Top2down();
void Down2top();
void *single_thread_bfs(void *threadarg);
void pthread_bfs();
void *single_thread_synch_bfs(void *threadarg);
void pthread_synch_bfs();
void *single_thread_synch_blocked_bfs(void *threadarg);
void pthread_synch_blocked_bfs();
 
int main() {
    //used to create easily readable outputs
    option_desc[0] = "1. Simple BFS...................";
    option_desc[1] = "2. Parallel Simple BFS..........";
    option_desc[2] = "3. Synch BFS....................";
    option_desc[3] = "4. Parallel Synch BFS...........";
	option_desc[4] = "5. Parallel Blocked Synch BFS...";

	printf("\n\n----------------------------------------------------------------\n");
	printf("INFORMATION \n");
	printf("----------------------------------------------------------------\n");
	printf("Number of nodes:           %d\n", N);
	//printf("Size of adjacency matrix:  %d bytes\n", N*N);
	printf("Number of threads:         %d\n", NUM_THREADS);
	printf("Block size:                %d\n", BLOCK_SIZE);
	printf("Alpha value:               %d\n", ALPHA);
	printf("Beta value:                %d\n", BETA);
	printf("Iterations:                %d\n", ITERS);
	printf("Options:                   %d\n", OPTIONS);
	

	unsigned long long total_size = ((((unsigned long long)(N))*(((unsigned long long)(N))-1))/2);

    // make all vertex unvisited    
    graph = (char *) malloc(total_size);
	
	printf("\n\n----------------------------------------------------------------\n");
	printf("TESTING \n");
	printf("----------------------------------------------------------------\n");

    initialize_graph(50);
    if(check_graph())
        printf("Graph is connected!\n");
	else
	{
		printf("Graph is not connected... exiting.\n\n");
		return 0;
	}

    //for each option clear memory and time results
    int i, j;
    check_count = 0;
    OPTION = 0;
    for (i = 0; i < ITERS; i++) 
    {
        front = 0;
        rear = 0;
        memset(visited, 0, sizeof(visited));
        memset(parents, 0, sizeof(parents));
        memset(q, 0, sizeof(q));
        clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &time1);
        bfs(0);
        clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &time2);
        time_stamp[OPTION][i] = diff(time1,time2);
        if(check_visited())
            check_count++;
    }

    //checks to see if each algorithm succusfully tranversed the graph 
    if(check_count == ITERS)
        printf("%s All tests passed!\n", option_desc[OPTION]);
    else
         printf("%s Only %d/%d tests passed.\n", option_desc[OPTION], check_count, ITERS);


    check_count = 0;
    OPTION++;
    for (i = 0; i < ITERS; i++) 
    {
        front = 0;
        rear = 0;
        memset(visited, 0, sizeof(visited));
        memset(parents, 0, sizeof(parents));
        memset(q, 0, sizeof(q));
        clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &time1);
        pthread_bfs();
        clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &time2);
        time_stamp[OPTION][i] = diff(time1,time2);
        if(check_visited())
            check_count++;
    }

    if(check_count == ITERS)
        printf("%s All tests passed!\n", option_desc[OPTION]);
    else
         printf("%s Only %d/%d tests passed.\n", option_desc[OPTION], check_count, ITERS);


    check_count = 0;
    OPTION++;
    for (i = 0; i < ITERS; i++) 
    {
        memset(visited, 0, sizeof(visited));
        memset(parents, 0, sizeof(parents));
        memset(where_in, 0, sizeof(where_in));
        memset(q, 0, sizeof(q));
        clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &time1);
        synch_bfs(0);
        clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &time2);
        time_stamp[OPTION][i] = diff(time1,time2);
        if(check_visited())
            check_count++;
    }

    if(check_count == ITERS)
        printf("%s All tests passed!\n", option_desc[OPTION]);
    else
         printf("%s Only %d/%d tests passed.\n", option_desc[OPTION], check_count, ITERS);


   	check_count = 0;
    OPTION++;
    for (i = 0; i < ITERS; i++) 
    {
        memset(visited, 0, sizeof(visited));
        memset(parents, 0, sizeof(parents));
        memset(where_in, 0, sizeof(where_in));
        memset(q, 0, sizeof(q));
        clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &time1);
        pthread_synch_bfs();
        clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &time2);
        time_stamp[OPTION][i] = diff(time1,time2);
        if(check_visited())
            check_count++;
    }

    if(check_count == ITERS)
        printf("%s All tests passed!\n", option_desc[OPTION]);
    else
         printf("%s Only %d/%d tests passed.\n", option_desc[OPTION], check_count, ITERS);

	check_count = 0;
    OPTION++;
    for (i = 0; i < ITERS; i++) 
    {
        memset(visited, 0, sizeof(visited));
        memset(parents, 0, sizeof(parents));
        memset(where_in, 0, sizeof(where_in));
        memset(q, 0, sizeof(q));
        clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &time1);
        pthread_synch_blocked_bfs();
        clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &time2);
        time_stamp[OPTION][i] = diff(time1,time2);
        if(check_visited())
            check_count++;
    }

    if(check_count == ITERS)
        printf("%s All tests passed!\n", option_desc[OPTION]);
    else
         printf("%s Only %d/%d tests passed.\n", option_desc[OPTION], check_count, ITERS);
	
	printf("\n\n----------------------------------------------------------------\n");
    printf("TIMING \n");
	printf("----------------------------------------------------------------\n");

    float sum;
    float average;
    float total = (float)(ITERS);
    // finds the average time for each of the options and prints them out
    for (j = 0; j < OPTIONS; j++) 
    {
        sum = 0;
        average = 0;
        printf("%s ", option_desc[j]);
        for (i = 0; i < ITERS; i++) 
        {
            float val = (float)((double)(CPG)*(double)(GIG * time_stamp[j][i].tv_sec + time_stamp[j][i].tv_nsec));
            sum += val;
        }
        average = (sum/total);
        printf("%0.0f cycles per element\n", ((average)/((float)(N))));
    }
	printf("\n\n");
    return 0;
}
 
void bfs(int v) {
 
    // make vertex v visited
    visited[v] = 1;
    parents[v] = v;
    // enqueue vertex v
    q[rear] = v; // insert at rear
    rear++; // increment rear
 
    while (rear != front) // condition for empty queue
    {
        // dequeue
        int u = q[front];
        //printf("%d ", u);
        front++;
 
        // check adjacent nodes from u
        int j = 0;
        for (j = 0; j < N; j++) {
            // if there is adjacent vertex enqueue it
            if (!visited[j] && get(u,j)) {
                q[rear] = j;
                rear++;
                visited[j] = 1;
                parents[j] = u;
            }
        }
    }
}

void synch_bfs(int v)
{
    //visits the root node and makes it the inital frontier 
    visited[v] = 1;
    visit_count=1;
    in[0] = v;
    in_place++;
    where_in[0] = 1;
    parents[v]=v;
    int i;
    char state = 'A';

    while(in_place != 0) //loop until there is no longer a frontier 
    {

        out_place =0;
        int case1 = (N - visit_count)/ALPHA; // determine the cases for state changes 
        int case2 = N/BETA;
        switch(state) //based off of given state 
        {
            case 'A':
            {
                Top2down();
                if(in_place >= case1)
                    state= 'B';
                break;
            }
            case 'B':
            {
                Down2top();
                if (in_place <= case2)
                    state ='C';
                break;
            }
            case 'C':
            {
                Top2down();
                break;
            }

        }
        ///given step
        memset(where_in, 0, sizeof(where_in)); // clear where_in aray
        in_place = out_place; // make in list the size of out list
        for (i = 0; i < out_place; ++i) // transfer data from out array to in array 
        {
            in[i] = out[i];
            where_in[out[i]] = 1; // if a given element is in in mark it as high 
        }
    }

}

void Top2down()
{
    int vertex, neigh;
    for (vertex = 0; vertex < in_place; ++vertex) // loops through current frontier
    {
        int current = in[vertex]; // gets exact node for makeshift list
        for (neigh = 0; neigh < N; ++neigh) // loops through all nodes
        {
            if (!visited[neigh] && get(current, neigh))  // checks if a connection exist and the node has not been visited 
            {
                parents[neigh] = current; //updates parent away
                out[out_place] = neigh; // updates the out list
                out_place++; // keeps track of length of out list
                visited[neigh] = 1; // updates visited array 
                visit_count++; // keeps track of the amount of nodes visited for Alpha and beta
            }
        }
    }

}

void Down2top()
{
    int node, neigh;
    for ( node = 1; node < N; ++node) // loops through all nodes
    {
        if (!visited[node] ) //checks if it had been visited 
        {
            for (neigh = 0; neigh < node; ++neigh) // looks for all other nodes for parents && only need to check values less than due to halved memory structure
            {  
    		    if (graph[((node*(node-1))/2)+neigh] && !where_in[neigh] ) // checks if there is a connection and parent is not in the current frontier
    		    {
    		        parents[node] = neigh; //updates parent away
    		        out[out_place] = neigh; // updates the out list
    		        out_place++; // keeps track of length of out list
    		        visited[node] = 1; // updates visited array 
    		        visit_count++; // keeps track of the amount of nodes visited for Alpha and beta
    		        break; // breaks if parent is found
    		    }
            }
        }
    }
}

void *single_thread_synch_bfs(void *threadarg)
{
      struct thread_data *my_data;
      my_data = (struct thread_data *) threadarg;
      int pid = my_data->thread_id;
      int v = my_data->v;

    int i;
	 int vertex, neigh;
	char state = 'A';
	if(pid == 0) //the first thread visits the root and its neighbors 
    {
        visited[0] = 1;
        visit_count=1;
        in[0] = 0;
        in_place++;
        where_in[0] = 1;
        parents[0]=0;
        out_place =0;
        
        for (vertex = 0; vertex < in_place; ++vertex)
        {
            int current = in[vertex];
            for (neigh = 0; neigh < N; ++neigh)
            {
                if (!visited[neigh] && get(current,neigh))
                {
                    parents[neigh] = current;
                    visited[neigh] = 1;
                    out[out_place] = neigh;
                    out_place++;
                    visit_count++;
                }
            }
        }
    }
    int rc1 = pthread_barrier_wait(&barr); // all threads wait before going to the next frontier 
    if (rc1 != 0 && rc1 != PTHREAD_BARRIER_SERIAL_THREAD) {
      printf("Could not wait on barrier\n");
      exit(-1);
    }


	int case1, case2;
    while(in_place != 0)
    {
        
        switch(state)
        {

            case 'A':
            {
                int vertex, neigh;
                for (vertex = pid; vertex < in_place; vertex+=NUM_THREADS) // each threads attempts to find a node to analyze
				{
				    int current = in[vertex];
				    for (neigh = 0; neigh < N; ++neigh)
				    {
						if (!visited[neigh] && get(current,neigh)) 
	                    {
							if (pthread_mutex_lock(&lock)) printf("\nERROR on lock\n"); //locks whats checked between threads
	                        parents[neigh] = current;
	                        visited[neigh] = 1;
	                        out[out_place] = neigh;
	                        out_place++;
	                        visit_count++;
							if (pthread_mutex_unlock(&lock)) printf("\nERROR on unlock\n");
	                    }
                    }
                }
                if(in_place >= case1)
                    state= 'B';
                break;
            }
            case 'B':
            {

                    int node, neigh;
					int start = (v == 0) ? 1 : v;
                    for (node = start; node < (v + (int)(N/NUM_THREADS)); ++node) // each thread gets a range of nodes to analyze better for fixed size
                    {
                        if (!visited[node] )
                        {
                            for (neigh = 0; neigh < node; neigh += 1)
                            {
                                if (get(node,neigh) && !where_in[neigh]) 
                                {
									if (pthread_mutex_lock(&lock)) printf("\nERROR on lock\n");//locks whats checked between threads
                                    parents[node] = neigh;
                                    visited[node] = 1;
                                    out[out_place] = neigh;
                                    out_place++;
                                    visit_count++;
									if (pthread_mutex_unlock(&lock)) printf("\nERROR on unlock\n");
                                    break;
									
                                }

                            }
                        }
                    }

                if (in_place <= case2)
                    state ='C';
                break;
            }
            case 'C':
            {
                int vertex, neigh;
                for (vertex = pid; vertex < in_place; vertex+=NUM_THREADS) // each threads attempts to find a node to analyze
				{
				    int current = in[vertex];
				    for (neigh = 0; neigh < N; ++neigh)
				    {
						if (!visited[neigh] && get(current,neigh)) 
	                    {
	                        if (pthread_mutex_lock(&lock)) printf("\nERROR on lock\n"); //locks whats checked between threads
							parents[neigh] = current;
	                        visited[neigh] = 1;
	                        out[out_place] = neigh;
	                        out_place++;
	                        visit_count++;
	                        if (pthread_mutex_unlock(&lock)) printf("\nERROR on unlock\n");
	                    }
                    }
                }
                break;
            }

        }

		int rc3 = pthread_barrier_wait(&barr); // wait before updating the next frontier 
        if (rc3 != 0 && rc3 != PTHREAD_BARRIER_SERIAL_THREAD) {
          printf("Could not wait on barrier\n");
          exit(-1);
        }

        ///the first thread handles data tranfer and calculating state transitions 
        if(pid == 0){
            memset(where_in, 0, sizeof(where_in));
            in_place = out_place;
            for (i = 0; i < out_place; ++i)
            {
                in[i] = out[i];
                where_in[out[i]] = 1;
            }
			out_place = 0;
			case1 = (N - visit_count)/ALPHA;
			case2 = N/BETA;
        }

		int rc2 = pthread_barrier_wait(&barr); // wait before moving on
		if (rc2 != 0 && rc2 != PTHREAD_BARRIER_SERIAL_THREAD) {
		  printf("Could not wait on barrier\n");
		  exit(-1);
		}
    }

}

void pthread_synch_bfs() { //initilizes threads and components needed for them 
 
  pthread_t threads[NUM_THREADS];
  struct thread_data thread_data_array[NUM_THREADS];
  int rc;
  long t;

  if(pthread_barrier_init(&barr, NULL, NUM_THREADS)) {
    printf("Could not create a barrier\n");
  } 


  if (pthread_mutex_init(&lock, NULL)) {
    printf("Could not create a lock\n");
  }

  //printf("-------------------\n");

  for (t = 0; t < NUM_THREADS; t++) {
    thread_data_array[t].thread_id = t;
    thread_data_array[t].v = ((int)(N/NUM_THREADS))*t;
    rc = pthread_create(&threads[t], NULL, single_thread_synch_bfs, (void*) &thread_data_array[t]);
    if (rc) {
      printf("ERROR; return code from pthread_create() is %d\n", rc);
      exit(-1);
    }
  }

  for (t = 0; t < NUM_THREADS; t++) {
    if (pthread_join(threads[t],NULL)){ 
      printf("ERROR; code on return from join is %d\n", rc);
      exit(-1);
    }
  }
}

void *single_thread_synch_blocked_bfs(void *threadarg) // same as previous function but with blocking for case B
{
      struct thread_data *my_data;
      my_data = (struct thread_data *) threadarg;
      int pid = my_data->thread_id;
      int v = my_data->v;

    visited[0] = 1;
    visit_count=1;
    in[0] = 0;
    in_place++;
    where_in[0] = 1;
    parents[0]=0;
    int i;
    char state = 'A';

	int case1, case2;
    while(in_place != 0)
    {
        //printf("state %c\n",state);

        
        switch(state)
        {

            case 'A':
            {
                int vertex, neigh;
                for (vertex = pid; vertex < in_place; vertex+=NUM_THREADS)
				{
				    int current = in[vertex];
				    for (neigh = 0; neigh < N; ++neigh)
				    {
						if (!visited[neigh] && get(current,neigh)) 
	                    {
	                        parents[neigh] = current;
	                        visited[neigh] = 1;
	                        out[out_place] = neigh;
	                        out_place++;
	                        visit_count++;
	                    }
                    }
                }
                if(in_place >= case1)
                    state= 'B';
                break;
            }
            case 'B':
            {

				int block_i, block_j;
				for(block_i = pid; block_i < N/BLOCK_SIZE; block_i+=NUM_THREADS) // Case B has a fixed size loop, thus only place we can put blocking 
				{
					for(block_j = 0; block_j <= block_i; block_j++)
					{
						if(block_i != block_j)
						{
							int node, neigh;
							for (node = (block_i*BLOCK_SIZE); node < ((block_i*BLOCK_SIZE) + BLOCK_SIZE); ++node)
							{
							    if (!visited[node] )
							    {
							        for (neigh = (block_j*BLOCK_SIZE); neigh < ((block_j*BLOCK_SIZE) + BLOCK_SIZE); ++neigh)
							        {
							            if (graph[((node*(node-1))/2)+neigh] && !where_in[neigh]) // need to find good way to see if in
							            {
							                parents[node] = neigh;
							                visited[node] = 1;
							                out[out_place] = neigh;
							                out_place++;
							                visit_count++;
							                break;
							            }

							        }
							    }
							}
						}
						else
						{
							int node, neigh;
							for (node = (block_i*BLOCK_SIZE) + 1; node < ((block_i*BLOCK_SIZE) + BLOCK_SIZE); ++node)
							{
							    if (!visited[node] )
							    {
							        for (neigh = (block_j*BLOCK_SIZE); neigh < node; ++neigh)
							        {
							            if (graph[((node*(node-1))/2)+neigh] && !where_in[neigh]) // need to find good way to see if in
							            {
							                parents[node] = neigh;
							                visited[node] = 1;
							                out[out_place] = neigh;
							                out_place++;
							                visit_count++;
							                break;
							            }

							        }
							    }
							}
						}
					}
				}
                
                if (in_place <= case2)
                    state ='C';
                break;
            }
            case 'C':
            {
                int vertex, neigh;
                for (vertex = pid; vertex < in_place; vertex+=NUM_THREADS)
				{
				    int current = in[vertex];
				    for (neigh = 0; neigh < N; ++neigh)
				    {
						if (!visited[neigh] && get(current,neigh)) 
	                    {
	                        parents[neigh] = current;
	                        visited[neigh] = 1;
	                        out[out_place] = neigh;
	                        out_place++;
	                        visit_count++;
	                    }
                    }
                }
                break;
            }

        }

		int rc1 = pthread_barrier_wait(&barr);
		if (rc1 != 0 && rc1 != PTHREAD_BARRIER_SERIAL_THREAD) {
		  printf("Could not wait on barrier\n");
		  exit(-1);
		}

        ///given step
        if(pid == 0){
            memset(where_in, 0, sizeof(where_in));
            in_place = out_place;
            for (i = 0; i < out_place; ++i)
            {
                in[i] = out[i];
                where_in[out[i]] = 1;
            }
			out_place = 0;
			case1 = (N - visit_count)/ALPHA;
			case2 = N/BETA;
        }

		int rc2 = pthread_barrier_wait(&barr);
		if (rc2 != 0 && rc2 != PTHREAD_BARRIER_SERIAL_THREAD) {
		  printf("Could not wait on barrier\n");
		  exit(-1);
		}
    }

}

void pthread_synch_blocked_bfs() {  //initilizes threads and components needed for them 
 
  pthread_t threads[NUM_THREADS];
  struct thread_data thread_data_array[NUM_THREADS];
  int rc;
  long t;

  if(pthread_barrier_init(&barr, NULL, NUM_THREADS)) {
    printf("Could not create a barrier\n");
  } 


  if (pthread_mutex_init(&lock, NULL)) {
    printf("Could not create a lock\n");
  }

  //printf("-------------------\n");

  for (t = 0; t < NUM_THREADS; t++) {
    thread_data_array[t].thread_id = t;
    thread_data_array[t].v = ((int)(N/NUM_THREADS))*t;
    rc = pthread_create(&threads[t], NULL, single_thread_synch_blocked_bfs, (void*) &thread_data_array[t]);
    if (rc) {
      printf("ERROR; return code from pthread_create() is %d\n", rc);
      exit(-1);
    }
  }

  for (t = 0; t < NUM_THREADS; t++) {
    if (pthread_join(threads[t],NULL)){ 
      printf("ERROR; code on return from join is %d\n", rc);
      exit(-1);
    }
  }
}

void *single_thread_bfs(void *threadarg)
{
  struct thread_data *my_data;
  my_data = (struct thread_data *) threadarg;
  int pid = my_data->thread_id;
  int v = my_data->v;
    // make vertex v visited
    visited[0] = 1;
    //printf("THREAD (%d) visited NODE (%d)\n", pid,v);
    parents[0] = 0;
    // enqueue vertex v
    if(pid == 0) // the first thread checks the given root
	{
		q[rear] = 0; // insert at rear
		rear++; // increment rear
    }
	
	int rc1 = pthread_barrier_wait(&barr);
		if (rc1 != 0 && rc1 != PTHREAD_BARRIER_SERIAL_THREAD) {
		  printf("Could not wait on barrier\n");
		  exit(-1);
		}

	int myfront = front;
    while (rear != myfront) // condition for empty queue
    {
        // dequeue
        int u = q[myfront];
		myfront++;

        // check adjacent nodes from u
        int i;
        for (i = v; i < (v + (int)(N/NUM_THREADS)); ++i) //each thread gets a range of neighbors they should check 
		{
            // if there is adjacent vertex enqueue it
            if (!visited[i] && get(u, i)) {
                if (pthread_mutex_lock(&lock)) printf("\nERROR on lock\n");
                q[rear] = i;
                rear++;
                if (pthread_mutex_unlock(&lock)) printf("\nERROR on unlock\n");
				visited[i] = 1;
                parents[i] = u;
            }
        }
    }

  pthread_exit(NULL);
}


void pthread_bfs() {  //initilizes threads and components needed for them 
 
  pthread_t threads[NUM_THREADS];
  struct thread_data thread_data_array[NUM_THREADS];
  int rc;
  long t;

  if(pthread_barrier_init(&barr, NULL, NUM_THREADS)) {
    printf("Could not create a barrier\n");
  } 


  if (pthread_mutex_init(&lock, NULL)) {
    printf("Could not create a lock\n");
  }

  //printf("-------------------\n");

  for (t = 0; t < NUM_THREADS; t++) {
    thread_data_array[t].thread_id = t;
    thread_data_array[t].v = ((int)(N/NUM_THREADS))*t;
    rc = pthread_create(&threads[t], NULL, single_thread_bfs, (void*) &thread_data_array[t]);
    if (rc) {
      printf("ERROR; return code from pthread_create() is %d\n", rc);
      exit(-1);
    }
  }

  for (t = 0; t < NUM_THREADS; t++) {
    if (pthread_join(threads[t],NULL)){ 
      printf("ERROR; code on return from join is %d\n", rc);
      exit(-1);
    }
  }

  //printf("-------------------\n\n");
}

void initialize_graph(int w)// initializes graph with 1's and 0's with the given percentage w
{
    int i, j;
    srand(time(0));
    for(i = 1; i < N; i++)
    {
        for(j = 0; j < i; j++)
        {
                int dist = rand() % 100; // 0 -> 100
        int r = (dist < w) ? 1 : 0;
                graph[((i*(i-1))/2)+j] = r;
        }
    }
}

int check_graph() // checks to see if all nodes can be visited 
{
    int i, j;
    int counter = 0;
    int non_zero_column = 1; 
    for(i = 1; i < N; i++)
    {
        for(j = 0; j < i; j++)
        {
            if(graph[((i*(i-1))/2)+j] == 0)
                counter++;
            else
                break;
        }

        if(counter == N)
            return 0;

        counter = 0;
    }

    return 1;
}

int get(int i, int j) // makes sure that given two coordinates that our graph is poperly indexed 
{
	int val;

	if(i == j)
		return 0;
	else
		val = (i > j) ? graph[((i*(i-1))/2)+j] : graph[((j*(j-1))/2)+i];

	return val;
}

int check_visited() // checks to see if all nodes are visited 
{
    int i;
    for(i = 0; i < N; i++)
    {
        if(visited[i] == 0)
            return 0;
    }

    return 1;
}


struct timespec diff(struct timespec start, struct timespec end)
{
  struct timespec temp;
  if ((end.tv_nsec-start.tv_nsec)<0) {
    temp.tv_sec = end.tv_sec-start.tv_sec-1;
    temp.tv_nsec = 1000000000+end.tv_nsec-start.tv_nsec;
  } else {
    temp.tv_sec = end.tv_sec-start.tv_sec;
    temp.tv_nsec = end.tv_nsec-start.tv_nsec;
  }
  return temp;
}

