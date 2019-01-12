#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <mpi.h>
//time
struct timespec tstart, tover;
long trio = 0, tcomp = 0, twio = 0, ttemp = 0;
void calc_time(long* target, struct timespec a, struct timespec b)
{
    int sec = a.tv_sec - b.tv_sec;
    int nsec = a.tv_nsec - b.tv_nsec;
    *target = ((long) sec) * 1000000000 + nsec;
}
int i;
int main(int argc,char* argv[])
{
    int rc;
    rc = MPI_Init (&argc,&argv);
    if (rc != MPI_SUCCESS)
    {
        printf ("Error starting MPI program. Terminating.\n");
        fflush(stdout);
        MPI_Abort (MPI_COMM_WORLD, rc);
    }
    int myrank, numprocs;
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
    MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
    unsigned int N = atoi(argv[1]); //the size of the list n
    char* input_file = argv[2];//the input file name//n 32-bit floats
    char* output_file = argv[3];//the output file name
    float* buf = (float *) malloc (N * sizeof (float));
    
	clock_gettime(CLOCK_REALTIME, &tstart);
	MPI_File thefile,thefileout;
    ///MPI_File_open(MPI_Comm comm, char *filename,int amode, MPI_Info info, MPI_File *fh)
    MPI_File_open(MPI_COMM_WORLD,input_file,MPI_MODE_RDONLY, MPI_INFO_NULL, &thefile);
    MPI_Status status;
    ///MPI_File_read/write_all(MPI_File fh, void *buf, int count, MPI_Datatype datatype, MPI_Status *status)
    MPI_File_read_all(thefile,buf,N, MPI_FLOAT,&status);
    MPI_File_close(&thefile);
	clock_gettime(CLOCK_REALTIME, &tover);
	//calc_time(&trio, tover, tstart);	
	
	clock_gettime(CLOCK_REALTIME, &tstart);
    int task = N/numprocs;
	if(task%2!=0)
		task--;
    int r = N % numprocs;
    int idx = myrank*task ;
    int nextidx = task*(myrank+1);
    if(myrank == numprocs-1)
    {
        nextidx =  N;
    }
	clock_gettime(CLOCK_REALTIME, &tover);
	//calc_time(&ttemp, tover, tstart);
	//tcomp = ttemp +tcomp;
	
    int sorted = 0;
    while(!sorted)
    {
        sorted=1;        
        clock_gettime(CLOCK_REALTIME, &tstart);
		for(i = idx; i < nextidx; i += 2)
        {
            if(i+1!=nextidx)
            {
                if(buf[i] > buf[i+1])
                {
                    //swap(buf, i, i+1);
                    float temp=buf[i];
                    buf[i]=buf[i+1];
                    buf[i+1]=temp;
                    sorted = 0;
                }
            }
        }
		clock_gettime(CLOCK_REALTIME, &tover);
		//calc_time(&ttemp, tover, tstart);
		//tcomp = ttemp +tcomp;
		float last_num=buf[idx],next_num=buf[nextidx-1];
        if(numprocs>1)
        {
			MPI_Barrier(MPI_COMM_WORLD);
			if(myrank!=0 && myrank%2==0)
				MPI_Send(buf+idx, 1, MPI_FLOAT, myrank-1, myrank, MPI_COMM_WORLD);
            if(myrank+1<numprocs && myrank%2==1)
				MPI_Recv(&next_num, 1, MPI_FLOAT, myrank+1, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
            if(myrank!=0 && myrank%2==1)
				MPI_Send(buf+idx, 1, MPI_FLOAT, myrank-1, myrank, MPI_COMM_WORLD);
            if(myrank+1<numprocs && myrank%2==0)
				MPI_Recv(&next_num, 1, MPI_FLOAT, myrank+1, MPI_ANY_TAG, MPI_COMM_WORLD, &status);

            if(myrank+1<numprocs && myrank%2==0)
				MPI_Send(buf+nextidx-1, 1, MPI_FLOAT, myrank+1, myrank, MPI_COMM_WORLD);
            if(myrank!=0 &&  myrank%2==1)
				MPI_Recv(&last_num, 1, MPI_FLOAT, myrank-1, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
            if(myrank+1<numprocs && myrank%2==1)
				MPI_Send(buf+nextidx-1, 1, MPI_FLOAT, myrank+1, myrank, MPI_COMM_WORLD);
            if(myrank!=0 && myrank%2==0)
				MPI_Recv(&last_num, 1, MPI_FLOAT, myrank-1, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
			clock_gettime(CLOCK_REALTIME, &tstart);
			if(idx!=0)
			{
				if(last_num > buf[idx])
				{
					buf[idx]=last_num;
					sorted = 0;
				}
			}
			clock_gettime(CLOCK_REALTIME, &tover);
			//calc_time(&ttemp, tover, tstart);
			//tcomp = ttemp +tcomp;
        }
		clock_gettime(CLOCK_REALTIME, &tstart);
        for(i = idx+1; i < nextidx; i += 2)
        {
            if(i+1!=nextidx)
            {
                if(buf[i] > buf[i+1])
                {
                    //swap(buf, i, i+1);
                    float temp=buf[i];
                    buf[i]=buf[i+1];
                    buf[i+1]=temp;
                    sorted = 0;
                }
            }
            else if(i+1<N)
            {
                if(buf[i] > next_num)
                {
                    buf[i]=next_num;
                    sorted = 0;
                }
            }
        }
		clock_gettime(CLOCK_REALTIME, &tover);
		//calc_time(&ttemp, tover, tstart);
		//tcomp = ttemp +tcomp;
		MPI_Barrier(MPI_COMM_WORLD);
		int temp;
		//MPI_Allreduce(&sorted, &temp, 1, MPI_INT, MPI_PROD, MPI_COMM_WORLD);
		MPI_Allreduce(&sorted, &temp, 1, MPI_INT, MPI_LAND, MPI_COMM_WORLD);
		sorted = temp;
		MPI_Barrier(MPI_COMM_WORLD);
    }
	
	clock_gettime(CLOCK_REALTIME, &tstart);
    MPI_File_open(MPI_COMM_WORLD,output_file,MPI_MODE_CREATE | MPI_MODE_WRONLY, MPI_INFO_NULL, &thefileout);
	if(numprocs==1)
    {
        MPI_File_write(thefileout,buf,N, MPI_FLOAT,&status);
    }else
    {
        MPI_Offset offset = (MPI_Offset) idx*sizeof(float);
        MPI_File_seek(thefileout, offset, MPI_SEEK_SET);
        if(myrank==numprocs-1)
        {
            MPI_File_write(thefileout, buf+idx, N-idx, MPI_FLOAT, &status);
        }else if(task>0)
        {
            MPI_File_write(thefileout, buf+idx, task, MPI_FLOAT, &status);
        }
    }
	clock_gettime(CLOCK_REALTIME, &tover);
	//calc_time(&twio, tover, tstart);
	if(myrank == numprocs-1)
	{
		printf("IO time = %lf\n",((trio+twio)/1000000000.0));
		printf("computing time = %lf\n",(tcomp/1000000000.0));
	}
    MPI_File_close(&thefileout);
    free( buf );
    MPI_Finalize();
    return 0;
}
