#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <mpi.h>
#include <string.h>
float* merge_tmp;
float* buf;
int sorted = 0;
int rtc = 0;
//time
struct timespec tstart, tover;
long trio = 0, tcomp = 0, twio = 0, ttemp = 0;
void calc_time(long* target, struct timespec a, struct timespec b)
{
    int sec = a.tv_sec - b.tv_sec;
    int nsec = a.tv_nsec - b.tv_nsec;
    *target = ((long) sec) * 1000000000 + nsec;
}
int compare(const void *a, const void *b)//這函式是 qsort 所需的比較函式
{
    float c = *(float *)a;
    float d = *(float *)b;
    if(c < d)
		return -1;   //傳回 -1 代表 a < b
    else if (c == d)
		return 0;   //傳回   0 代表 a = b
    else 
		return 1;    //傳回  1 代表 a>b
}
void merge(float *a, int n, int len)
{
    int idx1 = 0;
    int idx2 = n;
    int i = 0;
    while(idx1<n || idx2 <len)
	{
        if(idx1<n && idx2<len)
		{
            if(a[idx1]<=a[idx2])
			{
                merge_tmp[i] = a[idx1];
                i++;
                idx1++;
            }else if(a[idx1]>a[idx2])
			{
                merge_tmp[i] = a[idx2];
                i++;
                idx2++;
            }
        }else if(idx2 <len){
            memcpy(merge_tmp+i, a+idx2, (len-idx2)*sizeof(float));
			break;
        }else if(idx1 <n){
			memcpy(merge_tmp+i, a+idx1, (n-idx1)*sizeof(float));
			break;
        }
    }
    memcpy(a, merge_tmp, len*sizeof(float));
}
int main(int argc,char* argv[])
{
    //初始化
    int rc = MPI_Init (&argc,&argv);
    if (rc != MPI_SUCCESS)
    {
        printf ("Error starting MPI program. Terminating.\n");
        fflush(stdout);
        MPI_Abort (MPI_COMM_WORLD, rc);
    }
    int myrank, numprocs;
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
    MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
    //load value
    unsigned int N = atoi(argv[1]); //the size of the list n
    //char* input_file = argv[2];//the input file name//n 32-bit floats
    //char* output_file = argv[3];//the output file name
      
    //計算 分配工作
	//clock_gettime(CLOCK_REALTIME, &tstart);
    int task = N/numprocs;
	int r = N%numprocs;
    int idx = myrank*task ;
    int nextidx = task*(myrank+1);
    if(myrank == numprocs-1)
		nextidx =  N;
	int mytask = nextidx-idx;
	int nexttask = (myrank == numprocs-2)?task+r:task;
	//clock_gettime(CLOCK_REALTIME, &tover);
	//calc_time(&ttemp, tover, tstart);
	//tcomp = ttemp;
	
	//取得動態記憶體
    buf = (float *) malloc ((mytask+nexttask) * sizeof (float));  
	merge_tmp = (float*)malloc((mytask+nexttask)*sizeof(float));
	//讀檔
    MPI_File thefile,thefileout;
	//clock_gettime(CLOCK_REALTIME, &tstart);
    MPI_File_open(MPI_COMM_WORLD,argv[2],MPI_MODE_RDONLY, MPI_INFO_NULL, &thefile);
    MPI_Status status;
	MPI_Offset offset = (MPI_Offset) idx*sizeof(float);
	MPI_File_seek(thefile, offset, MPI_SEEK_SET);
	MPI_File_read(thefile, buf, mytask, MPI_FLOAT, &status);
    MPI_File_close(&thefile); 
	//clock_gettime(CLOCK_REALTIME, &tover);
	//calc_time(&trio, tover, tstart);	
	
	//clock_gettime(CLOCK_REALTIME, &tstart);
	qsort(buf,mytask,sizeof(float),compare);
	//clock_gettime(CLOCK_REALTIME, &tover);
	//calc_time(&ttemp, tover, tstart);
	//tcomp = ttemp + tcomp;
    //odd even sort
	if(numprocs!=1 && task>0)
	{
		while(!sorted&&rtc<48)
		{
			int temp;
			rtc++;
			if(myrank<numprocs && myrank%2==1)
			{
				MPI_Send(buf, mytask, MPI_FLOAT, myrank-1, myrank, MPI_COMM_WORLD);
				MPI_Recv(buf, mytask, MPI_FLOAT, myrank-1, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
			}else if(myrank+1<numprocs && myrank%2==0)
			{
				MPI_Recv(buf+mytask, nexttask, MPI_FLOAT, myrank+1, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
				//clock_gettime(CLOCK_REALTIME, &tstart);
				if(buf[mytask-1]>buf[mytask])
				{
					//qsort(buf,mytask+nexttask,sizeof(float), compare);
					temp=0;
					merge(buf,mytask,mytask+nexttask);
				}
				//clock_gettime(CLOCK_REALTIME, &tover);
				//calc_time(&ttemp, tover, tstart);
				//tcomp = ttemp + tcomp;				
				MPI_Send(buf+mytask, nexttask, MPI_FLOAT, myrank+1, myrank, MPI_COMM_WORLD);
			}			
			if(myrank!=0 && myrank<numprocs && myrank%2==0)
			{
				MPI_Send(buf, mytask, MPI_FLOAT, myrank-1, myrank, MPI_COMM_WORLD);
				MPI_Recv(buf, mytask, MPI_FLOAT, myrank-1, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
			}else if(myrank+1<numprocs && myrank%2==1)
			{
				MPI_Recv(buf+mytask, nexttask , MPI_FLOAT, myrank+1, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
				//clock_gettime(CLOCK_REALTIME, &tstart);
				if(buf[mytask-1]>buf[mytask])
				{
					//qsort(buf,mytask+task,sizeof(float), compare);
					temp=0;
					merge(buf,mytask,mytask+nexttask);
				}
				//clock_gettime(CLOCK_REALTIME, &tover);
				//calc_time(&ttemp, tover, tstart);
				//tcomp = ttemp + tcomp;		
				MPI_Send(buf+mytask, nexttask, MPI_FLOAT, myrank+1, myrank, MPI_COMM_WORLD);
			}			
			//MPI_Allreduce(&sorted, &temp, 1, MPI_INT, MPI_PROD, MPI_COMM_WORLD);
			MPI_Allreduce(&temp, &sorted, 1, MPI_INT, MPI_LAND, MPI_COMM_WORLD);
		}
	}
	
    //寫檔
    //clock_gettime(CLOCK_REALTIME, &tstart);
	MPI_File_open(MPI_COMM_WORLD,argv[3],MPI_MODE_CREATE | MPI_MODE_WRONLY, MPI_INFO_NULL, &thefileout);
    if(numprocs==1)
    {
        MPI_File_write(thefileout,buf,N, MPI_FLOAT,&status);
    }else if(task==0&&myrank==numprocs-1)
	{
        MPI_File_write(thefileout,buf,N, MPI_FLOAT,&status);
    }else
    {
		MPI_File_seek(thefileout, offset, MPI_SEEK_SET);
		MPI_File_write(thefileout, buf, nextidx-idx, MPI_FLOAT, &status);
    }
    MPI_File_close(&thefileout);
	//clock_gettime(CLOCK_REALTIME, &tover);
	//calc_time(&twio, tover, tstart);
    //free( buf );
	/*if(myrank == numprocs-1)
	{
		printf("IO time = %lf\n",(trio+twio/1000000000.0));
		printf("computing time = %lf\n",(tcomp/1000000000.0));
	}*/
    MPI_Finalize();
    return 0;
}
