#include<iostream>
#include<math.h>
#include<cstring>
#include <omp.h>
#include "mpi.h"
#include "sort.h"

int size, rank; long ndata; long offset;
pSort::dataType* f_buf;

long part(pSort::dataType* rbuf, long low, long high);
void QUICk(pSort::dataType* rbuf, long l, long h);
void INSERTION(pSort::dataType* rbuf, long b, long c);
void QUICk3(pSort::dataType* rbuf, long l, long h);
void MERGe(pSort::dataType* rbuf, long len);
void merge(pSort::dataType* result, pSort::dataType* left, pSort::dataType* right, long leftLen, long rightLen);

void merge(pSort::dataType* rbuf, long l, long m, long r);
void mergeSort(pSort::dataType* rbuf, long n);
long min(long x, long y);


void QUICK33(long* recvd_indices, pSort::dataType* rbuf, long count); //count is the data wit each processor
void QUICK2(pSort::dataType* rbuf, long ll);  //leng is the 
void QUICKint(long* buff, long l, long h);
long partint(long* buf, long low, long high);


void pSort::init()
{
	MPI_Init(NULL, NULL);
	MPI_Comm_size(MPI_COMM_WORLD, &size);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
}

int pSort::sort(pSort::dataType** rbuf, long* l, pSort::SortType type)	//=> *l, l= &ndata;   
{

	int procs = 20;  
	if (type == 1)
	{
#pragma omp parallel num_threads(procs)
		{	
#pragma omp single 
			{
				QUICk(*rbuf, 0, (*l) - 1);
			}
		}
	}

	else if (type == 2)
	{
	
#pragma omp parallel num_threads(procs)
		{
#pragma omp single
			{
				QUICk3(*rbuf, 0, *l - 1);
			}

		}
			
	}
		
	
	else if (type == 3)
	{
		mergeSort(*rbuf, *l);
	}

	for (int i = 0; i < size; i++)
	{
		if (rank == i)
		{
			QUICK2(*rbuf, *l);
		}
	}
	
	*rbuf = f_buf;
	*l = offset;
	MPI_Barrier(MPI_COMM_WORLD);
	
	return 0;

}


void QUICKint(long* buff, long l, long h)
{
	if (l < h)
	{
		long p = partint(buff, l, h);

		QUICKint(buff, l, p - 1);
		QUICKint(buff, p + 1, h);

	}
}
long partint(long* buf, long low, long high)
{
	long pivot = buf[high];    // pivot 
	long i = (low - 1);  // Index of smaller element 

	for (long j = low; j <= (high - 1); j++)
	{
		// If current element is smaller than or 
		// equal to pivot 
		if (buf[j] <= pivot)
		{
			i++;
			long t = buf[i];
			buf[i] = buf[j];
			buf[j] = t;
		}
	}
	long t = buf[i + 1];
	buf[i + 1] = buf[high];
	buf[high] = t;
	return (i + 1);
}


void QUICK2(pSort::dataType* rbuf, long ll)  //leng is the 
{
	long* local_indices = new long[size];
	long leng = ll * size;
	for (int i = 0; i < size; i++)
	{
		if (rank == i)
		{

			for (int z = 0; z < size; z++)
			{
				int p = size * size;
				int y = floor(z * leng / p);
				local_indices[z] = rbuf[y].key;
			}
		}
	}
	long* pivots_recvd = new long[(size * size)]; //for P0 to receive from others
	pivots_recvd[0] = 0;
	long* recvd_indices = new long[(size - 1)];    //for recving from P0
	recvd_indices[0] = 0;
	if (rank == 0)
	{
		for (int i = 0; i < size; i++)
		{
			pivots_recvd[i] = local_indices[i];
			
		}
	}


	for (int i = 1; i < size; i++)
	{
		if (rank == i)		//processes send their pivots to P0
		{
			MPI_Send(local_indices, size, MPI_LONG, 0, 2, MPI_COMM_WORLD);

		}
	}

	if (rank == 0)
	{
		for (int i = 1; i < size; i++)		//P0 recvs pivots frm other processes
		{
			MPI_Recv(pivots_recvd + (i * size), size, MPI_LONG, i, 2, MPI_COMM_WORLD, MPI_STATUSES_IGNORE);
		}

		int  l = 0; int h = (size * size) - 1;

		if (l < h)		//P0 sorts recvd pivots
		{
			long p = partint(pivots_recvd, l, h);

			QUICKint(pivots_recvd, l, p - 1);
			QUICKint(pivots_recvd, p + 1, h);
		}


		long* pivots_sent = new long[(size - 1)];
		long t = (floor(size / 2)) - 1;

		for (int i = 1; i < size; i++)
		{
			int y = i * size + t;
			pivots_sent[i - 1] = pivots_recvd[y];
		}
		for (int i = 0; i < size - 1; i++)
		{
			std::memcpy(recvd_indices + i, pivots_sent + i, sizeof(recvd_indices));
		}

		for (int i = 1; i < size; i++)
		{
			MPI_Send(pivots_sent, size - 1, MPI_LONG, i, 3, MPI_COMM_WORLD);
		}
		
	}

	for (int i = 1; i < size; i++)
	{
		if (rank == i)
		{
			
			MPI_Recv(recvd_indices, size - 1, MPI_LONG, 0, 3, MPI_COMM_WORLD, MPI_STATUSES_IGNORE);
		}
	}

	MPI_Barrier(MPI_COMM_WORLD);
	
	delete[] local_indices;		//ones we sent P0
	delete[] pivots_recvd;

	QUICK33(recvd_indices, rbuf, ll);


}

void QUICK33(long* recvd_indices, pSort::dataType* rbuf, long count) //count is the data wit each processor
{
	
	MPI_Datatype Type;

	f_buf = new pSort::dataType[2 * count];   //final buffer, to hold partitions
	f_buf[0].key = 0;
	int blocklengs[2]; MPI_Aint indices[2]; MPI_Datatype old_types[2];
	blocklengs[0] = 1; blocklengs[1] = LOADSIZE;
	old_types[0] = MPI_LONG; old_types[1] = MPI_CHAR;
	MPI_Get_address(&((*f_buf).key), &indices[0]);
	MPI_Get_address(&((*f_buf).payload), &indices[1]);
	indices[1] = indices[1] - indices[0];
	indices[0] = 0;
	MPI_Type_create_struct(2, blocklengs, indices, old_types, &Type);
	MPI_Type_commit(&Type);

	//nu[0].key = 0; //nu[5].key = 5;  delete[] * rbuf;  *rbuf = nu;

	int add = 0;
	MPI_Status status;
	long send_count = 0;
	int rem = 0;


	for (int i = 0; i < size; i++)
	{
		if (rank == i)
		{
			long c = 0;
			long begin = 0;
			for (int j = 0; j < size; j++)	 //partition and send
			{
				send_count = 0;
				pSort::dataType* inter;
				begin = c;

				if ((i < rem && c > count) || (i >= rem && c > (count - 1)))	//last part not found, end of items
				{
					send_count = 0;
				}

				else if (j == size - 1)		//exhausted the pivots which are one less than the processes, made it to the last partitioin
				{
					if (i < rem)
					{
						send_count = count + 1 - c;

					}
					else
					{
						send_count = count - c;
					}
				}
				else
				{
					while (true)
					{
						if ((i < rem && c == count) || (i >= rem && c == (count - 1)))  //end of items while checking
						{
							break;
						}
						if (rbuf[c].key <= recvd_indices[j])	//c will determine the partition count
						{
							c++;
							send_count++;

						}
						else
						{
							break;
						}
					}
				}
				inter= new pSort::dataType[send_count];
				inter[0].key = 0;
				if (send_count == 0) { inter = nullptr; }
				for (long g = 0; g < send_count; g++)		//copy the partition to send fwd
				{

					inter[g].key = rbuf[begin + g].key;
					for (int ii = 0; ii < LOADSIZE; ii++)
					{
						inter[g].payload[ii] = rbuf[g+begin].payload[ii];
					}
				}

				if (rank == j)
				{
					if (inter != nullptr)
					{
						for (long r = 0; r < send_count; r++)   //to self
						{
							f_buf[r+offset].key = rbuf[begin + r].key;
							for (int ii = 0; ii < LOADSIZE; ii++)
							{
								f_buf[r+offset].payload[ii] = rbuf[r + begin].payload[ii];
							}

						}
						offset += send_count;
					}
				}

				else  //senders' rank != receivers' rank
				{
					MPI_Send(inter, send_count, Type, j, i, MPI_COMM_WORLD);
				}
			}

		}

		else if (rank != i)
		{

			MPI_Recv(f_buf + offset, (2 * count), Type, i, i, MPI_COMM_WORLD, &status);
			MPI_Get_count(&status, Type, &add);
			offset += add;

		}

	}
	MPI_Barrier(MPI_COMM_WORLD);

	for (int i = 0; i < size; i++)
	{
		if (rank == i)
		{
			QUICk(f_buf, 0, offset - 1);
		}
	}

	


long part(pSort::dataType* rbuf, long low, long high)
{
	long pivot = rbuf[high].key;    // pivot
	long i = (low - 1);  // Index of smaller element

	for (long j = low; j <= (high - 1); j++)
	{
		// If current element is smaller than or equal to pivot
		if (rbuf[j].key <= pivot)
		{
			i++;
			pSort::dataType t = rbuf[i];
			rbuf[i] = rbuf[j];
			rbuf[j] = t;
		}
	}
	pSort::dataType t = rbuf[i + 1];
	rbuf[i + 1] = rbuf[high];
	rbuf[high] = t;
	return (i + 1);
}

void QUICk(pSort::dataType* rbuf, long l, long h)		//works fine
{
	if (l < h)
	{
		
			long p = part(rbuf, l, h);

#pragma omp task untied
			{
				QUICk(rbuf, l, p - 1);
			}

#pragma omp task untied
			{
				QUICk(rbuf, p + 1, h);
			}
		
	}
	
}

void QUICk3(pSort::dataType* rbuf, long l, long h)		//works slow and incorrect
{
	if (l < h)
	{
		if (h - l <= 10)
		{
			INSERTION(rbuf, l, h + 1);
		}
		else
		{
			long p = part(rbuf, l, h);

#pragma omp task untied
			{
				QUICk(rbuf, l, p - 1);
			}

#pragma omp task untied
			{
				QUICk(rbuf, p + 1, h);
			}
		}
	}

}

void INSERTION(pSort::dataType* rbuf, long b, long c)
{
	for (long i = b+1; i < c; i++)
	{
		pSort::dataType tmp = rbuf[i];  long j;

		for (j = i; j >= b+1 && tmp.key < rbuf[j - 1].key; j--)
			rbuf[j] = rbuf[j - 1];

		rbuf[j] = tmp;
	}
}



void pSort::close()
{
	MPI_Finalize();
	return;
}

	
pSort::dataType* slice(pSort::dataType* rbuf, long start, long end)
{
	pSort::dataType* result = new pSort::dataType[(end - start)];

	long i;

	for (i = start; i < end; i++)
	{
		result[i - start] = rbuf[i];
	}

	return result;
}

void merge(pSort::dataType* result, pSort::dataType* left, pSort::dataType* right, long leftLen, long rightLen)
{
	long i = 0, j = 0;
	while (i < leftLen && j < rightLen)
	{
		if (left[i].key < right[j].key)
		{
			result[i + j] = left[i];
			i++;
		}
		else
		{
			result[i + j] = right[j];
			j++;
		}
	}

	for (; i < leftLen; i++)
	{
		result[i + j] = left[i];
	}
	for (; j < rightLen; j++)
	{
		result[i + j] = right[j];
	}

}

void MERGe(pSort::dataType* rbuf, long len)
{
	if (len <= 1)
	{
		return;
	}
	if (len <= 10)
		INSERTION(rbuf, 0, len);
	pSort::dataType* left = slice(rbuf, 0, len / 2 + 1);
	pSort::dataType* right = slice(rbuf, len / 2, len);


#pragma omp task untied
		{
			MERGe(left, len / 2);
		}
#pragma omp task untied
		{
			MERGe(right, len - (len / 2));
		}
	
#pragma omp taskwait
	merge(rbuf, left, right, len / 2, len - (len / 2));

}


long min(long x, long y) { return (x < y) ? x : y; }


void mergeSort(pSort::dataType *rbuf, long n)
{
	long left_start=0;
#pragma omp parallel num_threads(20) shared(left_start)
	{
		
		for (long curr_size = 1; curr_size <= n - 1; curr_size = 2 * curr_size)
		{
			

#pragma omp for
				for (left_start = 0; left_start < n - 1; left_start += 2 * curr_size)
				{
					
					long mid = min(left_start + curr_size - 1, n - 1);

					long right_end = min(left_start + 2 * curr_size - 1, n - 1);

					merge(rbuf, left_start, mid, right_end);
				}
			
		}
	}
}

/* Function to merge the two haves arr[l..m] and arr[m+1..r] of array arr[] */
void merge(pSort::dataType *rbuf, long l, long m, long r)
{
	long i, j, k;
	long n1 = m - l + 1;
	long n2 = r - m;


	pSort::dataType* L = new pSort::dataType[n1];  L[0].key = 0;
	pSort::dataType* R = new pSort::dataType[n2];  R[0].key = 0;

	for (i = 0; i < n1; i++)
	{
		L[i].key = rbuf[l + i].key;
		for (long ii = 0; ii < LOADSIZE; ii++)
		{
			L[i].payload[ii] = rbuf[l+i].payload[ii];
		}
	}
	for (j = 0; j < n2; j++)
	{
		R[j].key = rbuf[m + 1 + j].key;
		for (long i = 0; i < LOADSIZE; i++)
		{
			R[j].payload[i] = rbuf[m + 1 + j].payload[i];
		}
	}

	i = 0;
	j = 0;
	k = l;
	while (i < n1 && j < n2)
	{
		if (L[i].key <= R[j].key)
		{
			rbuf[k].key = L[i].key;
			for (long ii = 0; ii < LOADSIZE; ii++)
			{
				rbuf[k].payload[ii] = L[i].payload[ii];
			}
			i++;
		}
		else
		{
			rbuf[k].key = R[j].key;
			for (long ii = 0; ii < LOADSIZE; ii++)
			{
				rbuf[k].payload[ii] = R[j].payload[ii];
			}
			j++;
		}
		k++;
	}

	while (i < n1)
	{
		rbuf[k].key = L[i].key;
		for (long ii = 0; ii < LOADSIZE; ii++)
		{
			rbuf[k].payload[ii] = L[i].payload[ii];
		}
		i++;
		k++;
	}

	while (j < n2)
	{
		rbuf[k] = R[j];
		for (long ii = 0; ii < LOADSIZE; ii++)
		{
			rbuf[k].payload[ii] = R[j].payload[ii];
		}
		j++;
		k++;
	}
}
