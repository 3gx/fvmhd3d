#ifndef __myMPI_H__
#define __myMPI_H__

#include <mpi.h>
#include "pboundary.h"
#include "mytimer.h"

namespace myMPI {

	template<class T>
		MPI_Datatype datatype();

	extern unsigned long long data_inout;
	extern double             data_inout_time;

#if 0
	template<class T>
		void  all2all(std::vector<T> sbuf[], std::vector<T> rbuf[], 
				const int myid, 
				const int nproc,
				bool DEBUG = false) 
		{

			MPI_Barrier(MPI_COMM_WORLD);
			double t0 = mytimer::get_wtime();
#if 1
			assert(myid >= 0);
			assert(myid < nproc);
			rbuf[myid] = sbuf[myid];
			for (int dist = 1; dist < nproc; dist++) {
				int src = (nproc + myid - dist) % nproc;
				int dst = (nproc + myid + dist) % nproc;
				int scount = sbuf[dst].size();
				int rcount = 0;
				MPI_Status stat;
				MPI_Sendrecv(&scount, 1, MPI_INT, dst, 0,
						&rcount, 1, MPI_INT, src, 0, MPI_COMM_WORLD, &stat);
				rbuf[src].resize(rcount);
				data_inout += sizeof(T)*rcount;
				MPI_Sendrecv(&sbuf[dst][0], scount, datatype<T>(), dst, 1,
						&rbuf[src][0], rcount, datatype<T>(), src, 1, MPI_COMM_WORLD, &stat);
			}
#else
			static MPI_Status  stat[NMAXPROC  ];
			static MPI_Request  req[NMAXPROC*2];
			static int        nrecv[NMAXPROC  ];

#if 0
			int nreq = 0;
			for (int p = 0; p < nproc; p++) 
			{
				int nsend = sbuf[p].size();
				MPI_Isend(&nsend, 1, MPI_INT, p, myid, MPI_COMM_WORLD, &req[nreq++]);
				data_inout += sizeof(T)*nsend;
			}
			for (int p = 0; p < nproc; p++) 
			{
				MPI_Irecv(&nrecv[p], 1, MPI_INT, p, p, MPI_COMM_WORLD, &req[nreq++]);
			}
			MPI_Waitall(nreq, req, stat);
#else
      static int nsend[NMAXPROC];
      for (int p = 0; p < nproc; p++) 
      {
        nsend[p] = sbuf[p].size();
				data_inout += sizeof(T)*nsend[p]*2;
      }
      MPI_Alltoall(nsend, 1, MPI_INT, nrecv, 1, MPI_INT, MPI_COMM_WORLD);
#endif

			int nreq = 0;
			for (int p = 0; p < nproc; p++) 
      {
				int nsend = sbuf[p].size();
        rbuf[p].resize(nrecv[p]);
				if (nsend    > 0) MPI_Isend(&sbuf[p][0], nsend,    datatype<T>(), p, 0, MPI_COMM_WORLD, &req[nreq++]);
				if (nrecv[p] > 0) MPI_Irecv(&rbuf[p][0], nrecv[p], datatype<T>(), p, 0, MPI_COMM_WORLD, &req[nreq++]);
			}
			MPI_Waitall(nreq, req, stat);

#endif
			MPI_Barrier(MPI_COMM_WORLD);
			data_inout_time += mytimer::get_wtime() - t0;


			if (0)
      {
				int nsend_recv_loc[2] = {0,0};
				for (int p = 0; p < nproc; p++) {
					nsend_recv_loc[0] += sbuf[p].size();
					nsend_recv_loc[1] += rbuf[p].size();
				}

				int nsend_recv_glob[2];
				MPI_Allreduce(&nsend_recv_loc, &nsend_recv_glob, 2, MPI_INT, MPI_SUM, MPI_COMM_WORLD);

				assert(nsend_recv_glob[0] == nsend_recv_glob[1]);    
			}

		}
#endif

  template<bool DO_ALL2ALL, class T>
    void  all2all(std::vector<T> sbuf[], std::vector<T> rbuf[], 
        const int myid, 
        const int nproc,
        const int scale,
        int nsend[NMAXPROC],
        int nrecv[NMAXPROC])
    {
			assert(myid >= 0);
			assert(myid < nproc);

#if 1
      if (DO_ALL2ALL)
      {
        for (int p = 0; p < nproc; p++) 
        {
          assert(sbuf[p].size() % scale == 0);
          nsend[p]    = sbuf[p].size()/scale;
        }
        MPI_Alltoall(nsend, 1, MPI_INT, nrecv, 1, MPI_INT, MPI_COMM_WORLD);
      }

			for (int p = 0; p < nproc; p++) 
			{
				assert(scale*nsend[p] == (int)sbuf[p].size());
				data_inout += sizeof(T) * (nsend[p] + nrecv[p]) * scale;
			}
#else
			for (int p = 0; p < nproc; p++) 
			{
				assert(sbuf[p].size() % scale == 0);
				nsend[p]    = sbuf[p].size()/scale;
				data_inout += sizeof(T) * (nsend[p] + nrecv[p]) * scale;
			}
			MPI_Alltoall(nsend, 1, MPI_INT, nrecv, 1, MPI_INT, MPI_COMM_WORLD);
#endif
      
			double t0 = mytimer::get_wtime();

#if 0 
			// NOT SAFE
			static MPI_Status  stat[NMAXPROC  ];
			static MPI_Request  req[NMAXPROC*2];
			int nreq = 0;
			for (int p = 0; p < nproc; p++) 
				if (p != myid)
				{
					int nsend = sbuf[p].size();
					rbuf[p].resize(scale*nrecv[p]);
					if (nsend    > 0) MPI_Isend(&sbuf[p][0], nsend,          datatype<T>(), p, 1, MPI_COMM_WORLD, &req[nreq++]);
					if (nrecv[p] > 0) MPI_Irecv(&rbuf[p][0], scale*nrecv[p], datatype<T>(), p, 1, MPI_COMM_WORLD, &req[nreq++]);
				}
			rbuf[myid] = sbuf[myid];
			MPI_Waitall(nreq, req, stat);
#endif

#if 0 
			// NOT SAFE
			double t0 = mytimer::get_wtime();
			static MPI_Status stat;
			for (int p = 0; p < nproc; p++) 
				if (p != myid)
				{
					int scount = nsend[p]*scale;
					int rcount = nrecv[p]*scale;
					rbuf[p].resize(rcount);
					if (scount + rcount > 0)
						MPI_Sendrecv(&sbuf[p][0], scount, datatype<T>(), p, 1,
								&rbuf[p][0], rcount, datatype<T>(), p, 1, MPI_COMM_WORLD, &stat);
				}
			rbuf[myid] = sbuf[myid];
#endif

#if 0
			// SAFE but slow, due to hand-shaking between procs even if no data is sent
			static MPI_Status stat;
			for (int dist = 1; dist < nproc; dist++) 
			{
				int src = (nproc + myid - dist) % nproc;
				int dst = (nproc + myid + dist) % nproc;
				int scount = nsend[dst]*scale;
				int rcount = nrecv[src]*scale;
				rbuf[src].resize(rcount);
				MPI_Sendrecv(&sbuf[dst][0], scount, datatype<T>(), dst, 1,
						&rbuf[src][0], rcount, datatype<T>(), src, 1, MPI_COMM_WORLD, &stat);
			}
			rbuf[myid] = sbuf[myid];
#endif
#if 0
			// Should be safe
			static MPI_Status  stat[NMAXPROC  ];
			static MPI_Request  req[NMAXPROC*2];
			for (int dist = 1; dist < nproc; dist++) 
			{
				const int src = (nproc + myid - dist) % nproc;
				const int dst = (nproc + myid + dist) % nproc;
				const int scount = nsend[dst]*scale;
				const int rcount = nrecv[src]*scale;
				rbuf[src].resize(rcount);
				int nreq = 0;
				if (scount > 0) MPI_Isend(&sbuf[dst][0], scount, datatype<T>(), dst, 1, MPI_COMM_WORLD, &req[nreq++]);
				if (rcount > 0) MPI_Irecv(&rbuf[src][0], rcount, datatype<T>(), src, 1, MPI_COMM_WORLD, &req[nreq++]);
				MPI_Waitall(nreq, req, stat);
			}
			rbuf[myid] = sbuf[myid];
#endif
#if 1
			// Must be safe
			static MPI_Status stat;
			for (int dist = 1; dist < nproc; dist++) 
			{
				const int src = (nproc + myid - dist) % nproc;
				const int dst = (nproc + myid + dist) % nproc;
				const int scount = nsend[dst]*scale;
				const int rcount = nrecv[src]*scale;
				rbuf[src].resize(rcount);
				if ((myid/dist) & 1)
				{
					if (scount > 0) MPI_Send(&sbuf[dst][0], scount, datatype<T>(), dst, 1, MPI_COMM_WORLD);
					if (rcount > 0) MPI_Recv(&rbuf[src][0], rcount, datatype<T>(), src, 1, MPI_COMM_WORLD, &stat);
				}
				else
				{
					if (rcount > 0) MPI_Recv(&rbuf[src][0], rcount, datatype<T>(), src, 1, MPI_COMM_WORLD, &stat);
					if (scount > 0) MPI_Send(&sbuf[dst][0], scount, datatype<T>(), dst, 1, MPI_COMM_WORLD);
				}
			}
			rbuf[myid] = sbuf[myid];
#endif
			data_inout_time += mytimer::get_wtime() - t0;
		}


	template<class T>
		void allgather(T &sbuf, std::vector<T> &rbuf, 
				const int myid, 
				const int nproc) {

			rbuf.resize(nproc);
			MPI_Allgather(&sbuf,    1, datatype<T>(), 
					&rbuf[0], 1, datatype<T>(), 
					MPI_COMM_WORLD);

		}

	template<class T>
		void Bcast(std::vector<T> &buf,
				const int myid, 
				const int nproc) {

			int n = buf.size();
			MPI_Bcast(&n, 1, MPI_INT, myid, MPI_COMM_WORLD);
			buf.resize(n);
			MPI_Bcast(&buf[0], n, datatype<T>(), myid, MPI_COMM_WORLD);
		}

	template<class T>
		void Bcast(T &el, const int myid, const int nproc) {
			MPI_Bcast(&el, 1, datatype<T>(), myid, MPI_COMM_WORLD);
		}

	void free_type();


}

#endif // _myMPI_H_
