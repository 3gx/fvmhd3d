#ifndef _DISTRIBUTE_H_
#define _DISTRIBUTE_H_

#include <algorithm>
#include "pboundary.h"
#include "pOctree.h"

struct Distribute 
{
	struct int3 
	{
		int x, y, z;
		int3() {}
		int3(const int _x, const int _y, const int _z) : x(_x), y(_y), z(_z) {}
		int3(const int _x) : x(_x), y(_x), z(_x) {}
		~int3() {}
	};

  pfloat3 pmin_domain, pmax_domain;

	enum {NLEAF = 8, NIMP = 4};
	pOctree::Tree<NLEAF, NIMP> tree;

	std::vector<pBoundary> tiles;
	
	unsigned int domain_fac;

	int  nproc, ntile;
	int3 nt;

	Distribute() : nproc(-1), ntile(0) {}
	~Distribute() {};

	//////////

	void set(const int _nproc, const int3 _nt, const pBoundary &_global_domain) 
	{
		pmin_domain   = _global_domain.rmin();
		pmax_domain   = _global_domain.rmax();
		nproc         = _nproc;
		nt            = _nt;
		ntile         = nt.x * nt.y * nt.z;

		assert(ntile == nproc);

		tiles.resize(ntile);
	}

	///////////

	template<int ch>
		void sort_coord_array(pfloat3 *r, int lo, int up) const {
			int i, j;
			pfloat3 tempr;
			while ( up>lo ) {
				i = lo;
				j = up;
				tempr = r[lo];

				/*** Split file in two ***/
				while ( i<j ) {
					switch(ch) {
						case 0:
							for ( ; r[j].x.to_double() > tempr.x.to_double(); j-- );
							for ( r[i]=r[j]; i<j && r[i].x.to_double() <= tempr.x.to_double(); i++ );
							break;
						case 1:
							for ( ; r[j].y.to_double() > tempr.y.to_double(); j-- );
							for ( r[i]=r[j]; i<j && r[i].y.to_double() <= tempr.y.to_double(); i++ );
							break;
						case 2:
							for ( ; r[j].z.to_double() > tempr.z.to_double(); j-- );
							for ( r[i]=r[j]; i<j && r[i].z.to_double() <= tempr.z.to_double(); i++ );
							break;
						default: assert(ch >=0 && ch <= 2);
					}
					r[j] = r[i];
				}
				r[i] = tempr;
				/*** Sort recursively, the smallest first ***/
				if ( i-lo < up-i ) { sort_coord_array<ch>(r,lo,i-1);  lo = i+1; }
				else               { sort_coord_array<ch>(r,i+1,up);  up = i-1; }
			}
		}


	template<int ch>
		void calculate_boxdim(const int np, pfloat3 pos[], const int istart, const int iend, 
				pfloat3 &rlow, pfloat3 &rhigh) const {
			switch(ch) {
#if 0
				case 0: rlow.x.set(pfloat<ch>::xmin); rhigh.x.set(pfloat<ch>::xmax); break;
				case 1: rlow.y.set(pfloat<ch>::xmin); rhigh.y.set(pfloat<ch>::xmax); break;
				case 2: rlow.z.set(pfloat<ch>::xmin); rhigh.z.set(pfloat<ch>::xmax); break;
#else
				case 0: rlow.x = pmin_domain.x; rhigh.x = pmax_domain.x; break;
				case 1: rlow.y = pmin_domain.y; rhigh.y = pmax_domain.y; break;
				case 2: rlow.z = pmin_domain.z; rhigh.z = pmax_domain.z; break;
#endif
				default: assert(ch >=0 && ch <= 2);
			}
			if(istart > 0) 
			{
#if 0
				pfloat3 r1 = pos[istart  ];
				pfloat3 r2 = pos[istart-1];
				switch(ch) {
					case 0: r1.x.half(); r2.x.half(); r1.x.add(r2.x); rlow.x = r1.x;break;
					case 1: r1.y.half(); r2.y.half(); r1.y.add(r2.y); rlow.y = r1.y;break;
					case 2: r1.z.half(); r2.z.half(); r1.z.add(r2.z); rlow.z = r1.z;break;
					default: assert(ch >=0 && ch <= 2);
				}
#else
				const pfloat3 r1 =                pos[istart    ];
				const pfloat3 r2 = (istart > 0) ? pos[istart - 1] : pmin_domain;
				switch(ch) {
					case 0: rlow.x = add_half(r1.x, r2.x); break;
					case 1: rlow.y = add_half(r1.y, r2.y); break;
					case 2: rlow.z = add_half(r1.z, r2.z); break;
					default: assert(ch >=0 && ch <= 2);
				}
#endif

			}

			if (iend < np-1) 
			{
#if 0
				pfloat3 r1 = pos[iend  ];
				pfloat3 r2 = pos[iend+1];
				switch(ch) {
					case 0: r1.x.half(); r2.x.half(); r1.x.add(r2.x); rhigh.x = r1.x; break;
					case 1: r1.y.half(); r2.y.half(); r1.y.add(r2.y); rhigh.y = r1.y; break;
					case 2: r1.z.half(); r2.z.half(); r1.z.add(r2.z); rhigh.z = r1.z; break;
					default: assert(ch >= 0 && ch <= 2);
				}
#else
				const pfloat3 r1 =                 pos[iend    ];
				const pfloat3 r2 = (iend < np-1) ? pos[iend + 1] : pmax_domain;
				switch(ch) { 
					case 0: rhigh.x = add_half(r1.x, r2.x); break;
					case 1: rhigh.y = add_half(r1.y, r2.y); break;
					case 2: rhigh.z = add_half(r1.z, r2.z); break;
					default: assert(ch >=0 && ch <= 2);
				}
#endif
			}
		}

	////////////

	void determine_division(std::vector<pfloat3> &pos) {

		std::vector<int> istart(ntile);
		std::vector<int> iend  (ntile);

		const int n = pos.size();
		assert(ntile < n);

		sort_coord_array<0>(&pos[0], 0, n-1);

		for (int i = 0; i < ntile; i++) {
			istart[i] = (i * n)/ntile;
			if (i > 0) iend[i-1] = istart[i] - 1;
		}
		iend[ntile - 1] = n - 1;

		std::vector<pfloat3> rlow(ntile), rhigh(ntile);

		pfloat3 r0, r1;
		for (int ix = 0; ix < nt.x; ix++) {
			const int ix0 =  ix   *nt.y*nt.z;
			const int ix1 = (ix+1)*nt.y*nt.z;
			calculate_boxdim<0>(n, &pos[0], istart[ix0], iend[ix1-1], r0, r1);
			assert(r0.x.uval != r1.x.uval);
			if (r1.x.uval != pmax_domain.x.uval)
				r1.x.dec();
			for (int i = ix0; i < ix1; i++){
				assert(i < ntile);
				rlow [i].x = r0.x;
				rhigh[i].x = r1.x;
			}
		}

		for(int ix = 0; ix < nt.x; ix++) {
			const int ix0 =  ix   *nt.y*nt.z;
			const int ix1 = (ix+1)*nt.y*nt.z;
			const int npy = iend[ix1-1] - istart[ix0] + 1;
			sort_coord_array<1>(&pos[0], istart[ix0], iend[ix1-1]);
			for (int iy = 0; iy < nt.y; iy++) {
				const int iy0 = ix0 +  iy   *nt.z;
				const int iy1 = ix0 + (iy+1)*nt.z;
				calculate_boxdim<1>(npy, &pos[0]+istart[ix0], istart[iy0]-istart[ix0], iend[iy1-1]-istart[ix0], r0,r1);
				assert(r0.y.uval != r1.y.uval);
				if (r1.y.uval != pmax_domain.y.uval)
					r1.y.dec();
				for (int i = iy0; i < iy1; i++){
					assert(i < ntile);
					rlow [i].y = r0.y;
					rhigh[i].y = r1.y;
				}
			}
		}

		for (int ix = 0; ix < nt.x; ix++) {
			const int ix0 = ix*nt.y*nt.z;
			for (int iy = 0; iy < nt.y; iy++) {
				const int iy0 = ix0 +  iy   *nt.z;
				const int iy1 = ix0 + (iy+1)*nt.z;
				const int npz = iend[iy1-1] - istart[iy0] + 1;
				sort_coord_array<2>(&pos[0], istart[iy0], iend[iy1-1]);
				for (int iz = 0; iz < nt.z; iz++) {
					const int iz0 = iy0 + iz;
					calculate_boxdim<2>(npz, &pos[0]+istart[iy0], istart[iz0]-istart[iy0], iend[iz0]-istart[iy0], r0, r1);
					assert(r0.z.uval != r1.z.uval);
					if (r1.z.uval != pmax_domain.z.uval)
						r1.z.dec();
					assert(iz0 < ntile);
					rlow [iz0].z = r0.z;
					rhigh[iz0].z = r1.z;
				}
			}
		}

		tiles.resize(ntile);
		for (int tile = 0; tile < ntile; tile++) 
			tiles[tile] = pBoundary(rlow[tile], rhigh[tile]);

	}

	void build_tree()
	{
		tree.clear();
		const int node_n = (int)(tiles.size() * 100.0/NLEAF);	
		tree.set_domain(pBoundary(pmin_domain, pmax_domain));
		tree.insert(&tiles[0], tiles.size(), 0, node_n);
		tree.root.inner_boundary();
	}


	/////////

	const pBoundary&        get_tile         (const int tile) const {return tiles[tile];}
	const int               get_proc         (const int tile) const {return       tile ;}


	void tiles_overlap(const pBoundary &bnd, std::vector<int> &tiles) const
	{
		tiles.clear();
		tree.root.walk_pBoundary(bnd, tiles);
	}
	void tiles_overlap(const pfloat3 &ppos, std::vector<int> &tiles) const
	{
		tiles.clear();
		tree.root.walk_pfloat3(ppos, tiles);
	}

#if 0
	bool isinproc(const pfloat3 &ppos, const int proc) const {
#if 0
		const int n = procs_tiles[proc].size();
		for (int i = 0; i < n; i++) 
			if (tiles[procs_tiles[proc][i]].isinbox(ppos)) return true;
		return false;
#else
		std::vector<int> tlist;
		tiles_overlap(ppos, tlist);
		assert(tlist.size() > 0);
		if (!(tlist.size() == 1)) {
			fprintf(stderr, " tlist.size()= %d\n", (int)tlist.size());
			for (size_t i = 0; i < tlist.size(); i++) 
			{
				fprintf(stderr, "tile= %d: bnd= ", tlist[i]);
				tiles[tlist[i]].dump(stderr, true);
			}
			const int nt = tlist.size() - 1;
			fprintf(stderr, "tile= %d: bnd= ", tlist[nt]+1);
			tiles[tlist[nt]+1].dump(stderr, true);
		}
		assert(tlist.size() == 1);

		const int nt = tlist.size();
		for (int i = 0; i < nt; i++) 
			if (tiles2proc[tlist[i]] == proc) return true;

		return false;
#endif
	}
#endif

	int which_proc(const pfloat3 &bnd) const {
#if 0
		for (int proc = 0; proc < nproc; proc++) 
			if (isinproc(bnd, proc)) return proc;
		return -1;
#else
		std::vector<int> tlist;
		tiles_overlap(bnd, tlist);
		assert(tlist.size() > 0);
		assert(tlist.size() == 1);
		

		return tlist[0];
#endif
	}

#if 0
	bool overlap(const pBoundary &bnd, const int proc) const {
		const int n = procs_tiles[proc].size();
		for (int i = 0; i < n; i++) 
			if (tiles[procs_tiles[proc][i]].overlap(bnd)) return true;
		return false;
	}
	bool which_overlap(const pBoundary &bnd, std::vector<int> &procs) {
		procs.clear();
		procs.reserve(32);
		for (int proc = 0; proc < nproc; proc++) {
			if (overlap(bnd, proc)) procs.push_back(proc);
		}
		return ((int)procs.size() > 0);
	};
#endif

	void which_tiles_overlap(
			const pBoundary &bnd, 
			std::vector<int> &out_tiles)
	{
		tiles_overlap(bnd, out_tiles);
	};

#if 0
	std::pair<int, int> which_tile(const pfloat3 &pos) const {
#if 0
		for (int proc = 0; proc < nproc; proc++) {
			const int n = procs_tiles[proc].size();
			for (int i = 0; i < n; i++) 
				if (tiles[procs_tiles[proc][i]].isinbox(pos))
					return std::pair<int, int>(proc, procs_tiles[proc][i]);
		}
		return std::pair<int, int>(-1, -1);
#else
		std::vector<int> tlist;
		tiles_overlap(pos, tlist);
		assert(tlist.size() > 0);
		assert(tlist.size() == 1);

#if 0
		const int nt = tlist.size();

		int tile = tlist[0];
		int proc = tiles2proc[tile];
		for (int i = 1; i < nt; i++)
		{
			const int itile = tlist[i];
			const int iproc = tiles2proc[itile];
			if (iproc < proc && itile < tile)
			{
				tile = itile;
				proc = iproc;
			}
		}
		return std::pair<int, int>(proc, tile);
#else
		return std::pair<int, int>(tiles2proc[tlist[0]], tlist[0]);
#endif
#endif
	}
#endif


};

#endif // _DISTRIBUTE_H_

