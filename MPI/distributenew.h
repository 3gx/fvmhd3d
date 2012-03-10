#ifndef __DISTRIBUTENEW_H__
#define __DISTRIBUTENEW_H__

#include <algorithm>
#include <vector>
#include "peano.h"

template<typename real, typename vec3, typename boundary>
struct DistributeNew
{

	int  nproc;
	boundary domain;
	std::vector<boundary> tiles;
  std::vector<int     > procs;

	DistributeNew() {}

	DistributeNew(const int _nproc, const boundary &_domain) : nproc(_nproc), domain(_domain) 
	{
    assert(
        nproc == 1 ||
        nproc == 2 ||
        nproc == 4 ||
        nproc == 8 ||
        nproc == 16 ||
        nproc == 32 ||
        nproc == 64 ||
        nproc == 128 ||
        nproc == 256 ||
        nproc == 512 ||
        nproc == 1024);
	};

	struct node
	{
		int np;
		boundary bnd;
		node() {}
		node(const int _np, const boundary &_bnd) : np(_np), bnd(_bnd) {}
	};

	std::vector<node> node_list;

	struct sort_x 
	{
		bool operator() (const vec3 &lhs, const vec3 &rhs)
		{
			return lhs.x < rhs.x;
		}
	};
	struct sort_y
	{
		bool operator() (const vec3 &lhs, const vec3 &rhs)
		{
			return lhs.y < rhs.y;
		}
	};
	struct sort_z
	{
		bool operator() (const vec3 &lhs, const vec3 &rhs)
		{
			return lhs.z < rhs.z;
		}
	};

	int lev_max;
	void build_nodes(std::vector<vec3> &pos_list, const boundary &bnd, const int lev)
	{
		assert(pos_list.size() > 2);
		if (lev == lev_max)
		{
			node_list.push_back(node(pos_list.size(), bnd));
			return;
		}

		vec3 xmin(pos_list[0]), xmax(pos_list[0]);
		const int n = pos_list.size();
		for (int i = 1; i < n; i++)
		{
			const vec3 &ipos = pos_list[i];
			xmin = mineach(xmin, ipos); 
			xmax = maxeach(xmax, ipos); 
		}
		const vec3 l = xmax - xmin;

		int split = -1;
		if      (l.x >= l.y && l.x >= l.z) split = 0;
		else if (l.y >= l.x && l.y >= l.z) split = 1;
		else if (l.z >= l.x && l.z >= l.y) split = 2;
		else assert(false);

		assert(split >= 0 && split <= 2);
		
		const int median = n >> 1;
    vec3 rmax = bnd.get_rmax();
		vec3 rmin = bnd.get_rmin();
#if 0	
		if      (split == 0) std::sort(pos_list.begin(), pos_list.end(), sort_x());
		else if (split == 1) std::sort(pos_list.begin(), pos_list.end(), sort_y());
		else if (split == 2) std::sort(pos_list.begin(), pos_list.end(), sort_z());
		else assert(false);
		rmax[split] = rmin[split] = (pos_list[median-1][split] + pos_list[median][split])*0.5;
#else
		if      (split == 0) std::nth_element(pos_list.begin(), pos_list.begin() + median, pos_list.end(), sort_x());
		else if (split == 1) std::nth_element(pos_list.begin(), pos_list.begin() + median, pos_list.end(), sort_y());
		else if (split == 2) std::nth_element(pos_list.begin(), pos_list.begin() + median, pos_list.end(), sort_z());
		else assert(false);
		rmax[split] = rmin[split] = pos_list[median][split];
		if      (split == 0) std::nth_element(pos_list.begin(), pos_list.begin() + median -1, pos_list.end(), sort_x());
		else if (split == 1) std::nth_element(pos_list.begin(), pos_list.begin() + median -1, pos_list.end(), sort_y());
		else if (split == 2) std::nth_element(pos_list.begin(), pos_list.begin() + median -1, pos_list.end(), sort_z());
		else assert(false);
		rmax[split] = rmin[split] = (rmax[split] + pos_list[median-1][split])*0.5;
#endif
		boundary bnd_left(bnd.get_rmin(), rmax);
		boundary bnd_rght(rmin, bnd.get_rmax());

		std::vector<vec3> pos_left, pos_rght;
		pos_left.reserve(  median);
		pos_rght.reserve(n-median);
		std::copy(pos_list.begin(),        pos_list.begin() + median, std::back_inserter(pos_left));
		std::copy(pos_list.begin()+median, pos_list.end(),            std::back_inserter(pos_rght));

		build_nodes(pos_left, bnd_left, lev+1);
		build_nodes(pos_rght, bnd_rght, lev+1);
	}
  
	void determine_division(const std::vector<vec3> &pos_list, const int ntile_in)
	{
		assert(ntile_in >= nproc);
    assert(ntile_in % nproc == 0);
		lev_max = (int)(std::log((real)(ntile_in-1))/std::log((real)2.0)) + 1;

		std::vector<vec3> pos(pos_list);

		node_list.clear();
		build_nodes(pos, domain, 0);

    const int ntile = node_list.size();
    assert(ntile == ntile_in);
    assert(ntile % nproc == 0);

#if 0
    std::random_shuffle(node_list.begin(), node_list.end());
    tiles.resize(ntile);
    procs.resize(ntile);
    for (int i = 0; i < ntile; i++)
    {
      tiles[i] = node_list[i].bnd;
      procs[i] = i % nproc;
    }
#else
    const real size = 2.0*std::max(domain.hsize().x, std::max(domain.hsize().y, domain.hsize().z));
    const int domain_fac = 1.0 / size * (((peano_hilbert::peanokey)1) << (BITS_PER_DIMENSION));

    std::vector<peano_hilbert::peano_struct> keys(ntile);
    for (int tile = 0; tile < ntile; tile++) 
    {
      keys[tile].idx = tile;
      const vec3 pos = node_list[tile].bnd.centre();
      const int x = (int)((pos.x - domain.get_rmin().x) * domain_fac);
      const int y = (int)((pos.y - domain.get_rmin().y) * domain_fac);
      const int z = (int)((pos.z - domain.get_rmin().z) * domain_fac);
      keys[tile].key = peano_hilbert::peano_hilbert_key(x, y, z, BITS_PER_DIMENSION);
    }
    std::sort(keys.begin(), keys.end(), peano_hilbert::cmp_peanokey_index());
    
    tiles.resize(ntile);
    procs.resize(ntile);
    for (int i = 0; i < ntile; i++)
    {
      tiles[i] = node_list[keys[i].idx].bnd;
      procs[i] = i % nproc;
    }
#endif



  }



};

#endif /* DISTRIBUTE */
