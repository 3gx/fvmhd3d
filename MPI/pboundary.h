#ifndef _PBOUNDARY_H_
#define _PBOUNDARY_H_

#include "pfloat.h"

#if 0
#define SMALLB (1.0 + 1.0e-6)
#else
#define SMALLB  1.0
#endif

struct pBoundary {

	pfloat3 pcentre;
	pfloat3  phsize;
	bool empty;

	const bool isempty() const {return empty;}

	pBoundary() : phsize(0.0), empty(true) {}
	pBoundary(const pfloat3 pos) : pcentre(pos), phsize(0.0), empty(false) {}
	pBoundary(const pfloat3 pos, const double h) : pcentre(pos), phsize(h), empty(false) {}

#if 0
	pBoundary(const pfloat3 r1, const pfloat3 r2) {set(r1,   r2);}
	pBoundary& set(pfloat3 r1, pfloat3 r2) {
		r1.x.half();
		r2.x.half();
		pcentre.x = r1.x; pcentre.x.add(r2.x);
		phsize.x  = r2.x;  phsize.x.sub(r1.x);

		r1.y.half();
		r2.y.half(); 
		pcentre.y = r1.y; pcentre.y.add(r2.y);
		phsize.y  = r2.y;  phsize.y.sub(r1.y);

		r1.z.half();
		r2.z.half();  
		pcentre.z = r1.z; pcentre.z.add(r2.z);
		phsize.z  = r2.z;  phsize.z.sub(r1.z);

		empty = false;

		return *this;
	}
#else
	pBoundary(pfloat3 r1, pfloat3 r2) 
	{
#if 0
		r1.x.half();
		r2.x.half();
		pcentre.x = r1.x; pcentre.x.add(r2.x);
		phsize.x  = r2.x;  phsize.x.sub(r1.x);

		r1.y.half();
		r2.y.half(); 
		pcentre.y = r1.y; pcentre.y.add(r2.y);
		phsize.y  = r2.y;  phsize.y.sub(r1.y);

		r1.z.half();
		r2.z.half();  
		pcentre.z = r1.z; pcentre.z.add(r2.z);
		phsize.z  = r2.z;  phsize.z.sub(r1.z);
#else

		pcentre.x = add_half(r1.x, r2.x);
		phsize.x  = r2.x;
		phsize.x.sub(r1.x);
		phsize.x.half();

		pcentre.y = add_half(r1.y, r2.y);
		phsize.y  = r2.y;
		phsize.y.sub(r1.y);
		phsize.y.half();
		
		pcentre.z = add_half(r1.z, r2.z);
		phsize.z  = r2.z;
		phsize.z.sub(r1.z);
		phsize.z.half();
#endif

		empty = false;
	}
	pBoundary(const dvec3 &centre, const dvec3 &hsize) : pcentre(centre), phsize(hsize), empty(false) {}
#endif

	pBoundary& set_x(pfloat<0> x1, pfloat<0> x2) {
		x1.half();
		x2.half();
		pcentre.x = x1; pcentre.x.add(x2);
		phsize.x  = x2;  phsize.x.sub(x1);
		empty = false;    
		return *this;
	}
	pBoundary& set_y(pfloat<1> x1, pfloat<1> x2) {
		x1.half();
		x2.half();
		pcentre.y = x1; pcentre.y.add(x2);
		phsize.y  = x2;  phsize.y.sub(x1);
		empty = false;    
		return *this;
	}
	pBoundary& set_z(pfloat<2> x1, pfloat<2> x2) {
		x1.half();
		x2.half();
		pcentre.z = x1; pcentre.z.add(x2);
		phsize.z  = x2;  phsize.z.sub(x1);
		empty = false;    
		return *this;
	}


	/***************/

	const pfloat3 rmin() const 
	{
		assert(!isempty());
		pfloat3 r = pcentre;
		r.x.sub(phsize.x);
		r.y.sub(phsize.y);
		r.z.sub(phsize.z);
		return r;
	}
	const pfloat3 rmax() const 
	{
		assert(!isempty());
		pfloat3 r = pcentre;
		r.x.add(phsize.x);
		r.y.add(phsize.y);
		r.z.add(phsize.z);
		return r;
	}


	const dvec3 centre() const 
	{
		assert(!isempty());
		return pcentre.to_dvec3();
	}
	const dvec3 hsize() const {
		assert(!isempty());
		return phsize.to_dvec3();
	}

	///////////

	pBoundary& merge(const pBoundary &b) 
	{
		if (b.isempty()) return *this;
		if (isempty()) {
			pcentre = b.pcentre; 
			phsize  = b.phsize;
			empty  = b.empty;
			return *this;
		}

		pfloat_merge(pcentre.x, phsize.x, b.pcentre.x, b.phsize.x, 
				pcentre.x, phsize.x);
		pfloat_merge(pcentre.y, phsize.y, b.pcentre.y, b.phsize.y, 
				pcentre.y, phsize.y);
		pfloat_merge(pcentre.z, phsize.z, b.pcentre.z, b.phsize.z, 
				pcentre.z, phsize.z);

		//     pfloat3 rmin1 =   rmin();
		//     pfloat3 rmax1 =   rmax();
		//     pfloat3 rmin2 = b.rmin();
		//     pfloat3 rmax2 = b.rmax();


		//     rmin1.x.min(rmin2.x);
		//     rmin1.y.min(rmin2.y);
		//     rmin1.z.min(rmin2.z);

		//     rmax1.x.max(rmax2.x);
		//     rmax1.y.max(rmax2.y);
		//     rmax1.z.max(rmax2.z);

		//     return set(rmin1, rmax1);
		return *this;
	}

	/////////

#if 1
	const bool overlap(const pBoundary &b) const 
	{
		if (isempty()) return false;
		//     assert(!isempty());
		pfloat3 dr = pcentre;
		dr.x.psub(b.pcentre.x);
		dr.y.psub(b.pcentre.y);
		dr.z.psub(b.pcentre.z);

		pfloat3 ds = phsize;
		ds.x.add(b.phsize.x);
		ds.y.add(b.phsize.y);
		ds.z.add(b.phsize.z);

		//     fprintf(stderr, "dr= %g %g %g\n", 
		// 	    dr.x.uval,
		// 	    dr.y.uval,
		// 	    dr.z.uval);
		//     fprintf(stderr, "ds= %g %g %g\n", 
		// 	    ds.x.uval,
		// 	    ds.y.uval,
		// 	    ds.z.uval);

		//     bool bx = dr.x.aleq(ds.x);
		//     bool by = dr.y.aleq(ds.y);
		//     bool bz = dr.z.aleq(ds.z);

		//     fprintf(stderr, "bx= %d\n", bx);
		//     fprintf(stderr, "by= %d\n", by);
		//     fprintf(stderr, "bz= %d\n", bz);

#if 0
		fprintf(stderr, " s1  %s \n", dr.x.aleq(ds.x) ? "true" : "false");
		fprintf(stderr, " s2  %s \n", dr.y.aleq(ds.y) ? "true" : "false");
		fprintf(stderr, " s3  %s \n", dr.z.aleq(ds.z) ? "true" : "false");
#endif

		return (
				dr.x.aleq(ds.x) && 
				dr.y.aleq(ds.y) &&
				dr.z.aleq(ds.z));
	} 
#else
	friend const bool overlap(const pBoundary &b1, const pBoundary &b2)
	{
		if (b1.isempty() || b2.isempty()) return false;

		pfloat3 dr = b1.pcentre;
		dr.x.psub(b2.pcentre.x);
		dr.y.psub(b2.pcentre.y);
		dr.z.psub(b2.pcentre.z);

		pfloat3 ds = b1.phsize;
		ds.x.add(b2.phsize.x);
		ds.y.add(b2.phsize.y);
		ds.z.add(b2.phsize.z);

		return (
				dr.x.aleq(ds.x) && 
				dr.y.aleq(ds.y) &&
				dr.z.aleq(ds.z));
	} 
#endif

	const bool isinbox(const pfloat3 pos) const 
	{
		assert(!isempty());
		pfloat3 dr = pcentre;
		dr.x.psub(pos.x);
		dr.y.psub(pos.y);
		dr.z.psub(pos.z);

		pfloat3 ds = phsize;

		return (
				dr.x.aleq(ds.x) && 
				dr.y.aleq(ds.y) &&
				dr.z.aleq(ds.z));
	}

	const bool isinbox(const pBoundary &b) const 
	{
		assert(!isempty());
		pfloat3 dr = pcentre;
		dr.x.psub(b.pcentre.x);
		dr.y.psub(b.pcentre.y);
		dr.z.psub(b.pcentre.z);

		pfloat3 ds = phsize;
		ds.x.sub(b.phsize.x);
		ds.y.sub(b.phsize.y);
		ds.z.sub(b.phsize.z);

#if 1
		const bool flag1 = 
			dr.x.aleq(ds.x) &&  
			dr.y.aleq(ds.y) &&  
			dr.z.aleq(ds.z);
		const bool flag2 = 
			b.phsize.x.aleq(phsize.x) && 
			b.phsize.y.aleq(phsize.y) && 
			b.phsize.z.aleq(phsize.z);
		return flag1 && flag2;
#else
		return (
				dr.x.aleq(ds.x) &&  
				dr.y.aleq(ds.y) &&  
				dr.z.aleq(ds.z));
#endif
	}

	void dump(FILE *fout, bool flag = false) const {
		assert(!isempty());
#if 0 
		fprintf(fout, "c= %20.16lg %20.16lg %20.16lg; hs= %20.16lg %20.16lg %20.16lg", 
				pcentre.x.to_double(), pcentre.y.to_double(), pcentre.z.to_double(),
				phsize.x.to_double(),  phsize.y.to_double(),  phsize.z.to_double());
#else
		const dvec3 min = rmin().to_dvec3();
		const dvec3 max = rmax().to_dvec3();
		fprintf(fout, "min= %20.16lg %20.16lg %20.16lg; max= %20.16lg %20.16lg %20.16lg",
				min.x, min.y, min.z,
				max.x, max.y, max.z);
#endif
		if (flag) fprintf(fout, "\n");
	}
};

#endif // _PBOUNDARY_H_

#ifdef _TOOLBOX_PBOUNDARY_
int main(int argc, char *argv[]) {
	pBoundary b1;
	pBoundary b2;

	b1.pcentre = (pfloat3){0.3f, 0.0f, 0.0f};
	b2.pcentre = (pfloat3){0.8f, 0.0f, 0.0f};

	b1.phsize = (float3){0.3f/2, 0.0f, 0.0f};
	b2.phsize = (float3){0.3f/2, 0.0f, 0.0f};

	pBoundary b3 = b1;
	pBoundary b4 = b2;

	fprintf(stderr, "b1= "); b1.dump(stderr); fprintf(stderr, "\n");
	fprintf(stderr, "b2= "); b2.dump(stderr); fprintf(stderr, "\n");
	fprintf(stderr, "b3= "); 
	//   b1.pmerge(b2);
	b3.merge(b4);

	b1.dump(stderr); 
	fprintf(stderr, "\nb4= "); 
	b3.dump(stderr); 
	fprintf(stderr, "\n");

	return 0;
}
#endif
