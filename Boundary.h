#ifndef __BOUNDARY_H__
#define __BOUNDARY_H__


#ifndef HUGE
#define HUGE 1e30
#endif

template<class fvec3> 
struct Boundary 
{
  protected:
    fvec3 rmin, rmax;
    fvec3  _hsize, _centre;

  public:

    void pup(PUP::er &p)
    {
      p|rmin;
      p|rmax;
      _centre = (rmax + rmin)*0.5;
      _hsize  = (rmax - rmin)*0.5;
    }

    Boundary() : rmin(fvec3(HUGE, HUGE, HUGE)), rmax(fvec3(-HUGE, -HUGE, -HUGE)) {}
    Boundary(const fvec3 &pos) : rmin(pos), rmax(pos) {_centre = pos; _hsize = 0.0;}
    Boundary(const fvec3 &_rmin, const fvec3 &_rmax) : rmin(_rmin), rmax(_rmax) {
      _centre = (rmax + rmin)*0.5;
      _hsize  = (rmax - rmin)*0.5;
    }

    const fvec3 centre() const { return _centre;} 
    const fvec3 hsize()  const { return _hsize; }
    const fvec3& get_rmin() const { return  rmin;}
    const fvec3& get_rmax() const { return  rmax;}


    Boundary& merge(const Boundary &b) 
    {
      rmin.x = (b.rmin.x < rmin.x) ? b.rmin.x : rmin.x;
      rmin.y = (b.rmin.y < rmin.y) ? b.rmin.y : rmin.y;
      rmin.z = (b.rmin.z < rmin.z) ? b.rmin.z : rmin.z;
      rmax.x = (b.rmax.x > rmax.x) ? b.rmax.x : rmax.x;
      rmax.y = (b.rmax.y > rmax.y) ? b.rmax.y : rmax.y;
      rmax.z = (b.rmax.z > rmax.z) ? b.rmax.z : rmax.z;
      _centre = (rmax + rmin)*0.5;
      _hsize  = (rmax - rmin)*0.5;
      return *this;
    }

    const bool overlap(const Boundary &b, const fvec3 &Lbox) const 
    {
      const fvec3 s =   hsize() + b.hsize();
      const	fvec3 l = (centre() - b.centre()).abseach();
      const fvec3 lp(
          (l.x < 0.5*Lbox.x) ? l.x: std::abs(Lbox.x - l.x),
          (l.y < 0.5*Lbox.y) ? l.y: std::abs(Lbox.y - l.y),
          (l.z < 0.5*Lbox.z) ? l.z: std::abs(Lbox.z - l.z));
      return lp.x <= s.x && lp.y <= s.y && lp.z <= s.z;
    }

    const bool overlap(const Boundary &b, const fvec3 &Lbox, bool &flag) const 
    {
      const fvec3 s =   hsize() + b.hsize();
      const	fvec3 l = (centre() - b.centre()).abseach();
      const bool fx = l.x > 0.5*Lbox.x;
      const bool fy = l.y > 0.5*Lbox.y;
      const bool fz = l.z > 0.5*Lbox.z;

      flag |= fx | fy | fz;
      const fvec3 lp(
          l.x*(1-fx) + fx * std::abs(Lbox.x - l.x),
          l.y*(1-fy) + fy * std::abs(Lbox.y - l.y),
          l.z*(1-fz) + fz * std::abs(Lbox.z - l.z));
      return lp.x <= s.x && lp.y <= s.y && lp.z <= s.z;
    }


    // periodic isinbox
    const bool isinbox(const Boundary &b, const fvec3 &Lbox) const 
    {
      const fvec3 s =   hsize() - b.hsize();
      const	fvec3 l = (centre() - b.centre()).abseach();
      const fvec3 lp(
          (l.x < 0.5*Lbox.x) ? l.x: Lbox.x - l.x,
          (l.y < 0.5*Lbox.y) ? l.y: Lbox.y - l.y,
          (l.z < 0.5*Lbox.z) ? l.z: Lbox.z - l.z);
      return lp.x <= s.x && lp.y <= s.y && lp.z <= s.z;
    }

    const bool not_overlap(const Boundary &b) const {
      return 
        (rmax.x < b.rmin.x) || (b.rmax.x < rmin.x)
        || (rmax.y < b.rmin.y) || (b.rmax.y < rmin.y)
        || (rmax.z < b.rmin.z) || (b.rmax.z < rmin.z);

    }
    const bool overlap(const Boundary &b) const {return !not_overlap(b);}

    const bool isinbox(const fvec3 pos) const {
      return 
        rmin.x <= pos.x && pos.x <= rmax.x &&
        rmin.y <= pos.y && pos.y <= rmax.y &&
        rmin.z <= pos.z && pos.z <= rmax.z;
    }

    const bool isinbox(const Boundary &b) const {
      return
        rmin.x <= b.rmin.x && b.rmax.x <= rmax.x &&
        rmin.y <= b.rmin.y && b.rmax.y <= rmax.y &&
        rmin.z <= b.rmin.z && b.rmax.z <= rmax.z;
    }

    void dump(FILE *fout, bool flag = false) const {
      fprintf(fout, "min= %20.16lg %20.16lg %20.16lg; max= %20.16lg %20.16lg %20.16lg", 
          rmin.x, rmin.y, rmin.z,
          rmax.x,  rmax.y,  rmax.z);
      if (flag) fprintf(fout, "\n");
    }
};

#endif // __BOUNDARY_H__
