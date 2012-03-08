#ifndef __MESHPOINT_H__
#define __MESHPOINT_H__

struct MeshPoint
{
	enum {VIRTUAL_BIT = (1 << 28)};
//	enum {ACTIVE_BIT  = (1 << 27)};
	enum {REFINE_BIT  = (1 << 30), DEREFINE_BIT = (1 << 29)}; 
  enum {REMOVE_BIT  = (1 << 28)};
	enum {NO_BOUNDARY = 0, REFLECTING  = 1, OUTFLOW = 2, DIOD = 3, INFLOW = 4};

  vec3 pos_orig, vel_orig;
  vec3 pos, vel;
  vec3 acc0, acc1;
	real tbeg, tend;
  real dt_new;
	real Volume;
	int rung;
	int boundary;
  unsigned int idx;
  unsigned int status;
  real etaJ;

  void pup(PUP::er &p)
  {
    p|pos_orig;
    p|vel_orig;
    p|pos;
    p|vel;
    p|acc0;
    p|acc1;
    p|tbeg;
    p|tend;
    p|dt_new;
    p|Volume;
    p|rung;
    p|boundary;
    p|idx;
    p|status;
    p|etaJ;
  }

  MeshPoint() : idx(0), status(0), etaJ(0.0) {}

  void set_pos(const vec3 &_pos) {pos = pos_orig = _pos;}
  void set_vel(const vec3 &_vel) {vel = vel_orig = _vel;}

#if 0
	const bool is_active() const {return (status  & ACTIVE_BIT) == ACTIVE_BIT;}
	void      set_active() {status |=  ACTIVE_BIT;}
	void    unset_active() {status &= ~ACTIVE_BIT;}
#endif
	
  const bool is_remove() const {return (status & REMOVE_BIT) == REMOVE_BIT;}
	void      set_remove() {status |=  REMOVE_BIT;}
	void    unset_remove() {status &= ~REMOVE_BIT;}

	const bool is_virtual() const {return false;} // (idx & VIRTUAL_BIT) == VIRTUAL_BIT;}
	void      set_virtual() {} ;//idx |=  VIRTUAL_BIT;}
	void    unset_virtual() {} ;//idx &= ~VIRTUAL_BIT;}

  const bool is_boundary() const {return boundary != NO_BOUNDARY;}

};

#endif // __MESHPOINT_H__

