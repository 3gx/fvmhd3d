#ifndef __PARTICLE_H__
#define __PARTICLE_H__

struct Particle
{
	enum {VIRTUAL_BIT = (1 << 28)};
	enum {ACTIVE_BIT  = (1 << 27)};
	enum {REFINE_BIT  = (1 << 30), DEREFINE_BIT = (1 << 29)}; 
  enum {REMOVE_BIT  = (1 << 28)};
	enum {NO_BOUNDARY = 0, REFLECTING  = 1, OUTFLOW = 2, DIOD = 3};

  vec3 orig_pos, orig_vel;
  vec3 pos, vel;
  vec3 acc0, acc1;
  real pot;
	real tlast, tend, new_dt;
	real volume, volume_new;
	float rmax;
  long long idx;
  int local_id;
	int rung;
	int boundary;
  unsigned int status;

  real etaOhm;

  Particle() : idx(0), status(0), etaOhm(0.0) {}

	void set_pos(const vec3 &_pos) {pos = orig_pos = _pos;}

	const bool is_active() const {return (status  & ACTIVE_BIT) == ACTIVE_BIT;}
	void      set_active() {status |=  ACTIVE_BIT;}
	void    unset_active() {status &= ~ACTIVE_BIT;}
	
  const bool is_remove() const {return (status & REMOVE_BIT) == REMOVE_BIT;}
	void      set_remove() {status |=  REMOVE_BIT;}
	void    unset_remove() {status &= ~REMOVE_BIT;}

	const bool is_virtual() const {return false;} // (idx & VIRTUAL_BIT) == VIRTUAL_BIT;}
	void      set_virtual() {} ;//idx |=  VIRTUAL_BIT;}
	void    unset_virtual() {} ;//idx &= ~VIRTUAL_BIT;}

  const bool is_boundary() const {return boundary != NO_BOUNDARY;}

};

#endif // __PARTICLE_H__

