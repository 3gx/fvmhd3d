#include "fvmhd3d.h"

namespace fvmhd3d
{

  globalDomains::globalDomains() {};
  globalDomains::globalDomains(CkMigrateMessage *msg) {}

  void globalDomains::setDomains(const CkVec<boundary> &bnd, CkCallback &cb)
  {
    assert(numElements == (int)bnd.size());
    proc_domains.resize(numElements);
    for (int i = 0; i < numElements; i++)
      proc_domains[i] = bnd[i];

    const int node_n = (int)(numElements * 100.0/NLEAF);

    proc_tree.clear();
    proc_tree.set_domain(
        boundary(
          global_domain.centre() - global_domain.hsize()*1.5,
          global_domain.centre() + global_domain.hsize()*1.5));
    proc_tree.insert(&proc_domains[0], numElements, 0, node_n);
    proc_tree.root.inner_boundary();

    contribute(cb);
  }
  
  void globalDomains::pup(PUP::er &p)
  {
    CBase_globalDomains::pup(p);

    int n = proc_domains.size();
    p|n;

    if (p.isUnpacking())
    {
      proc_domains.resize(n);
      assert(n == numElements);
    }


    PUParray(p, &proc_domains[0], numElements);


    if (p.isUnpacking())
    {
      const int node_n = (int)(numElements * 100.0/NLEAF);

      proc_tree.clear();
      proc_tree.set_domain(
          boundary(
            global_domain.centre() - global_domain.hsize()*1.5,
            global_domain.centre() + global_domain.hsize()*1.5));
      proc_tree.insert(&proc_domains[0], numElements, 0, node_n);
      proc_tree.root.inner_boundary();
    }
  }


}
