#include <stdio.h>
#include "hocdec.h"
extern int nrnmpi_myid;
extern int nrn_nobanner_;
modl_reg(){
  if (!nrn_nobanner_) if (nrnmpi_myid < 1) {
    fprintf(stderr, "Additional mechanisms from files\n");

    fprintf(stderr," TCa_d.mod");
    fprintf(stderr," cadecay.mod");
    fprintf(stderr," hpg.mod");
    fprintf(stderr," kA.mod");
    fprintf(stderr," kamt.mod");
    fprintf(stderr," kca3.mod");
    fprintf(stderr," kdrmt.mod");
    fprintf(stderr," kfasttab.mod");
    fprintf(stderr," kslowtab.mod");
    fprintf(stderr," lcafixed.mod");
    fprintf(stderr," nafast.mod");
    fprintf(stderr," naxn.mod");
    fprintf(stderr, "\n");
  }
  _TCa_d_reg();
  _cadecay_reg();
  _hpg_reg();
  _kA_reg();
  _kamt_reg();
  _kca3_reg();
  _kdrmt_reg();
  _kfasttab_reg();
  _kslowtab_reg();
  _lcafixed_reg();
  _nafast_reg();
  _naxn_reg();
}
