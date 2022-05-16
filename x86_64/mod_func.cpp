#include <stdio.h>
#include "hocdec.h"
extern int nrnmpi_myid;
extern int nrn_nobanner_;
#if defined(__cplusplus)
extern "C" {
#endif

extern void _BK_reg(void);
extern void _CaL_reg(void);
extern void _CaN_reg(void);
extern void _CaT_reg(void);
extern void _ccanl_reg(void);
extern void _HCN_reg(void);
extern void _ichan2_reg(void);
extern void _Ka_reg(void);
extern void _Kir_reg(void);
extern void _netstim125_reg(void);
extern void _netstimbox_reg(void);
extern void _SK_reg(void);

void modl_reg() {
  if (!nrn_nobanner_) if (nrnmpi_myid < 1) {
    fprintf(stderr, "Additional mechanisms from files\n");
    fprintf(stderr, " \"BK.mod\"");
    fprintf(stderr, " \"CaL.mod\"");
    fprintf(stderr, " \"CaN.mod\"");
    fprintf(stderr, " \"CaT.mod\"");
    fprintf(stderr, " \"ccanl.mod\"");
    fprintf(stderr, " \"HCN.mod\"");
    fprintf(stderr, " \"ichan2.mod\"");
    fprintf(stderr, " \"Ka.mod\"");
    fprintf(stderr, " \"Kir.mod\"");
    fprintf(stderr, " \"netstim125.mod\"");
    fprintf(stderr, " \"netstimbox.mod\"");
    fprintf(stderr, " \"SK.mod\"");
    fprintf(stderr, "\n");
  }
  _BK_reg();
  _CaL_reg();
  _CaN_reg();
  _CaT_reg();
  _ccanl_reg();
  _HCN_reg();
  _ichan2_reg();
  _Ka_reg();
  _Kir_reg();
  _netstim125_reg();
  _netstimbox_reg();
  _SK_reg();
}

#if defined(__cplusplus)
}
#endif
