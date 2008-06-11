#include <config.h>

#include <devdefs.h>

#include "compitf.h"
#include "compext.h"
#include "compinit.h"


SPICEdev COMPinfo = {
    {
        "Comparator",
        "Simple comparator",

        &COMPnSize,
        &COMPnSize,
        COMPnames,

        &COMPpTSize,
        COMPpTable,

        &COMPmPTSize,
        COMPmPTable,

#ifdef XSPICE
/*----  Fixed by SDB 5.2.2003 to enable XSPICE/tclspice integration  -----*/
        NULL,  /* This is a SPICE device, it has no MIF info data */

        0,     /* This is a SPICE device, it has no MIF info data */
        NULL,  /* This is a SPICE device, it has no MIF info data */

        0,     /* This is a SPICE device, it has no MIF info data */
        NULL,  /* This is a SPICE device, it has no MIF info data */

        0,     /* This is a SPICE device, it has no MIF info data */
        NULL,  /* This is a SPICE device, it has no MIF info data */
/*---------------------------  End of SDB fix   -------------------------*/
#endif

	0
    },

    DEVparam      : COMPparam,
    DEVmodParam   : COMPmParam,
    DEVload       : COMPload,
    DEVsetup      : COMPsetup,
    DEVunsetup    : NULL,
    DEVpzSetup    : COMPsetup,
    DEVtemperature: COMPtemp,
    DEVtrunc      : NULL,
    DEVfindBranch : NULL,
    DEVacLoad     : COMPacload,  /* ac load and normal load are identical */
    DEVaccept     : NULL,
    DEVdestroy    : COMPdestroy,
    DEVmodDelete  : COMPmDelete,
    DEVdelete     : COMPdelete,
    DEVsetic      : NULL,
    DEVask        : COMPask,
    DEVmodAsk     : COMPmodAsk,
    DEVpzLoad     : COMPpzLoad,
    DEVconvTest   : NULL,     /* COMPconvTest, XXXX experimental */
    DEVsenSetup   : COMPsSetup,
    DEVsenLoad    : COMPsLoad,
    DEVsenUpdate  : NULL,
    DEVsenAcLoad  : COMPsAcLoad,
    DEVsenPrint   : COMPsPrint,
    DEVsenTrunc   : NULL,
    DEVdisto      : NULL,
    DEVnoise      : COMPnoise,
#ifdef CIDER
    DEVdump       : NULL,
    DEVacct       : NULL,
#endif                        
    DEVinstSize   : &COMPiSize,
    DEVmodSize    : &COMPmSize

};


SPICEdev *
get_comp_info(void)
{
    return &COMPinfo;
}
