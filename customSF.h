#ifndef CUSTOMSF_H
#define CUSTOMSF_H

#include "samSz.h"
#include <spr.h>
#include <sam.h>
#include <apfMDS.h>
#include <math.h>
#include <cstdio>

namespace {

  apf::Field* getSprSF(apf::Mesh2* m);

  apf::Field* multipleSF(apf::Mesh* m, apf::Field* sf, double factor);

  apf::Field* cylRefineSF(apf::Mesh2* m, double* box, double ref, double cor);

}

#endif



