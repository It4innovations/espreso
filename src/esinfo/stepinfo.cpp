
#include "stepinfo.h"

int espreso::step::loadstep = 0;
int espreso::step::substep = 0;
int espreso::step::iteration = 0;

espreso::step::TYPE espreso::step::type = espreso::step::TYPE::TIME;

int espreso::step::duplicate::instances = 1;
int espreso::step::duplicate::offset = 0;
int espreso::step::duplicate::size = 1;
int espreso::step::duplicate::totalsize = 1;

double espreso::step::time::start = 0;
double espreso::step::time::current = 0;
double espreso::step::time::shift = 1;
double espreso::step::time::final = 1;
double espreso::step::time::precision = 1e-8;

double espreso::step::frequency::start = 0;
double espreso::step::frequency::current = 0;
double espreso::step::frequency::shift = 1;
double espreso::step::frequency::final = 1;
double espreso::step::frequency::precision = 1e-8;
double espreso::step::frequency::angular = 0;

int espreso::step::ftt::step = 0;
int espreso::step::ftt::steps = 1;
double espreso::step::ftt::time = 0;
double espreso::step::ftt::period = 1;
