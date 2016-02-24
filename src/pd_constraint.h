#ifndef PD_CONSTRAINT_H
#define PD_CONSTRAINT_H


#include <stdint.h>


/* TODO: DOD, layout for now optimized by solver */
struct PdConstraintAttachment {
        uint32_t i;
        float position[3];
};


struct PdConstraintSpring {
        uint32_t i[2]; /* indices of spring's endpoints */
        float rest_length;
};


#endif
