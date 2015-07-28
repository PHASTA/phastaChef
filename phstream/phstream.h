#ifndef PHSTREAM_H_
#define PHSTREAM_H_
#include <stdio.h>
struct RStream;
struct GRStream;
/** @brief make restart stream */
RStream* makeRStream();
/** @brief clear restart stream */
void clearRStream(RStream* rs);
/** @brief detach output stream */
void destroyRStream(RStream* rs);

/** @brief make geom-restart stream */
GRStream* makeGRStream();
/** @brief clear geom-restart stream */
void clearGRStream(GRStream* grs);
/** @brief destroy geom-restart stream */
void destroyGRStream(GRStream* grs);

/** @brief open restart stream for reading*/
FILE* openRStreamRead(RStream* rs);
/** @brief open restart stream for writing*/
FILE* openRStreamWrite(RStream* rs);

/** @brief open named stream in geom-restart stream for reading*/
FILE* openGRStreamRead(GRStream* grs, const char* named);
/** @brief open named stream in geom-restart stream for writing*/
FILE* openGRStreamWrite(GRStream* grs, const char* named);

/** @brief dev function */
RStream* attachRStream(GRStream* grs);
#endif 
