
#ifndef SRC_ANALYSIS_COMPOSER_IJ_H_
#define SRC_ANALYSIS_COMPOSER_IJ_H_

namespace espreso {

struct IJ { esint row, column; };

inline bool operator==(const IJ &left, const IJ &right) { return left.row == right.row && left.column == right.column; }
inline bool operator!=(const IJ &left, const IJ &right) { return !(left == right); }
inline bool operator <(const IJ &left, const IJ &right) { return left.row == right.row ? left.column < right.column : left.row < right.row; }

}

#endif /* SRC_ANALYSIS_COMPOSER_IJ_H_ */
