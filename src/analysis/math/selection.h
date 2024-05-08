
#ifndef SRC_ANALYSIS_MATH_SELECTION_H_
#define SRC_ANALYSIS_MATH_SELECTION_H_

struct Selection {
    esint offset, size, step;

    Selection(): offset(0), size(1), step(1) {}
    Selection(esint offset, esint size, esint step): offset(offset), size(size), step(step) {}

    bool operator==(const Selection &other) const { return offset == other.offset && size == other.size && step == other.step; }
};

#endif /* SRC_ANALYSIS_MATH_SELECTION_H_ */
