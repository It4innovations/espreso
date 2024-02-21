
#ifndef SRC_MATH_PRIMITIVES_SLICE_H_
#define SRC_MATH_PRIMITIVES_SLICE_H_

struct Slice {
    esint start, end, step;

    Slice(): start(0), end(0), step(1) {}
    Slice(esint start, esint end = 0, esint step = 1): start(start), end(end), step(step) {}

    Slice evaluate(esint size)
    {
        Slice s;
        if (end == 0) {
            end = size;
        } else if (end < 0) {
            end = size + end;
        }
        if (start < 0) {
            start = size + start;
        }
        return s;
    }
};

#endif /* SRC_MATH_PRIMITIVES_SLICE_H_ */
