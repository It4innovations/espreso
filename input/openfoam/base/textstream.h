#ifndef TEXTSTREAM_H
#define TEXTSTREAM_H

#include <fstream>


class TextStream : public std::ofstream {

public:
    TextStream(const char *name) : std::ofstream(name), indent(0) {}

    void incIndent() {
        indent += 1;
    }

    void decIndent() {
        indent -= 1;
    }

    int getIndent() {
        return indent;
    }

    void beginLine() {
        for (int i = 0; i < indent; i++) {
            (*this) << "\t";
        }
    }

protected:
    int indent;

};

template<typename T> void write(TextStream &ts, const T * value) {
    ts << value;
}

template<typename T> void write(TextStream &ts, const T &value) {
    value.write(ts);
}

template<> inline void write(TextStream &ts, const bool &value) {
    if (value) {
        ts << "true";
    } else {
        ts << "false";
    }
}

template<> inline void write(TextStream &ts, const char *value) {
    ts << value;
}

template<> inline void write(TextStream &ts, const int &value) {
    ts << value;
}

template<> inline void write(TextStream &ts, const unsigned int &value) {
    ts << value;
}

template<> inline void write(TextStream &ts, const double &value) {
    ts << value;
}

template<> inline void write(TextStream &ts, const std::string &value) {
    ts << value;
}

template<> inline void write(TextStream &ts, const mesh::Point &value) {
    ts << "(" << value.x << " " << value.y << " " << value.z << ")";
}

template<typename T> void write(TextStream &ts, const std::vector<T*> &value) {
    ts << value.size() << "(\n";
    ts.incIndent();
    for (int i = 0; i < value.size(); i++) {
        ts.beginLine();
        write(ts, *value[i]);
        ts << "\n";
    }
    ts.decIndent();
    ts.beginLine();
    ts << ")";
}

template<typename T> void write(TextStream &ts, const std::vector<T> &value) {
    ts << value.size() << "(\n";
    ts.incIndent();
    for (int i = 0; i < value.size(); i++) {
        ts.beginLine();
        write(ts, value[i]);
        ts << "\n";
    }
    ts.decIndent();
    ts.beginLine();
    ts << ")";
}

#endif // TEXTSTREAM_H
