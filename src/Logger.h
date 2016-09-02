#include <iostream>
#include <fstream>
//#include <iomanip>

/* TeeBuf
 * ==========
 * From http://wordaligned.org/articles/cpp-streambufs */
class TeeBuf: public std::streambuf {
public:
    // Construct a streambuf which tees output to both input
    // streambufs.
    TeeBuf(std::streambuf* sb1, std::streambuf* sb2)
        : sb1(sb1) , sb2(sb2)
    	{}

private:
    std::streambuf* sb1;
    std::streambuf* sb2;
    // This tee buffer has no buffer. So every character "overflows"
    // and can be put directly into the teed buffers.
    virtual int overflow(int c) {
        if (c == EOF) {
            return !EOF;
        } else {
            int r1 = sb1->sputc(c);
            int r2 = sb2->sputc(c);
            return r1 == EOF || r2 == EOF ? EOF : c;
        }
    }

    // Sync both teed buffers.
    virtual int sync() {
        int r1 = sb1->pubsync();
        int r2 = sb2->pubsync();
        return r1 == 0 && r2 == 0 ? 0 : -1;
    }
};

/* TeeOStream
 * ==========
 * From http://wordaligned.org/articles/cpp-streambufs */
class TeeOStream : public std::ostream {
public:
    // Construct an ostream which tees output to the supplied
    // ostreams.
    TeeOStream(std::ostream& o1, std::ostream& o2)
        : std::ostream(&buf), buf(o1.rdbuf(), o2.rdbuf())
        {}

private:
    TeeBuf buf;
};

/* Logger
 * ========== */
struct Logger {
    std::ofstream l;
    std::ostream& o;
    std::ostream& e;

    TeeOStream lo;
    TeeOStream le;

    Logger(std::ofstream&& logf)
        : l(std::move(logf))
        , o(std::cout)
        , e(std::cerr)
        , lo(l, o)
        , le(l, e)
        {}
};
