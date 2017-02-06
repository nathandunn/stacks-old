#include <string>

#include "constants.h"
#include "MetaPopInfo.h"
#include "locus.h"
#include "BamI.h"

using namespace std;

int main(int argc, char** argv) {
    MetaPopInfo mpopi;
    string bam_path;
    Bam* bam_f = new Bam(bam_path.c_str());
    CLocReadSet loc (mpopi);
}
