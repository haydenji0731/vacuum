#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include "tmerge.h"
#include <gclib/GArgs.h>
#include <gclib/GVec.hh>
#include <chrono>
#include <unordered_map>
#include <gclib/GHashMap.hh>

#define VERSION "0.0.1"

const char* USAGE="Vacuum v" VERSION "\n"
                  "==================\n"
                  "The Vacuum utility can take a BAM file and a BED file containing coordinates of spurious junctions.\n"
                  "Junctions of a spliced read are compared against the spurious junctions in the input BED file.\n"
                  "If a BAM record contains >= 1 spurious junctions, then it is removed.\n"
                  "==================\n"
                  "\n"
                  " usage: ./vacuum [-o BAM output file] input.BAM input.BED\n"
                  "\n"
                  " Input arguments (required): \n"
                  "  input.BAM\t\talignment file in SAM/BAM/CRAM format\n"
                  "  input.BED\t\tlist of spurious junctions in BED format\n"
                  "       "
                  "\n"
                  " Optional arguments (-o must be specified):\n"
                  "  -h,--help\tShow this help message and exit\n"
                  "  --version\tShow program version and exit\n"
                  "  -o\tFile for BAM output\n";

GStr inbamname;
GStr inbedname;
GStr outfname;
GStr outfname_removed;
GSamWriter* outfile=NULL;
GSamWriter* removed_outfile=NULL;
std::unordered_map<std::string, int> ht;
GArray<CJunc> spur_juncs;
std::map<std::tuple<std::string, std::string, int, int>, std::vector<PBRec*>> removed_brecs;
int spliced_alignments=0;
bool verbose=false;


struct CJunc {
    int start, end;
    char strand;
    const chr*;
    CJunc(int vs=0, int ve=0, char vstrand='.', const char* vchr='*'):
            start(vs), end(ve), strand(vstrand), chr(vchr){ }

    // overload operators
    bool operator==(const CJunc& a) {
        return (start==a.start && end==a.end && strcmp(chr, a.chr) == 0);
    }

    bool operator<(const CJunc& a) {
        int chr_cmp = strverscmp(chr, a.chr);
        if (chr_cmp == 0) {
            if (start == a.start) {
                return (end < a.end);
            } else {
                return (start < a.start);
            }
        } else {
            return (chr_cmp < 0); //version order 
        }
    }
};


struct PBRec {
    GSamRecord* r;
    PBRec(GSamRecord *rec=NULL):
    r(rec){ }
};


void processOptions(int argc, char **argv);


void loadBed(GStr inbedname) {
    std::ifstream bed_f(inbedname);
    std::string line;
    while (getline(bed_f, line)) {
        GStr gline = line.c_str();
        GVec<GStr> junc;
        int cnt = 0;
        while (cnt < 6) {
            GStr tmp = gline.split("\t");
            junc.Add(gline);
            gline=tmp;
            cnt++;
        }
        char* chr =junc[0].detach();
        //char chr = chrname[strlen(chrname) - 1];
        CJunc j(junc[1].asInt(), junc[2].asInt(), *junc[5].detach(), chr);
        spur_juncs.Add(j);
    }
}


void flushBrec(GVec<PBRec> &pbrecs) {
    if (pbrecs.Count()==0) return;
    for (int i=0; i < pbrecs.Count(); i++) {
        std::string kv = pbrecs[i].r->name();
        std::string tmp = std::to_string(pbrecs[i].r->pairOrder());
        kv += ";";
        kv += tmp;
        if (ht.find(kv) != ht.end()) {
            int new_nh = pbrecs[i].r->tag_int("NH", 0) - ht[kv];
            pbrecs[i].r->add_int_tag("NH", new_nh);
        }
        outfile->write(pbrecs[i].r);
    }
}

int main(int argc, char *argv[]) {
    processOptions(argc, argv);
    loadBed(inbedname);
    GSamReader bamreader(inbamname.chars(), SAM_QNAME|SAM_FLAG|SAM_RNAME|SAM_POS|SAM_CIGAR|SAM_AUX);
    outfile=new GSamWriter(outfname, bamreader.header(), GSamFile_BAM);

    if (outfname_removed.is_empty()) {
        removed_outfile = NULL;
    } else {
        removed_outfile=new GSamWriter(outfname_removed, bamreader.header(), GSamFile_BAM);
    }

    std::cout << "brrrm! identifying alignment records with spurious splice junctions" << std::endl;
    auto start=std::chrono::high_resolution_clock::now();
    int spur_cnt = 0;
    GSamRecord brec;
    std::tuple<std::string, std::string, int, int> prev_key;
    std::tuple<std::string, std::string, int, int> curr_key;
    
    while (bamreader.next(brec)) {
        if (brec.isUnmapped()) {
            continue;
        }

        bam1_t* in_rec = brec.get_b();
        curr_key = std::make_tuple(brec.name(), brec.refName(), brec.get_b()->core.pos, brec.get_b()->core.mpos);

        if (brec.exons.Count() > 1) {
            const char* chrn=brec.refName();
            char strand = brec.spliceStrand();
            bool spur = false;
            for (int i = 1; i < brec.exons.Count(); i++) {
                CJunc j(brec.exons[i-1].end, brec.exons[i].start-1, strand, chr);
                if (spur_juncs.Exists(j)) {
                    spur = true;
                    spur_cnt++;
                    break;
                }
            }
            if (!spur) {
                GSamRecord *rec = new GSamRecord(brec);
                PBRec *newpbr = new PBRec(rec);
                kept_brecs.Add(newpbr);
            } else {
                spur_cnt++;
                std::string kv = brec.name();
                std::string tmp = std::to_string(brec.pairOrder());
                kv += ";";
                kv += tmp;
                // key not present
                if (ht.find(kv) == ht.end()) {
                    ht[kv] = 1;
                } else {
                    int val = ht[kv];
                    val++;
                    ht[kv] = val;
                }
            }
// using GArray
//            } else {
//                CRead cr(brec.pairOrder(), brec.name());
//                int idx;
//                if (spur_reads.Found(cr, idx)) {
//                    spur_reads[idx].incSpurCnt();
//                } else {
//                    spur_reads.Add(cr);
//                }
//            }
        } else {
            GSamRecord *rec = new GSamRecord(brec);
            PBRec* newpbr = new PBRec(rec);
            kept_brecs.Add(newpbr);
        }
    }
    std::cout << "vacuuming completed. writing only clean bam records to the output file." << std::endl;
    flushBrec(kept_brecs);
    bamreader.bclose();
    delete outfile;
    auto end =std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::seconds>(end - start);
    std::cout << spur_cnt << " spurious alignment records were removed." << std::endl;
    std::cout << "Vacuuming completed in " << duration.count() << " seconds" << std::endl;
}

void processOptions(int argc, char* argv[]) {
    GArgs args(argc, argv, "help;verbose;version;SMLPEDVho:r:");
    args.printError(USAGE, true);

    if (args.getOpt('h') || args.getOpt("help")) {
        fprintf(stdout,"%s",USAGE);
        exit(0);
    }

    if (args.getOpt("version")) {
        fprintf(stdout,"%s\n", VERSION);
        exit(0);
    }

    // ifn = input file name
    bool set = false;
    const char* ifn=NULL;

    while ((ifn=args.nextNonOpt()) != NULL) {
        if (!set) {
            inbamname = ifn;
            set=true;
        } else {
            inbedname = ifn;
        }
    }

    if (inbamname == NULL || inbedname == NULL) {
        GMessage(USAGE);
        GMessage("\nError: no input BAM/BED file provided!\n");
        exit(1);
    }

    outfname=args.getOpt('o');
    if (outfname.is_empty()) {
        GMessage(USAGE);
        GMessage("\nError: output filename must be provided.");
        exit(1);
    }

    outfname_removed=args.getOpt('r');
    verbose=(args.getOpt("verbose")!=NULL || args.getOpt('V')!=NULL);
    if (verbose) {
        fprintf(stderr, "Running Vacuum " VERSION ". Command line:\n");
        args.printCmdLine(stderr);
    }

}