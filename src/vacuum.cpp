#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <chrono>
#include <unordered_set>
#include <tuple>
#include <map>
#include <unordered_map>
#include <gclib/GArgs.h>
#include <gclib/GStr.h>
#include "GSam.h"

#define VERSION "0.0.1"

const char* USAGE="Vacuum v" VERSION "\n"
                  "==================\n"
                  "The Vacuum utility filters out spurious spliced alignments from SAM/BAM/CRAM files.\n"
                  "Junctions in a spliced alignment are compared against the spurious junctions in the input BED file.\n"
                  "If a BAM record contains >= 1 spurious junctions, then it is removed.\n"
                  "==================\n"
                  "\n"
                  "Usage: ./vacuum [options] input.[SAM|BAM|CRAM] input.BED\n"
                  "\n"
                  "Input arguments (required):\n"
                  "  input.[SAM|BAM|CRAM]   Alignment file in SAM/BAM/CRAM format.\n"
                  "  input.BED              List of spurious junctions in BED format.\n"
                  "\n"
                  "Optional arguments:\n"
                  "  -h,--help              Show this help message and exit.\n"
                  "  --version              Show program version and exit.\n"
                  "  -o                     File for output in BAM/CRAM format. Must be specified.\n"
                  "  -r                     File for removed alignments in BAM format.\n"
                  "  --ref                  Reference genome file for CRAM input.\n"
                  "  --remove_mate          If provided, removes the mate of a spurious read. If not, the mate is unpaired but not removed.\n"
                  "  -V,--verbose           Verbose output.\n";


GStr inbamname;
GStr inbedname;
GStr outfname;
GStr outfname_removed;
GStr cram_ref; //for cram files
GSamWriter* outfile=NULL;
GSamWriter* removed_outfile=NULL;
bool remove_mate=false;
bool verbose=false;
std::unordered_map<std::string, int> ht;
int num_mates=0;
int num_spur_removed=0;
int num_alns_output=0;
int num_spur_alns_both_mates=0;

struct CJunc {
    int start, end;
    char strand;
    const char* chr;
    CJunc(int vs=0, int ve=0, char vstrand='.', const char* vchr="*"):
            start(vs), end(ve), strand(vstrand), chr(vchr){ }

    // overload operators
    bool operator==(const CJunc& a) const {
        return (start==a.start && end==a.end && strcmp(chr, a.chr) == 0);
    }

    bool operator<(const CJunc& a) const {
        int chr_cmp = strcmp(chr, a.chr);
        if (chr_cmp == 0) {
            if (start == a.start) {
                if (end == a.end) {
                    return (strand < a.strand);
                } else {
                    return (end < a.end);
                }
            } else {
                return (start < a.start);
            }
        } else {
            return (chr_cmp < 0);
        }
    }
};

namespace std {
    template<>
    struct hash<CJunc> {
        std::size_t operator()(const CJunc& cj) const {
            std::hash<int> int_hasher;
            std::hash<char> char_hasher;
            std::size_t hash = 17;

            // Hashing start
            hash = hash * 31 + int_hasher(cj.start);

            // Hashing end
            hash = hash * 31 + int_hasher(cj.end);

            // Hashing strand
            hash = hash * 31 + char_hasher(cj.strand);

            // Hashing chr
            const char* chr_ptr = cj.chr;
            std::size_t str_hash = 0;
            while (*chr_ptr) {
                str_hash = str_hash * 31 + char_hasher(*chr_ptr);
                ++chr_ptr;
            }
            hash = hash * 31 + str_hash;

            return hash;
        }
    };
}

struct PBRec {
    GSamRecord* r;
    PBRec(GSamRecord *rec=NULL):
    r(rec){ }
};


void processOptions(int argc, char **argv);


std::unordered_set<CJunc> loadBed(GStr inbedname) {
    std::ifstream bed_f(inbedname);
    std::string line;
    std::unordered_set<CJunc> spur_juncs;
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
        const char* chr =junc[0].detach();
        CJunc j(junc[1].asInt(), junc[2].asInt(), *junc[5].detach(), chr);
        spur_juncs.insert(j);
    }
    return spur_juncs;
}

bool check_identical_cigar(bam1_t* rec1, bam1_t* rec2) {
    if (rec1->core.n_cigar == rec2->core.n_cigar &&
        memcmp(bam_get_cigar(rec1), bam_get_cigar(rec2), rec1->core.n_cigar * sizeof(uint32_t)) == 0) {
            return true;
        }
    return false;
}


void filter_bam(GSamWriter* outfile, GSamWriter* removed_outfile,
                 std::map<std::tuple<std::string, std::string, int, int>, std::vector<PBRec*>>& removed_brecs,
                 GSamReader* bamreader) {

    GSamRecord brec;
    std::map<std::tuple<std::string, std::string, int, int>, int> mates_unpaired; //keep track of count of mates that are unpaired

    //iterate over bam and filter out spurious junctions
    while (bamreader->next(brec)) {
        if (brec.isUnmapped()) {
            continue;
        }

        //check if the alignment needs to be removed
        std::tuple<std::string, std::string, int, int> key = std::make_tuple(brec.name(), brec.refName(),
                                                            brec.get_b()->core.pos, brec.get_b()->core.mpos);
        auto it = removed_brecs.find(key);
        if (it != removed_brecs.end()) {
            bool found = false;
            for (PBRec* item : it->second) {
                bam1_t* in_rec = brec.get_b();
                bam1_t* rm_rec = item->r->get_b();
                if( check_identical_cigar(in_rec, rm_rec) ) {
                    found = true;
                    num_spur_removed++;
                    if (removed_outfile != NULL) {
                        removed_outfile -> write(item->r);
                    }
                    break;  //escape for loop because alignment has been found
                }
            }
            if (found) {
                continue; //resume while loop (otherwise will be written to outfile)
            }
        }

        //write to outfile if alignment is not paired
        if (!brec.isPaired()) {
            outfile->write(&brec);
            num_alns_output++;
            continue;
        }

        //check if the alignment is paired with a removed alignment
        std::tuple<std::string, std::string, int, int> mate_key = std::make_tuple(brec.name(), brec.refName(),
                                                                brec.get_b()->core.mpos, brec.get_b()->core.pos);
        auto it_rem = removed_brecs.find(mate_key);
        if (it_rem != removed_brecs.end()) {
            int num_rem = it_rem->second.size();
            bool update_flag = true;
            //if more then 1 mate needs to be removed, check how many mates have already been unpaired:
            int &num_mts_seen = mates_unpaired[mate_key]; //if not seen, defaults to 0
            if (num_mts_seen == num_rem) {
                update_flag = false; // all mates have been unpaired
            } else {
                num_mts_seen++;
                }

            if (update_flag) {
                num_mates++;
            }

            //write to removed_outfile if remove_mate is true
            if (update_flag && remove_mate) {
                if (removed_outfile != NULL) {
                    removed_outfile->write(&brec);
                }
                continue;
            }

            //update NH tag:
            std::string kv = brec.name();
            std::string tmp = std::to_string(brec.pairOrder());
            kv += ";";
            kv += tmp;
            if (ht.find(kv) != ht.end()) {
            int new_nh = brec.tag_int("NH", 0) - ht[kv];
            brec.add_int_tag("NH", new_nh);
        }

            //update flag, tlen, mpos
            if (update_flag) {
                brec.get_b()->core.flag &= ~3;
                brec.get_b()->core.isize = 0; //set template len to zero
                brec.get_b()->core.mpos =  brec.get_b()->core.pos; //set mate pos to pos
            }
        }

        //write to outfile:
        outfile->write(&brec);
        num_alns_output++;
    }
}


int main(int argc, char *argv[]) {
    std::map<std::tuple<std::string, std::string, int, int>, std::vector<PBRec*>> removed_brecs;
    int spliced_alignments=0;
    processOptions(argc, argv);
    std::unordered_set<CJunc> spur_juncs = loadBed(inbedname);

    GSamReader* bamreader;

    if (inbamname.endsWith(".cram") && !cram_ref.is_empty()) {
        // Initialize for CRAM file with user-specified reference genome
        bamreader = new GSamReader(inbamname.chars(), SAM_QNAME|SAM_FLAG|SAM_RNAME|SAM_POS|SAM_CIGAR|SAM_AUX, cram_ref);
    } else {
        // Initialize for BAM file or CRAM file with header-specified reference genome
        bamreader = new GSamReader(inbamname.chars(), SAM_QNAME|SAM_FLAG|SAM_RNAME|SAM_POS|SAM_CIGAR|SAM_AUX);
    }

    if (outfname.endsWith(".cram")) {
        outfile = new GSamWriter(outfname.chars(), bamreader->header(), GSamFile_CRAM);
    } else {
        outfile = new GSamWriter(outfname.chars(), bamreader->header(), GSamFile_BAM);
    }


    if (outfname_removed.is_empty()) {
        removed_outfile = NULL;
    } else {
        removed_outfile=new GSamWriter(outfname_removed, bamreader->header(), GSamFile_BAM);
    }

    auto start_vacuum=std::chrono::high_resolution_clock::now();
    if (verbose) {
        std::cout << std::endl;
        std::cout << "brrrm! Vacuuming BAM file debris in: " << inbamname << std::endl;
        std::cout << std::endl;
    }

    int num_alignments = 0;
    int num_removed_spliced = 0;
    int num_total_spliced = 0;
    int num_unmapped = 0;
    GSamRecord brec;
    std::tuple<std::string, std::string, int, int> key;

    auto start_flagging=std::chrono::high_resolution_clock::now();
    while (bamreader->next(brec)) {
        num_alignments++;
        if (brec.isUnmapped()) {
            num_unmapped++;
            continue;
        }

        bam1_t* in_rec = brec.get_b();
        key = std::make_tuple(brec.name(), brec.refName(),
                                    brec.get_b()->core.pos, brec.get_b()->core.mpos);

        bool spur = false;
        if (brec.exons.Count() > 1) {
            num_total_spliced++;
            const char* chr=brec.refName();
            char strand = brec.spliceStrand();
            for (int i = 1; i < brec.exons.Count(); i++) {
                CJunc j(brec.exons[i-1].end, brec.exons[i].start-1, strand, chr);
                if (spur_juncs.find(j) != spur_juncs.end()) {
                    spur = true;
                    break;
                }
            }
            if (spur) {
                num_removed_spliced++;
                //add to hash table for NH tag:
                std::string kv = brec.name();
                std::string tmp = std::to_string(brec.pairOrder());
                kv += ";";
                kv += tmp;
                if (ht.find(kv) == ht.end()) { // key not present
                    ht[kv] = 1;
                } else {
                    int val = ht[kv];
                    val++;
                    ht[kv] = val;
                }

                //add spurs to removed_brecs:
                GSamRecord *rec = new GSamRecord(brec);
                PBRec *newpbr = new PBRec(rec);
                if (removed_brecs.find(key) == removed_brecs.end()) {
                    std::vector<PBRec*> v;
                    v.push_back(newpbr);
                    removed_brecs[key] = v;
                } else {
                    removed_brecs[key].push_back(newpbr);
                }

            }
        }
    }

    bamreader->rewind();
    auto end_flagging=std::chrono::high_resolution_clock::now();
    auto duration_flagging = std::chrono::duration_cast<std::chrono::seconds>(end_flagging - start_flagging).count();

    if (verbose) {
        std::cout << "Alignment removal identification completed in: " << duration_flagging << " second(s)" << std::endl;
        std::cout << "Total alignments processed: " << num_alignments << std::endl;
        std::cout << "Unmapped alignments: " << num_unmapped << std::endl;
        std::cout << "Spliced alignments identified: " << num_total_spliced << std::endl;
        std::cout << "Alignments flagged for removal: " << num_removed_spliced << std::endl;
    }


    auto begin_filtering=std::chrono::high_resolution_clock::now();
    filter_bam(outfile, removed_outfile, removed_brecs, bamreader);
    auto end_filtering=std::chrono::high_resolution_clock::now();

    delete outfile;
    if (removed_outfile != NULL) {
        delete removed_outfile;
    }

    auto end_vacuum =std::chrono::high_resolution_clock::now();
    auto duration_vacuum = std::chrono::duration_cast<std::chrono::seconds>(end_vacuum - start_vacuum);

    if (verbose) {
        if (remove_mate) {
            if (removed_outfile != NULL) {
                std::cout << "Mates of spliced alignments removed: " << num_mates << std::endl;
            }
            else {
                std::cout << "Mates of spliced alignments unpaired: " << num_mates << std::endl;}
        }
        std::cout << "Alignments written to output: " << num_alns_output << std::endl;
        std::cout << "Cleaning completed in: " << duration_vacuum.count() << " second(s)" << std::endl;
        std::cout << "Congratulations! Your vacuumed BAM file is now optimized for analysis." << std::endl;
    }
}

void processOptions(int argc, char* argv[]) {
    GArgs args(argc, argv, "help;verbose;version;remove_mate;ref=;SMLPEDVho:r:");
    args.printError(USAGE, true);

    verbose = (args.getOpt("verbose")!= NULL || args.getOpt('V')!= NULL);

    if (args.getOpt("ref") && verbose) {
        cram_ref = args.getOpt("ref");
        std::cout << "Captured --ref: " << cram_ref << std::endl;
    }

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

    if (verbose) {
        std::cout << std::endl;
        fprintf(stderr, "Running Vacuum " VERSION ". Command line:\n");
        args.printCmdLine(stderr);
    }

    remove_mate=(args.getOpt("remove_mate")!=NULL);

}
