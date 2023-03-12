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
#include <set>
#include <tuple>

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
std::unordered_map<std::string, int> ht;
GSamWriter* outfile=NULL;
GSamWriter* removed_outfile=NULL;
bool remove_mate=false;

struct CJunc {
    int start, end;
    char strand;
    const char* chr;
    CJunc(int vs=0, int ve=0, char vstrand='.', const char* vchr="*"):
            start(vs), end(ve), strand(vstrand), chr(vchr){ }

    // overload operators
    bool operator==(const CJunc& a) {
        return (start==a.start && end==a.end && strcmp(chr, a.chr) == 0);
    }

    bool operator<(const CJunc& a) const {
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


std::set<CJunc> loadBed(GStr inbedname) {
    std::ifstream bed_f(inbedname);
    std::string line;
    std::set<CJunc> spur_juncs;
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


// void flushBrec(GVec<PBRec> &pbrecs) {
//     if (pbrecs.Count()==0) return;
//     for (int i=0; i < pbrecs.Count(); i++) {
//         std::string kv = pbrecs[i].r->name();
//         std::string tmp = std::to_string(pbrecs[i].r->pairOrder());
//         kv += ";";
//         kv += tmp;
//         if (ht.find(kv) != ht.end()) {
//             int new_nh = pbrecs[i].r->tag_int("NH", 0) - ht[kv];
//             pbrecs[i].r->add_int_tag("NH", new_nh);
//         }
//         outfile->write(pbrecs[i].r);
//     }
// }

bool check_identical_cigar(bam1_t* rec1, bam1_t* rec2) {
    if (rec1->core.n_cigar == rec2->core.n_cigar && 
        memcmp(bam_get_cigar(rec1), bam_get_cigar(rec2), rec1->core.n_cigar * sizeof(uint32_t)) == 0) {
            return true;
        }
    return false;
}


void filter_bam(GSamReader &bamreader, int unmapped, GSamWriter* outfile, GSamWriter* removed_outfile,
                 std::map<std::tuple<std::string, std::string, int, int>, std::vector<PBRec*>>& removed_brecs) {
    
    GSamRecord prev_brec;
    GSamRecord brec;
    std::tuple<std::string, std::string, int, int> prev_key;
    std::map<std::tuple<std::string, std::string, int, int>, int> mates_unpaired;


    while (bamreader.next(brec)) {
        bool brec_to_outfile = true;

        if (brec.isUnmapped()) {
            continue;
        }
        
        std::tuple<std::string, std::string, int, int> key = std::make_tuple(brec.name(), brec.refName(), 
                                                            brec.get_b()->core.pos, brec.get_b()->core.mpos);
        
        
        auto it = removed_brecs.find(key);
        if (it != removed_brecs.end()) {
            for (PBRec* item : it->second) {
                bam1_t* in_rec = brec.get_b();
                bam1_t* rm_rec = item->r->get_b();
                if( check_identical_cigar(in_rec, rm_rec) ) { 
                    if (removed_outfile != NULL) {
                        removed_outfile -> write(item->r);
                        continue;
                    }   
                }
            } 
        }

        if (!brec.isPaired()) {
            continue;
        }



        std::tuple<std::string, std::string, int, int> mate_key = std::make_tuple(brec.name(), brec.refName(),
                                                                brec.get_b()->core.mpos, brec.get_b()->core.pos);
        
        auto it_rem = removed_brecs.find(mate_key);
        if (it_rem != removed_brecs.end()) {
            int num_rem = it_rem->second.size(); //count of mates that need to be unpaired or removed:
            bool update_flag = true;
            if (num_rem > 1) {
                //check how many mates have already been unpaired:
                auto it_mts = mates_unpaired.find(mate_key);
                int num_mts = it_mts->second; 
                if (num_mts == num_rem) {
                    update_flag = false;
                    break;
                }

                //add mate_key to mates_unpaired:
                if (mates_unpaired.find(mate_key) == mates_unpaired.end()) {
                    mates_unpaired[mate_key] = 1;
                } else {
                    int val = mates_unpaired[mate_key];
                    val++;
                    mates_unpaired[mate_key] = val;
                }
            }

            if (remove_mate) {
                removed_outfile->write(&brec);
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

            //update flag to be unpaired
            if (update_flag) {
                brec.get_b()->core.flag &= 3;
            }
        
        }

        //write to outfile:
        outfile->write(&brec);
    }
}


std::map<std::tuple<std::string, std::string, int, int>, std::vector<PBRec*>> removed_brecs;
int spliced_alignments=0;
bool verbose=false;

int main(int argc, char *argv[]) {
    processOptions(argc, argv);
    std::set<CJunc> spur_juncs = loadBed(inbedname);
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
    int spliced_alignments = 0;
    GSamRecord curr_brec;
    GSamRecord prev_brec;
    std::tuple<std::string, std::string, int, int> prev_key;
    std::tuple<std::string, std::string, int, int> curr_key;
    

    
    while (bamreader.next(curr_brec)) {
        if (curr_brec.isUnmapped()) {
            continue;
        }

        bam1_t* in_rec = curr_brec.get_b();
        curr_key = std::make_tuple(curr_brec.name(), curr_brec.refName(), 
                                    curr_brec.get_b()->core.pos, curr_brec.get_b()->core.mpos);

        bool spur = false;
        if (curr_brec.exons.Count() > 1) {
            spliced_alignments++;
            const char* chr=curr_brec.refName();
            char strand = curr_brec.spliceStrand();
            for (int i = 1; i < curr_brec.exons.Count(); i++) {
                CJunc j(curr_brec.exons[i-1].end, curr_brec.exons[i].start-1, strand, chr);
                if (spur_juncs.find(j) != spur_juncs.end()) {
                    spur = true;
                    break;
                }
            }
            if (spur) {
                spur_cnt++;
                std::string kv = curr_brec.name();
                std::string tmp = std::to_string(curr_brec.pairOrder());
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
            prev_key = curr_key;
        } 
    }

    std::cout << "vacuuming completed. writing only clean bam records to the output file." << std::endl;
    //flushBrec(kept_brecs);
    bamreader.bclose();
    delete outfile;
    if (removed_outfile != NULL) {
        delete removed_outfile;
    }
    auto end =std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::seconds>(end - start);

    if (verbose) {
        std::cout << "Input bam file: " << inbamname.chars() << std::endl;
        std::cout << "Count of spliced alignments: " << spliced_alignments << std::endl;
        std::cout << spur_cnt << " spurious alignment records were removed." << std::endl;
        std::cout << "Vacuuming completed in " << duration.count() << " seconds" << std::endl;
    }
}

void processOptions(int argc, char* argv[]) {
    GArgs args(argc, argv, "help;verbose;version;remove_mate;SMLPEDVho:r:");
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

    remove_mate=(args.getOpt("remove_mate")!=NULL);


}