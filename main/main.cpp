//
//  main.cpp
//  bloom
//
//  Created by Karel Brinda on 29.04.15.
//  and Maciek Sykulski
//

#include <vector>
#include <getopt.h>
#include <err.h>
#include <omp.h>

#include <queue>
#include <utility>

#include "main.hpp"

#include <boost/regex.hpp>

#define MIN_BLOOM_CAPACITY 1024.

int main(int argc,char** argv){
    if (argc<2){
        return usage();
    }
    
    if (strcmp(argv[1], "stats") == 0) {
        return main_stats(argc-1, argv+1);
    }
    else if (strcmp(argv[1], "create") == 0) {
        return main_create(argc-1, argv+1);
    }
    else if (strcmp(argv[1], "create_many") == 0) {
        return main_create_many(argc-1, argv+1);
    }
    else if (strcmp(argv[1], "bitwise") == 0) {
        return main_bitwise(argc-1, argv+1);
    }
    else if (strcmp(argv[1], "dump") == 0) {
        return main_dump(argc-1, argv+1);
    }
    else if (strcmp(argv[1], "shrink") == 0) {
        return main_shrink(argc-1, argv+1);
    }
    else if (strcmp(argv[1], "hamming") == 0) {
        return main_hamming(argc-1, argv+1);
    }
    else if (strcmp(argv[1], "symmdiffmat") == 0) {
        return main_symmdiffmat(argc-1, argv+1);
    }
    else if (strcmp(argv[1], "query") == 0) {
        return main_query(argc-1, argv+1);
    }
    else if (strcmp(argv[1], "query_and_split") == 0) {
        return main_query_and_split(argc-1, argv+1);
    }
    else if (strcmp(argv[1], "merge") == 0) {
        return main_merge(argc-1, argv+1);
    }
    else {
        fprintf(stderr, "Unrecognized command '%s'\n", argv[1]);
        return 1;
    }
    
    return 0;
}

int usage(){
    fprintf(stderr, "\n");
    fprintf(stderr, "Program: ulysses\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "Usage:   ulysses <command> [options]\n\n");
    fprintf(stderr, "Command: create        create a Bloom filter from a FASTA file\n");
    fprintf(stderr, "         create_many   create many Bloom filters from large FASTA file(s)\n");
    fprintf(stderr, "         bitwise       bitwise operations on Bloom filters\n");
    fprintf(stderr, "         shrink        get a shrinked Bloom filter\n");
    fprintf(stderr, "         stats         print statistics about a Bloom filter\n");
    fprintf(stderr, "         dump          dump of a Bloom filter\n");
    fprintf(stderr, "         hamming       print Hamming distance matrix for given Bloom filters\n");
    fprintf(stderr, "         symmdiffmat   print set symmetric distance matrix for given Bloom filters\n");
    fprintf(stderr, "         query         query a Bloom filter\n");
    fprintf(stderr, "         query_and_split    query a Bloom filter and split output k-mer wise to classified/nonclassified\n");
    fprintf(stderr, "         merge         merge two fasta files k-mer wise, with readnames of the format >read_name|rnum|5|ind|011001101001|\n");
    fprintf(stderr, "\n");
    return 1;
}

void estimate_bloom_size(bitvector::size_type size, bitvector::size_type ones, uint32_t k,
                         double & dens, double & est){
    dens=((double)ones)/size;
    est=-((double)size)*log1p(-dens)/k;
}


void compute_print_bloom_stats(const Bloom & bf){
    bitvector::size_type ones=bf.ones();
    double dens;
    double est;
    
    estimate_bloom_size(bf.array.size(),ones,bf.nh,dens,est);    
    
    printf("num_hash_functions\t%u\n",bf.nh);
    printf("Array_size_b\t%lu\n",bf.array.size());
    printf("Array_size_B\t%lu\n",bf.array.size()/8);
    printf("Array_size_MB\t%.2f\n",1.0*(bf.array.size()/8)/(1024*1024));
    printf("Ones\t%lu\n",ones);
    printf("Density\t%f\n",dens);
    printf("Est_num_items\t%ld\n",(long)round(est));
    if (bf.seed.span>0){
        printf("Seed\t%s\n",bf.seed.seedstr);
        printf("Seed_weight\t%d\n",bf.seed.weight);
        printf("Seed_span\t%d\n",bf.seed.span);
    }   
}

int main_stats(int argc,char** argv){

    if (argc != 2) {
        fprintf(stderr, "\n");
        fprintf(stderr, "Usage:   ulysses stats <in.bf>\n");
        fprintf(stderr, "\n");
        return 1;
    }

    Bloom bf(0,0,NULL);
    bf.load(argv[1]);
    
    compute_print_bloom_stats(bf);
    
    return 0;
}

bool seedstr_correct(const char * seedstr){
    for (const char * c = seedstr; *c; ++c)
        if ((*c != '#' )&&(*c!='-'))
            return false;
    return true;
}


int main_create(int argc,char** argv){
    
    int c;

    const char *seedstr=DEFAULT_SEED;
    coor as_B=1048576;
    int nh=6;
    int r=0;

    while ((c = getopt(argc, argv, "a:s:h:r")) >= 0) {
        switch (c) {
            case 'a':
                as_B=atol(optarg);
                break;
            case 's':
                seedstr = strdup(optarg);
                break;
            case 'h':
                nh=atol(optarg);
                break;
            case 'r':
                r=1;
                break;
            default:
                return 1;
        }
    }
    
    if ((optind + 2 != argc)|| (!seedstr_correct(seedstr))) {
        fprintf(stderr, "\n");
        fprintf(stderr, "Usage:   ulysses create [options] <in.fa> <out.bf>\n");
        fprintf(stderr, "\n");
        fprintf(stderr, "Options: -s STR seed [" DEFAULT_SEED "]\n");
        fprintf(stderr, "         -a INT size of the array in bytes [1048576]\n");
        fprintf(stderr, "         -h INT number of hash functions [6]\n");
        fprintf(stderr, "         -r     include also reverse complements of kmers\n");
        fprintf(stderr, "\n");
        return 1;
    }
    
    char *fa_fn=argv[optind];
    char *bf_fn=argv[optind+1];

    Bloom bf(as_B*8,nh,seedstr);
    bf.create(fa_fn,r);
    bf.save(bf_fn);
    return 0;
}

// Returns rounded up size of optimal bloom filter with k-hash functions and n elements
bitvector::size_type get_size_in_bytes(bitvector::size_type k,bitvector::size_type n){
    //k=m/n*ln2    
    //ln(p) = - m/n *(ln2)^2 = -k*ln(2)
    //p = 2^-k    
    bitvector::size_type m = ceil((double)(n*k)/log(2.));
    return (bitvector::size_type)pow(2.,ceil(log2(m)-3.));
}


int main_create_many(int argc,char** argv){
    int c;

    const char *seedstr=DEFAULT_SEED;
    coor as_B=1048576;
    int nh=6;
    int r=0;
    
    std::string Multi_fasta_filename = "/dev/fd/0"; //default use stdin
    std::string ID_to_taxon_map_filename;
    std::string exclude_bloom_filename;
    std::string include_bloom_filename;
    std::string initial_bloom_filename;
    
    bool print_stats = false;
    bool shrink = false;
    
    while ((c = getopt(argc, argv, "a:s:h:rF:m:e:n:tk")) >= 0) {
        switch (c) {
            case 'a':
                as_B=atol(optarg);
                break;
            case 's':
                seedstr = strdup(optarg);
                break;
            case 'h':
                nh=atol(optarg);
                break;
            case 'r':
                r=1;
                break;
            case 'F' :
                Multi_fasta_filename = optarg;
                break;
            case 'm' :
                ID_to_taxon_map_filename = optarg;
                break;	
            case 'i':
                initial_bloom_filename = optarg;
                break;            
            case 'e':
                exclude_bloom_filename = optarg;
                break;
            case 'n':
                include_bloom_filename = optarg;
                break;                
            case 't':
                print_stats = true;
                break;
            case 'k':
                shrink = true;
                break;
            default:
                return 1;
        }
    }
    
    if ((optind + 1 != argc)|| (!seedstr_correct(seedstr)) || 
	Multi_fasta_filename.empty() || ID_to_taxon_map_filename.empty()) {
        fprintf(stderr, "\n");
        fprintf(stderr, "Usage:   ulysses create_many [options] -m <file.map> <out.dir>\n");
        fprintf(stderr, "\n");
        fprintf(stderr, "Options: -s STR seed [" DEFAULT_SEED "]\n");
        fprintf(stderr, "         -a INT size of the array in bytes [1048576]\n");
        fprintf(stderr, "         -h INT number of hash functions [6]\n");
        fprintf(stderr, "         -F <input.fa> optional input fasta file [ default stdin ]\n");
        fprintf(stderr, "         -m <file.map> required map file mapping read names to taxonomic ids.\n");
        fprintf(stderr, "         -i BLM initial bloom filter to OR with (must have hash and size the same)\n");
        fprintf(stderr, "         -r     include also reverse complements of kmers\n");
        fprintf(stderr, "         -t       print statistics of the final BF\n");
        fprintf(stderr, "         -e BLM exclude k-mers present in BLM bloom filter\n");
        fprintf(stderr, "         -n BLM include only k-mers present in BLM bloom filter\n");
        fprintf(stderr, "         -k       shrink to optimal size before saving\n");
        fprintf(stderr, "\n");        
        fprintf(stderr, "This program reads fasta format from standard input as default\n");
        fprintf(stderr, "(or from a file given after -F parameter).\n");
        fprintf(stderr, "The mapping <file.map> should contain lines of the format:\n");
        fprintf(stderr, "<read_name> <taxonomic_id>\n");
        fprintf(stderr, "i.e. two strings separated with a tab character.\n");
        fprintf(stderr, "In the output directory <out.dir> bloom filter files will be created\n");
        fprintf(stderr, "with names of the format:\n");
        fprintf(stderr, "<taxonomic_id>_a<size_of_the_array_in_bytes>_h<number_of_hash_functions>_s<seed_string>_r<0|1>.bf\n");
        fprintf(stderr, "Only those reads which map their names in <file.map> are used,\n");
        fprintf(stderr, "and as many bloom filters are created as there are taxonomic_ids\n");    
        fprintf(stderr, "which have some reads map to them in the input file.\n");
        fprintf(stderr, "\n");
        return 1;
    }
    
    //char *fa_fn=argv[optind];
    char *bf_fn=argv[optind]; //Output directory
    
    fprintf(stderr, "Reading seqid_to_taxon map... ");
    read_ID_to_taxon_map(ID_to_taxon_map_filename);
    fprintf(stderr, "DONE.\n");
    
    Bloom initial_bf(0,0,NULL);
    if (initial_bloom_filename.size()>0){
        initial_bf.load(initial_bloom_filename.c_str());
        assert((initial_bf.nh==(bitvector::size_type)nh)&&(initial_bf.array.size()==(bitvector::size_type)(as_B*8)));
    }
    
    Bloom exclude_bf(0,0,NULL);
    if (exclude_bloom_filename.size()>0){
        exclude_bf.load(exclude_bloom_filename.c_str());
    }

    Bloom include_bf(0,0,NULL);
    if (include_bloom_filename.size()>0){
        include_bf.load(include_bloom_filename.c_str());
    }

    std::map<std::string,Bloom> * taxon_bloom_map = bloom_create_many_blooms(
        initial_bloom_filename.size()>0?&initial_bf:NULL,       
        exclude_bloom_filename.size()>0?&exclude_bf:NULL,
        include_bloom_filename.size()>0?&include_bf:NULL,
                             as_B*8,nh,seedstr,Multi_fasta_filename.c_str(),r);
    
    //save all of them to directory
    auto map_fun = [=]( std::pair<const std::string,Bloom>& p) {  
                    if (shrink){
                        double dens, est;        
                        estimate_bloom_size(p.second.array.size(),p.second.ones(),p.second.nh,dens,est);
                        est = est>MIN_BLOOM_CAPACITY?est:MIN_BLOOM_CAPACITY;
                        bitvector::size_type new_size = get_size_in_bytes(p.second.nh,round(est));
                        double factor = (p.second.array.size()/8)/new_size;
                        if ((long)(factor)>1){
                            fprintf(stderr,"Shrinking to size:%ld\n",(p.second.array.size()/8)/(long)(factor));
                            p.second.shrink((long)factor);
                        }
                    }
                    std::string fname = std::string(bf_fn)+"/"
                                      +p.first
                                      +std::string("_a")+std::to_string(p.second.array.size()/8)
                                      +std::string("_h")+std::to_string(p.second.nh)
                                      +std::string("_s")+std::string(seedstr)
                                      +std::string("_r")+std::to_string(r)                                      
                                      +std::string(".bf");
                    fprintf(stderr, "Saving bloom filter to file %s\n",fname.c_str());
                    if (print_stats){
                        printf("filename\t%s\n",fname.c_str());
                        compute_print_bloom_stats(p.second);
                    }
                    p.second.save(fname.c_str()); 
                   };
                   
    std::for_each(taxon_bloom_map->begin(), taxon_bloom_map->end(), map_fun);
    //bloom_save(&bf,bf_fn);  
        
    delete taxon_bloom_map;
    
    return 0;
}


int main_bitwise(int argc,char** argv){
    int c;
    
    char buf_ops[1000];
    char *buf_bfs[1000];
    
    int ops=0;
    bool print_stats = false;
    bool shrink = false;    
    
    while ((c = getopt(argc, argv, "a:o:x:tk")) >= 0) {
        switch (c) {
            case 'a':
                buf_ops[ops]='a';
                buf_bfs[ops]=strdup(optarg);
                ops++;
                break;
            case 'o':
                buf_ops[ops]='o';
                buf_bfs[ops]=strdup(optarg);
                ops++;
                break;
            case 'x':
                buf_ops[ops]='x';
                buf_bfs[ops]=strdup(optarg);
                ops++;
                break;
            case 't':
                print_stats = true;
                break;
            case 'k':
                shrink = true;
                break;
            default:
                return 1;
        }
    }
    
    
    if (optind + 2 != argc) {
        fprintf(stderr, "\n");
        fprintf(stderr, "Usage:   ulysses bitwise [options] <in.bf> <out.bf>\n");
        fprintf(stderr, "\n");
        fprintf(stderr, "Options: -a STR   bitwise AND of other BF\n");
        fprintf(stderr, "         -o STR   bitwise OR of other BF\n");
        fprintf(stderr, "         -x STR   bitwise XOR of other BF\n");
        fprintf(stderr, "         -t       print statistics of the final BF\n");
        fprintf(stderr, "         -k       shrink before saving to optimal size\n");
        fprintf(stderr, "Use '-' as <out.bf> not to save the final BF.\n");
        fprintf(stderr, "\n");
        return 1;
    }

    Bloom bf(0,0,NULL);
    bf.load(argv[optind]);
    
    for(int i=0;i<ops;i++){
        Bloom bf2(0,0,NULL);
        bf2.load(buf_bfs[i]);
        switch (buf_ops[i]) {
            case 'a':
                bf &= bf2;
                break;
            case 'o':
                bf |= bf2;
                break;
            case 'x':
                bf ^= bf2;
                break;
            default:
                return 1;
        }            
    }
    
    if (shrink){
        double dens, est;        
        estimate_bloom_size(bf.array.size(),bf.ones(),bf.nh,dens,est);
        est = est>MIN_BLOOM_CAPACITY?est:MIN_BLOOM_CAPACITY;
        bitvector::size_type new_size = get_size_in_bytes(bf.nh,round(est));        
        //fprintf(stderr,"Shrinking to size:%ld\n",new_size);
        bf.shrink((bf.array.size()/8)/new_size);
    }
    
    if (print_stats)
        compute_print_bloom_stats(bf);
           
    if (strcmp(argv[optind+1], "-")!=0)    
        bf.save(argv[optind+1]);
        
    return 0;
}


int main_dump(int argc,char** argv){

    int c;
    
    int mode=0;

    while ((c = getopt(argc, argv, "xb")) >= 0) {
        switch (c) {
            case 'x':
                mode=1;
                break;
            case 'b':
                mode=2;
                break;
            default:
                return 1;
        }
    }
    
    if (optind + 1 != argc) {
        fprintf(stderr, "\n");
        fprintf(stderr, "Usage:   ulysses dump <in.bf>\n");
        fprintf(stderr, "\n");
        fprintf(stderr, "Options: -x     hexadecimal mode\n");
        fprintf(stderr, "         -b     bit mode\n");
        fprintf(stderr, "\n");
        return 1;
    }

    Bloom bf(0,0,NULL);
    bf.load(argv[optind]);

    if(mode==0){
        for(bitvector::size_type i=0;i<bf.array.size();i++){
            if(bf.array.test(i)){
		printf("%lu\n",i);
            }            
        }
    }
    else if (mode==1){
        for(bitvector::size_type i=0;i<bf.array.size()/8;++i){
	    if (i % 20 == 0){
		    printf("\n");
		}
	    int hex = 0;
	    for(int j=0;j<8;++j){
		hex|=bf.array[8*i+j]<<j;
	    }		
            printf("%02X ",hex);
        }
    }
    else if(mode==2){
        for(bitvector::size_type i=0;i<bf.array.size();i++){
            if (i%80==0){
                printf("\n");
            }
            if(bf.array.test(i)){
                    printf("1");
            }
            else{
                 printf("0");
            }            
        }
    }
    
    printf("\n");
    
    return 0;
}

int main_shrink(int argc,char** argv){
    
    int c;
    
    long factor=2;
    
    while ((c = getopt(argc, argv, "f:")) >= 0) {
        switch (c) {
            case 'f':
                factor=atol(optarg);
                break;
            default:
                return 1;
        }
    }
    
    if (optind + 2 != argc) {
        fprintf(stderr, "\n");
        fprintf(stderr, "Usage:   ulysses shrink [options] <in.bf> <out.bf>\n");
        fprintf(stderr, "\n");
        fprintf(stderr, "Options: -f INT   factor [2]\n");
        fprintf(stderr, "\n");
        return 1;
    }
    
    Bloom bf(0,0,NULL);
    bf.load(argv[optind]);
    bf.shrink(factor);
    bf.save(argv[optind+1]);
        
    return 0;
}


int main_symmdiffmat(int argc,char** argv){
    
    int c;
    
    while ((c = getopt(argc, argv, "")) >= 0) {
        switch (c) {
            default:
                return 1;
        }
    }
    
    if (optind + 2 > argc) {
        fprintf(stderr, "\n");
        fprintf(stderr, "Usage:   ulysses symmdiffmat <in1.bf> <in2.bf> ...\n");
        fprintf(stderr, "\n");
        return 1;
    }
    
    int n=argc-optind;
    
    std::vector<Bloom> bfs_vec;
    std::vector<double> s_est;    
    std::vector< std::vector<double> > sdiff( n, std::vector<double>(n) ); 
    
    bfs_vec.resize(n);
    s_est.resize(n);
    double dens;
    for(int i=0;i<n;i++){
        bfs_vec[i].load(argv[optind+i]);                    
        estimate_bloom_size(bfs_vec[i].array.size(),bfs_vec[i].ones(),bfs_vec[i].nh,dens,s_est[i]); 
    }
    
    Bloom bf_tmp(0,bfs_vec[0].nh,NULL);

    for(int i=0;i<n;i++){
        for(int j=0;j<=i;j++){
            if (i==j){
                sdiff[i][i]=0.;
                continue;
            }            
            bf_tmp=bfs_vec[i] | bfs_vec[j];
            double dens;
            double est;
            estimate_bloom_size(bf_tmp.array.size(),bf_tmp.ones(),bf_tmp.nh,dens,est); 
            sdiff[i][j]=sdiff[j][i]= (est - s_est[i]) + (est - s_est[j]);
        }
    }
            
    
    for(int i=0;i<n;i++){
        for(int j=0;j<n;j++){
            printf("%f ",sdiff[i][j]);
        }
        printf("\n");
    }
    
    return 0;
}


int main_hamming(int argc,char** argv){
    
    int c;
    
    while ((c = getopt(argc, argv, "")) >= 0) {
        switch (c) {
            default:
                return 1;
        }
    }
    
    if (optind + 2 > argc) {
        fprintf(stderr, "\n");
        fprintf(stderr, "Usage:   ulysses hamming <in1.bf> <in2.bf> ...\n");
        fprintf(stderr, "\n");
        return 1;
    }
    
    int n=argc-optind;
    std::vector<Bloom> bfs_vec;
    std::vector< std::vector<int> > hamming( n, std::vector<int>(n) ); 
    
    bfs_vec.resize(n);
    for(int i=0;i<n;i++){
        //bf[i].init(0,0,NULL);
        bfs_vec[i].load(argv[optind+i]);
    }
    
    Bloom bf_tmp(0,bfs_vec[0].nh,NULL);

    for(int i=0;i<n;i++){
        for(int j=0;j<=i;j++){
            if (i==j){
                hamming[i][i]=0;
                continue;
            }
            bf_tmp=bfs_vec[i]^bfs_vec[j];
            hamming[i][j]=hamming[j][i]=bf_tmp.ones();
        }
    }
            
    
    for(int i=0;i<n;i++){
        for(int j=0;j<n;j++){
            printf("%10d ",hamming[i][j]);
        }
        printf("\n");
    }
    
    return 0;
}

KSEQ_INIT(gzFile, gzread)

int main_query(int argc,char** argv){
    
    int c;
    
    //int k=0;
    int d=0;
    
    while ((c = getopt(argc, argv, "r")) >= 0) {
        switch (c) {
//             case 'k':
//                 k=1;
//                 break;
            case 'r':
                d=1;
                break;
            default:
                return 1;
        }
    }
    
    if (optind + 2 != argc) {
        fprintf(stderr, "\n");
        fprintf(stderr, "Usage:   ulysses query [options] <in.bf> <in.fq>\n");
        fprintf(stderr, "\n");
//        fprintf(stderr, "Options: -k   reports for all kmers\n");
        fprintf(stderr, "Options: -r       test also reverse directions of reads\n");
        fprintf(stderr, "\n");
        return 1;
    }
    
    Bloom bf(0,0,NULL);
    bf.load(argv[optind]);
    
    gzFile fp;
    kseq_t *seq;
    int l;
    
    fp = gzopen(argv[optind+1], "r");
    assert(fp != Z_NULL);
    seq = kseq_init(fp);
        
    unsigned int bytes_compr_kmer=compressed_kmer_size(bf.seed.weight);
    
    while ((l = kseq_read(seq)) >= 0) { // STEP 4: read sequence
        
        for (int dir=0;dir<=d;dir++){
            printf("%s\t%s\t",seq->name.s, dir ? "r" : "f");
            for (int i=0;i<l-bf.seed.span+1;i++){
                int res=bf.query((uchar*)&(seq->seq.s[i]), dir, bytes_compr_kmer);
                if (res==1){
                    printf("1");
                }
                else if (res==0){
                    printf("0");
                }
                if (res==-1){
                    printf("x");
                }
            }
            printf("\n");
        }
    }
    kseq_destroy(seq);
    gzclose(fp);
    
    return 0;
}

#define RNUM_IND_REGEX "\\|rnum\\|([0-9]+)\\|ind\\|([01]+)\\|"


void classify_read(std::string &name, std::string &seqs, Bloom & bf, 
                   int d, unsigned int bytes_compr_kmer,
                   boost::regex & expression,
                   uint64_t & rnum,
                   std::ostringstream &class_oss,
                   std::ostringstream &unclass_oss){

    std::string::iterator start, end;    
    boost::match_results<std::string::iterator> match;

    const char * seqs_char = seqs.c_str();
    
    start = name.begin();
    end = name.end();            
    boost::regex_search(start,end, match, expression, boost::match_default);
    std::string::iterator kmer_it;
    if (match[0].matched){
        //std::string s_rnum(match[1].first,match[1].second);            
        //rnum = std::stol(s_rnum);
        kmer_it = match[2].first;
    } else {
        if ((long)seqs.length()-bf.seed.span+1<=0)  {
            //Skip read if its too short, but preserve numbering
            ++rnum;
            return;
        }
        std::string ones(seqs.length()-bf.seed.span+1,'1');
        name += "|rnum|"+std::to_string(rnum)+"|ind|"+ones+"|";
        start = name.begin()+ (end-start);
        end = name.end();
        boost::regex_search(start,end, match, expression, boost::match_default);
        kmer_it = match[2].first;
    }
    ++rnum;
    std::string name_notfound(name);
    std::string::iterator kmer_it_nf = name_notfound.begin() + (kmer_it - name.begin());
    bool kmer_found = false;
    bool kmer_notfound = false;
    for (long i=0;i<((long)seqs.length()-bf.seed.span+1);i++){
        if (*kmer_it=='1'){
            bool hit = false;
            bool nokmer = false;
            for (int dir=0;dir<=d;dir++){
                int res=bf.query((const uchar*)&(seqs_char[i]), dir, bytes_compr_kmer);
                if (res==-1) {
                    nokmer = true;
                    //break;
                }
                hit|=(res>0);
            }
            if (nokmer) {
                *kmer_it = '0';
                *kmer_it_nf = '0';                        
            } else if (hit){
                *kmer_it = '1'; kmer_found = true;
                *kmer_it_nf = '0';                                            
            } else { //no hit
                *kmer_it = '0';
                *kmer_it_nf = '1'; kmer_notfound = true;                        
            }                    
        }
        ++kmer_it;
        ++kmer_it_nf;
    }
    if (kmer_found){
        class_oss << ">" << name << "\n" 
                  << seqs << "\n";
    }
    if (kmer_notfound){
        unclass_oss << ">" << name_notfound << "\n" 
                  << seqs << "\n";
    }                    
}


const size_t DEF_WORK_UNIT_SIZE = 500000;

class OutputWorkUnit {
public:
    uint64_t priority;
    std::string classified;
    std::string unclassified;
    
    OutputWorkUnit(uint64_t p, std::string &&cs, std::string &&ucs):
            priority(p),classified(std::move(cs)),unclassified(std::move(ucs))
            {};
    ~OutputWorkUnit(){};
    
    OutputWorkUnit(const OutputWorkUnit& other):
      priority( other.priority ),
      classified( other.classified ),
      unclassified( other.unclassified ){} 
    
    OutputWorkUnit(OutputWorkUnit && other):
      priority( std::move(other.priority) ),
      classified( std::move(other.classified)),
      unclassified( std::move(other.unclassified)){}     
    
    OutputWorkUnit& operator=(OutputWorkUnit&& other){
        if (this != &other) {
            std::swap(priority,other.priority);
            std::swap(classified,other.classified);
            std::swap(unclassified,other.unclassified);
        }
        return *this;
    }
     
    OutputWorkUnit& operator=(const OutputWorkUnit&); // no implementation 
    
    bool operator>(const OutputWorkUnit & ou) const{
        return priority>ou.priority;
    };
};


int main_query_and_split(int argc,char** argv){
    
    int Num_threads = 1;
    #ifdef _OPENMP
    omp_set_num_threads(Num_threads);   
    #endif
    size_t Work_unit_size = DEF_WORK_UNIT_SIZE;        
    
    int c;    
    int d=0;
    long long sig;
    
    while ((c = getopt(argc, argv, "rt:u:")) >= 0) {
        switch (c) {
            case 'r':
                d=1;
                break;
            case 't' :
                sig = atoll(optarg);
                if (sig <= 0)
                    errx(EX_USAGE, "can't use nonpositive thread count");
                #ifdef _OPENMP
                if (sig > omp_get_num_procs())
                    errx(EX_USAGE, "thread count exceeds number of processors");
                Num_threads = sig;
                omp_set_num_threads(Num_threads);
                #endif
                break;
            case 'u' :
                sig = atoll(optarg);
                if (sig <= 0)
                    errx(EX_USAGE, "can't use nonpositive work unit size");
                Work_unit_size = sig;
                break;            
            default:
                return 1;
        }
    }
    
    if (optind + 4 != argc) {
        fprintf(stderr, "\n");
        fprintf(stderr, "Usage:   ulysses query_and_split [options] <in.bf> <in.fa> <out_found.fa> <out_notfound.fa>\n");
        fprintf(stderr, "\n");
        fprintf(stderr, "Options: -r       test also reverse directions of reads\n");
        fprintf(stderr, "         -t #     Number of threads\n");
       fprintf(stderr,  "         -u #     Thread work unit size (in bp, def.=%lu)\n", DEF_WORK_UNIT_SIZE); 
        fprintf(stderr, "\n");
        fprintf(stderr, "This command produces two fasta files with additional strings in read names\n");
        fprintf(stderr, "of the format: >read_name|rnum|5|ind|011001101001|\n");
        fprintf(stderr, "where a number after 'rnum' is a read number, either rewritten from the input file\n");
        fprintf(stderr, "or assigned here counting from 0. Bit flags after 'ind' indicate which kmers \n");
        fprintf(stderr, "are present in the bloom filter (the file <out_found.fa>) \n");
        fprintf(stderr, "and which kmers are not present in the bloom filter  (the file <out_notfound.fa>)\n");
        fprintf(stderr, "Each of these output files can be later used as an input file to this command\n");
        fprintf(stderr, "to further split sets of kmers in reads with respect of presence in other \n");
        fprintf(stderr, "bloom filters.\n");
        fprintf(stderr, "\n");
        return 1;
    }
    
    Bloom bf(0,0,NULL);
    bf.load(argv[optind]);
    
    gzFile fp;
    kseq_t *seq;
                
    fp = gzopen(argv[optind+1], "r");
    assert(fp != Z_NULL);
    seq = kseq_init(fp);
    
    FILE *ffp, *nffp;    
    ffp = fopen(argv[optind+2], "w");
    nffp = fopen(argv[optind+3], "w");
    
    if (ffp == NULL || nffp == NULL) {
        fprintf(stderr, "Can't open output files for writing!\n");
        exit(1);
    }
        
    unsigned int bytes_compr_kmer=compressed_kmer_size(bf.seed.weight);
        
    long l = 0;
    uint64_t total_workunits = 0;
    uint64_t workunit_num_output_now = 0;
    uint64_t total_sequences = 0;
    uint64_t total_sequences_processed = 0;
    uint64_t total_bases = 0;    
    
    std::priority_queue< OutputWorkUnit, 
                         std::vector< OutputWorkUnit >, 
                         std::greater< OutputWorkUnit > > out_pq;
    uint64_t recent_max_pq = 0;
                         
    boost::regex expression(RNUM_IND_REGEX);
    
    #pragma omp parallel
    {
        uint64_t rnum;
        uint64_t workunit_num;
        std::ostringstream classified_output_ss, unclassified_output_ss;
        std::vector<std::pair<std::string,std::string> > work_unit;
        
        while (l>=0) { //while reader is valid, non-critical check
            work_unit.clear();
            size_t total_nt = 0;
            bool empty_workunit = true;
            #pragma omp critical(get_input)
            {
                rnum = total_sequences; // rewrite reads numbering
                workunit_num = total_workunits;
                while (total_nt < Work_unit_size) {
                    l = kseq_read(seq);
                    if (! (l>=0))
                        break;
                    empty_workunit = false;
                    work_unit.emplace_back(seq->name.s, seq->seq.s);                    
                    total_nt += work_unit.back().second.length();
                    ++total_sequences;
                }
                if (!empty_workunit)
                    ++total_workunits;
            } // end critical
            
            if (empty_workunit)
                break;
            
            classified_output_ss.str("");
            unclassified_output_ss.str("");            
            for (size_t j = 0; j < work_unit.size(); j++)
                    classify_read(work_unit[j].first, work_unit[j].second, 
                                  bf, d, bytes_compr_kmer,
                                  expression, rnum,
                                  classified_output_ss, unclassified_output_ss);
                    
            #pragma omp critical(write_output)
            {
                total_sequences_processed += work_unit.size();
                total_bases += total_nt;
                std::cerr << "\rProcessed " << total_sequences_processed << " sequences (" << total_bases << " bp) ...";

#ifndef NDEBUG                
                if (!out_pq.empty()) {                    
                    if ((recent_max_pq+10)<out_pq.size()){
                        recent_max_pq=out_pq.size();
                        std::cerr << "Output queue of size " << out_pq.size() <<"\n";
                    }
                }
#endif
                if (workunit_num==workunit_num_output_now) {
                    fprintf(ffp, classified_output_ss.str().c_str());
                    fprintf(nffp, unclassified_output_ss.str().c_str());
                    ++workunit_num_output_now;
                                        
                    while (!out_pq.empty() && 
                            out_pq.top().priority==workunit_num_output_now) {                                        
                        fprintf(ffp, out_pq.top().classified.c_str());
                        fprintf(nffp,out_pq.top().unclassified.c_str()); 
                        out_pq.pop();
                        ++workunit_num_output_now;
                    }
                } else {
                    out_pq.emplace(workunit_num, classified_output_ss.str(),
                                               unclassified_output_ss.str());                    
                }                
            } //end critical
        }
    } //end parallel section
    
    assert(out_pq.empty());
    //fprintf(stderr,"Last kseq_read:%d",kseq_read(seq));
    kseq_destroy(seq);
    gzclose(fp);
    fclose(ffp);
    fclose(nffp);
    
    return 0;
}




int main_merge(int argc,char** argv){
    
    int c;
    
    //int k=0;
    bool intersection=false;
    
    std::string output_file = "/dev/fd/1"; //default use stdout
    
    while ((c = getopt(argc, argv, "iF:")) >= 0) {
        switch (c) {
            case 'i':
                intersection=true;
                break;
            case 'F' :
                output_file = optarg;
                break;
            default:
                return 1;
        }
    }
    
    if ((optind + 2 != argc)||output_file.empty()) {
        fprintf(stderr, "\n");
        fprintf(stderr, "Usage:   ulysses merge [options] <in1.fa> <in2.fa>\n");
        fprintf(stderr, "\n");
        fprintf(stderr, "Options: -i    output only intersection of two inputs\n");
        fprintf(stderr, "Options: -F <out.fa>   output to file (default stdout)\n");
        fprintf(stderr, "\n");     
        fprintf(stderr, "This option merges two fasta files with read names\n");
        fprintf(stderr, "of the format: >read_name|rnum|5|ind|011001101001|\n");
        fprintf(stderr, "The input files are assumed to be sorted by 'rnum'.\n");     
        fprintf(stderr, "The output is a fasta file sorted by rnum, where 'ind' bit flags\n");     
        fprintf(stderr, "in reads with the same 'rnum' are merged by OR-ing bits.\n");     
        fprintf(stderr, "When -i flag is used only common reads are output, with their bitflags\n");     
        fprintf(stderr, "merged by AND-ing bits.\n");     
        fprintf(stderr, "\n");     
        return 1;
    }
    
    gzFile fp1, fp2;
    kseq_t *seq1, *seq2;
    int l1, l2;
            
    fp1 = gzopen(argv[optind], "r");
    assert(fp1 != Z_NULL);
    seq1 = kseq_init(fp1);
    
    fp2 = gzopen(argv[optind+1], "r");
    assert(fp2 != Z_NULL);
    seq2 = kseq_init(fp2);
    
    FILE *outf;
    outf = fopen(output_file.c_str(), "w");    
    
    if (outf == NULL){ 
        fprintf(stderr, "Can't open output file for writing!\n");
        exit(1);
    }
            
    //unsigned int rnum = 0;
    
        
    std::string::iterator start, end;
    boost::regex expression(RNUM_IND_REGEX);
    boost::match_results<std::string::iterator> match1, match2;
    boost::match_flag_type flags = boost::match_default;        
    
    std::string name1,name2;
        
    l1 = kseq_read(seq1);
    l2 = kseq_read(seq2);
    long rnum1=-1,  rnum2=-1;
    bool isread1=false, isread2=false;
    do {                 
        if (l1>=0 && !isread1) {
            name1 = std::string(seq1->name.s);                        
            start = name1.begin();
            end = name1.end();            
            boost::regex_search(start,end, match1, expression, flags);     
            assert(match1[0].matched);
            std::string s_rnum(match1[1].first,match1[1].second);            
            rnum1 = std::stol(s_rnum);
            isread1 = true;
        }
        if (l2>=0 && !isread2) {
            name2 = std::string(seq2->name.s);
            start = name2.begin();
            end = name2.end();            
            boost::regex_search(start,end, match2, expression, flags);     
            assert(match2[0].matched);
            std::string s_rnum(match2[1].first,match2[1].second);            
            rnum2 = std::stol(s_rnum);
            isread2 = true;
        }
        
        if (isread1 && isread2 && rnum1==rnum2) {
            std::string::iterator kmer_it1 = match1[2].first;
            std::string::iterator kmer_it2 = match2[2].first;
            bool one_found = false;
            if (intersection){
                for(;kmer_it1!=match1[2].second;++kmer_it1,++kmer_it2){
                    *kmer_it1 = std::min(*kmer_it1,*kmer_it2);
                    one_found = *kmer_it1=='1';
                }
            } else {
                for(;kmer_it1!=match1[2].second;++kmer_it1,++kmer_it2){
                    *kmer_it1 = std::max(*kmer_it1,*kmer_it2);
                }
            }
            assert(kmer_it2==match2[2].second);
            
            if (!intersection || one_found){
              fprintf(outf,">");
              fprintf(outf,name1.c_str());
              fprintf(outf,"\n");
              fprintf(outf,seq1->seq.s);
              fprintf(outf,"\n");
            }
            
            l1 = kseq_read(seq1); isread1 = false;
            l2 = kseq_read(seq2); isread2 = false;       
        } else        
        if (isread1 && (!isread2 || rnum1<rnum2)) {
            if (!intersection) {
                fprintf(outf,">");
                fprintf(outf,name1.c_str());
                fprintf(outf,"\n");
                fprintf(outf,seq1->seq.s);
                fprintf(outf,"\n");
            }
            l1 = kseq_read(seq1); isread1 = false;            
        } else
        if (isread2 && (!isread1 || rnum1>rnum2)) {
            if (!intersection) {
                fprintf(outf,">");
                fprintf(outf,name2.c_str());
                fprintf(outf,"\n");
                fprintf(outf,seq2->seq.s);
                fprintf(outf,"\n");
            }
            l2 = kseq_read(seq2); isread2 = false;       
        }       
       
    } while ((l1>=0 && l2>=0) || (!intersection && (l1>=0 || l2>=0) ) );    
    
    kseq_destroy(seq1);kseq_destroy(seq2);
    gzclose(fp1);gzclose(fp2);
    fclose(outf);
    
    return 0;
}