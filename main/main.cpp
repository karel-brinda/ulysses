//
//  main.cpp
//  bloom
//
//  Created by Karel Brinda on 29.04.15.
//  and Maciek Sykulski
//

#include <vector>
#include <getopt.h>
#include "main.hpp"

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
        fprintf(stderr, "Usage:   ulysses create_many [options] <out.dir>\n");
        fprintf(stderr, "\n");
        fprintf(stderr, "Options: -s STR seed [" DEFAULT_SEED "]\n");
        fprintf(stderr, "         -a INT size of the array in bytes [1048576]\n");
        fprintf(stderr, "         -h INT number of hash functions [6]\n");
        fprintf(stderr, "         -i BLM initial bloom filter to OR with (must have hash and size the same)\n");
        fprintf(stderr, "         -r     include also reverse complements of kmers\n");
        fprintf(stderr, "         -t       print statistics of the final BF\n");
        fprintf(stderr, "         -e BLM exclude k-mers present in BLM bloom filter\n");
        fprintf(stderr, "         -n BLM include only k-mers present in BLM bloom filter\n");
        fprintf(stderr, "         -k       shrink to optimal size before saving\n");
        fprintf(stderr, "\n");
        return 1;
    }
    
    //char *fa_fn=argv[optind];
    char *bf_fn=argv[optind]; //Output directory
    
    fprintf(stderr, "Reading seqid_to_taxon map... ");
    read_ID_to_taxon_map(ID_to_taxon_map_filename);
    fprintf(stderr, "DONE.\n");
    
    Bloom initial_bf(0,0,nullptr);
    if (initial_bloom_filename.size()>0){
        initial_bf.load(initial_bloom_filename.c_str());
        assert((initial_bf.nh==(bitvector::size_type)nh)&&(initial_bf.array.size()==(bitvector::size_type)(as_B*8)));
    }
    
    Bloom exclude_bf(0,0,nullptr);
    if (exclude_bloom_filename.size()>0){
        exclude_bf.load(exclude_bloom_filename.c_str());
    }

    Bloom include_bf(0,0,nullptr);
    if (include_bloom_filename.size()>0){
        include_bf.load(include_bloom_filename.c_str());
    }

    std::map<std::string,Bloom> * taxon_bloom_map = bloom_create_many_blooms(
        initial_bloom_filename.size()>0?&initial_bf:nullptr,       
        exclude_bloom_filename.size()>0?&exclude_bf:nullptr,
        include_bloom_filename.size()>0?&include_bf:nullptr,
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
        
    unsigned int bytes_kmer=compressed_kmer_size(bf.seed.weight);
    
    while ((l = kseq_read(seq)) >= 0) { // STEP 4: read sequence
        
        for (int dir=0;dir<=d;dir++){
            printf("%s\t%s\t",seq->name.s, dir ? "r" : "f");
            for (int i=0;i<l-bf.seed.span+1;i++){
                int res=bf.query((uchar*)&(seq->seq.s[i]), dir, bytes_kmer);
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

