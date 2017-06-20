//
//  narrow_TEpeaks.cpp
//  TEToolkit_c++
//
//
//  Created by Ying Jin on 2/22/16.
//  Copyright (c) 2016 Ying Jin. All rights reserved.
//


#include <string>
#include <stdlib.h>

#include <algorithm>
#include <iostream>
#include <fstream>
#include <unistd.h>
#include <ctime>
#include <sstream>

//#include <boost/iostreams/filter/zlib.hpp>

//#include <boost/iostreams/filter/gzip.hpp>
#include <boost/iostreams/filtering_stream.hpp>
#include <boost/iostreams/device/file.hpp>
#include <boost/algorithm/string/predicate.hpp>

//#include "IntervalTree.h"
//#include "TEToolkit/Candidate_Peaks.h"

//#include "TEToolkit/EMestimate_reads.h"
//#include "zeroin.h"
#include "myLog.h"
#include "Parser.h"
//#include "TEToolkit/PeakDetect.h"
//#include "TEToolkit/PeakModel.h"
//#include "ShortRead.h"

#include "htslib/sam.h"

unsigned int tick_time =3000;

/*bool isgz(std::string filename)
{
    //first, check that it exists
    std::ifstream in(filename);
    if(! in)
    {
        return false;
    }
    in.close();
    
    //
     If the file is not .gz and we try
     to read from it using boost's filtering_istream
     classes, an exception will be thrown.
     
     Therefore, catching the exception tells
     you the file is not .gz.
     //
    boost::iostreams::filtering_istream gzin;
    gzin.push(boost::iostreams::gzip_compressor());
    gzin.push(boost::iostreams::file_source(filename),std::ios_base::in | std::ios_base::binary);
    
    char c;
    try
    {
        gzin >> c;
    }
    catch ( boost::iostreams::gzip_error & e )
    {
        gzin.pop();
        return false;
    }
    return true;
}
*/

std::vector<std::pair<int, int> > fetch_exon(int st,uint32_t * cigar,uint32_t n_cigar,std::string format)
{
    //''' fetch exon regions defined by cigar. st must be zero based return list of tuple of (st, end)'''
    //match = re.compile(r'(\d+)(\D)')
    int chrom_st = st;
    if (format == "BAM") { chrom_st += 1;}
    
    std::vector<std::pair<int,int> >exon_bound;
    
    for (unsigned int i=0; i < n_cigar; i++){ //code and size
        if (bam_cigar_op(cigar[i])==BAM_CMATCH) { exon_bound.push_back(std::pair<int, int> (chrom_st,chrom_st + bam_cigar_oplen(cigar[i])-1));} //match
        if (bam_cigar_op(cigar[i])==BAM_CINS) {continue;} //insertion
        if (bam_cigar_op(cigar[i])==BAM_CDEL) { chrom_st += bam_cigar_oplen(cigar[i]);} //deletion
        if (bam_cigar_op(cigar[i])==BAM_CREF_SKIP) { chrom_st += bam_cigar_oplen(cigar[i]); } //gap
        if (bam_cigar_op(cigar[i])==BAM_CSOFT_CLIP) { chrom_st += bam_cigar_oplen(cigar[i]);}// soft clipping
    }
    return exon_bound;
}

int rightMost_pos(int st,uint32_t * cigar,uint32_t n_cigar,std::string format)
{
    //''' fetch exon regions defined by cigar. st must be zero based return list of tuple of (st, end)'''
    //match = re.compile(r'(\d+)(\D)')
    int chrom_st = st;
    if (format == "BAM") { chrom_st += 1;}
    
    //std::vector<std::pair<int,int> >exon_bound;
    
    for (unsigned int i=0; i < n_cigar; i++){ //code and size
        if (bam_cigar_op(cigar[i])==BAM_CMATCH) { chrom_st += bam_cigar_oplen(cigar[i])-1;} //match
        if (bam_cigar_op(cigar[i])==BAM_CINS) {continue;} //insertion
        if (bam_cigar_op(cigar[i])==BAM_CDEL) { chrom_st += bam_cigar_oplen(cigar[i]);} //deletion
        if (bam_cigar_op(cigar[i])==BAM_CREF_SKIP) { chrom_st += bam_cigar_oplen(cigar[i]); } //gap
        if (bam_cigar_op(cigar[i])==BAM_CSOFT_CLIP) { chrom_st += bam_cigar_oplen(cigar[i]);}// soft clipping
    }
    return chrom_st;
}

template <typename T>
void process_aligned_fragment(global_context_t<T> * global_context, thread_context_t<T> * thread_context,read_pair_t<T> * read_pair)
{
    
   // int taglength = 0;
   // int fraglength = 0;
    
    bool multi_read = false;
    
    if (read_pair->first_reads.size() > 1 || read_pair->second_reads.size() >1) {
        multi_read = true;
    }
    
    if (global_context->format=="BAM" || global_context->format == "SAM") {
        
        bam1_t * cur_read1 = NULL;
        bam1_t * cur_read2 = NULL;
       // int readID = read_pair->readID;
        
        size_t alignement_cnt = 0;
        std::vector<std::pair<read_t, read_t> > multi_aligments;
        double w = 1.0;
        
        if (!multi_read ) { //unique read
            if (!global_context->isPE) { //SE
                cur_read1 = (bam1_t *)read_pair->first_reads[0];
                if((size_t)cur_read1->core.tid >= global_context->refnames.size()){
                    
                    error("Missing reference sequences with id " + std::to_string(cur_read1->core.tid));
                    std::exit(1);
                }
                else {
                    //reference sequence
                    std::string chrom = global_context->refnames[cur_read1->core.tid];
                    if (IS_REVERSE(cur_read1)) {
                        uint32_t *cigar = bam_get_cigar(cur_read1);
                        thread_context->track->add_loc(chrom,rightMost_pos(cur_read1->core.pos, cigar, cur_read1->core.n_cigar,global_context->format),1,w);
                        
                    }
                    else {
                        thread_context->track->add_loc(chrom,cur_read1->core.pos,0,w);
                    }
                    thread_context->taglength += cur_read1->core.l_qseq;
                    thread_context->n_tags += 1;
                    thread_context->track->inc_uniq_read();
                    
                }
                
            }
            else {//PE
                if (read_pair->first_reads.size()==1) {
                    cur_read1 = (bam1_t *)read_pair->first_reads[0];
                    
                }
                if (read_pair->second_reads.size()==1) {
                    cur_read2 = (bam1_t *)read_pair->second_reads[0];
                }
                
                if(cur_read1 != NULL && !IS_UNMAPPED(cur_read1))
                {
                    std::string chrom1 = global_context->refnames[cur_read1->core.tid];
                    
                    //second read is not mapped, add first read only
                    if (cur_read2 == NULL || IS_UNMAPPED(cur_read2)) {
                        uint32_t *cigar = bam_get_cigar(cur_read1);
                        thread_context->track->add_loc_pe(chrom1,cur_read1->core.pos,rightMost_pos(cur_read1->core.pos, cigar, cur_read1->core.n_cigar,global_context->format),w);
                    }
                    else {//both mates were mapped
                        std::string chrom2 = global_context->refnames[cur_read2->core.tid];
                        
                        int strand1 = IS_REVERSE(cur_read1)? -1 : 1;
                        int strand2 = IS_REVERSE(cur_read2)? -1 : 1;
                        
                        //same chromosome and different direction (considered proper pair)
                        if (chrom1 == chrom2 && (strand1 * strand2 ==-1) ) {
                            //tracking fragment length
                            thread_context->fraglength += cur_read1->core.isize;
                            thread_context->n_frags += 1;
                            //using the left most position as start
                            int start = std::min(cur_read1->core.pos,cur_read2->core.pos);
                            thread_context->track->add_loc_pe(chrom1,start,start + cur_read1->core.isize,w);
                        }
                        else { //add bad mates separately
                            uint32_t *cigar = bam_get_cigar(cur_read1);
                            
                            thread_context->track->add_loc_pe(chrom1,cur_read1->core.pos,rightMost_pos(cur_read1->core.pos, cigar, cur_read1->core.n_cigar,global_context->format),w*0.5);
                            
                            cigar = bam_get_cigar(cur_read2);
                            
                            thread_context->track->add_loc_pe(chrom2,cur_read2->core.pos,rightMost_pos(cur_read2->core.pos, cigar, cur_read2->core.n_cigar,global_context->format),w*0.5);
                            
                        }

                    }
                    
                                        
                }
                else { //read1 is not mapped
                    std::string chrom = global_context->refnames[cur_read2->core.tid];
                    uint32_t *cigar = bam_get_cigar(cur_read2);
                    
                    thread_context->track->add_loc_pe(chrom,cur_read2->core.pos,rightMost_pos(cur_read2->core.pos, cigar, cur_read2->core.n_cigar,global_context->format),w);
                }
                thread_context->track->inc_uniq_read();
            }
            
        }
        else {
            //multi-read
        //loop over aligments of the first segment
        for (size_t k =0;  k < read_pair->first_reads.size() ; k++) {
            
            alignement_cnt += 1;
            
            read_t algn1;
            read_t algn2;
            
            algn1.chrom = "";
            algn2.chrom = "";
            algn1.start = -1;
            algn2.start = -1;
            algn1.end = -1;
            algn2.end = -1;
            algn1.strand = ".";
            algn2.strand = ".";
            
            //get the first alignment of the first segment
            cur_read1 = (bam1_t *)read_pair->first_reads[k];
            
            if (read_pair->second_reads.size() > 0) {
                //the first alignment of the second segment
                cur_read2 = (bam1_t *)read_pair->second_reads[k];
            }
            
            if(cur_read1 != NULL ) {
                if(!IS_UNMAPPED(cur_read1)){
                    //chomosome appear in the header
                    if((size_t)cur_read1->core.tid >= global_context->refnames.size()){
                    
                        error("Missing reference sequences with id " + std::to_string(cur_read1->core.tid));
                        std::exit(1);
                    }
                    else {
                        //reference sequence
                        algn1.chrom = global_context->refnames[cur_read1->core.tid];
                    }
                    algn1.strand = IS_REVERSE(cur_read1) ? "-" : "+";
                    //left most position
                    algn1.start = cur_read1->core.pos;
                    
                    uint32_t *cigar = bam_get_cigar(cur_read1);
                    //save the right most position
                    algn1.end = rightMost_pos(cur_read1->core.pos, cigar, cur_read1->core.n_cigar,global_context->format);
                    
                    bam_destroy1(cur_read1);

                    //std::cout << "cur read " << readID  << " fpos " << fpos1 << std::endl;
                }
            }
            
            if(cur_read2 != NULL){
                if(! IS_UNMAPPED(cur_read2)){
                    //uint32_t *cigar = bam_get_cigar(cur_read2);
                    
                    algn2.strand =  IS_REVERSE(cur_read2) ? "-" : "+";
                    
                    if((size_t)cur_read2->core.tid >= global_context->refnames.size()){
                        
                        error("Missing reference sequences with id " + std::to_string(cur_read2->core.tid));
                        std::exit(1);
                    }
                    else {
                        algn2.chrom = global_context->refnames[cur_read2->core.tid];
                    }

                    algn2.start = cur_read2->core.pos;
                    uint32_t *cigar = bam_get_cigar(cur_read2);
                    //save the right most position
                    algn2.end = rightMost_pos(cur_read2->core.pos, cigar, cur_read2->core.n_cigar,global_context->format);
                    
                    bam_destroy1(cur_read2);
                }
            }
            
            if (algn1.chrom != "" || algn2.chrom != "") {
                multi_aligments.push_back(std::pair<read_t, read_t> (algn1,algn2));
            }

        }
        
        //there are remaining alignments of the second segment, e.g., only the second segment is mapped
        if (alignement_cnt < read_pair->second_reads.size()) {
            
            read_t algn2;
            algn2.chrom = "";
            algn2.start = -1;
            algn2.end = -1;
            algn2.strand = ".";
            
            for (size_t i = alignement_cnt; i < read_pair->second_reads.size(); i++) {
                
                read_t algn2;
                read_t algn1;
                
                algn2.chrom = "";
                algn2.start = -1;
                algn2.end = -1;
                algn2.strand = ".";
                
                algn1.chrom = "";
                algn1.start = -1;
                algn1.end = -1;
                algn1.strand = ".";
                
                cur_read2 = (bam1_t *)read_pair->second_reads[i];
                
                if((size_t)cur_read2->core.tid >= global_context->refnames.size()){
                        error("Missing reference sequences with id " + std::to_string(cur_read2->core.tid));
                        std::exit(1);
                    }
                algn2.chrom = global_context->refnames[cur_read2->core.tid];
                algn2.start = cur_read2->core.pos;
                uint32_t *cigar = bam_get_cigar(cur_read2);
                //save the right most position
                algn2.end = rightMost_pos(cur_read2->core.pos, cigar, cur_read2->core.n_cigar,global_context->format);
                
                algn2.strand =  IS_REVERSE(cur_read2) ? "-" : "+";
                
                bam_destroy1(cur_read2);
                
                if (algn1.chrom != "" || algn2.chrom != "") {
                    multi_aligments.push_back(std::pair<read_t, read_t> (algn1,algn2));
                }

            }
            
        }
            
            //debug("before inc_multi_read");
        if (multi_aligments.size() >0) {
            thread_context->track->inc_multi_read();
        }
        if (global_context->flag !=MULTICOUNTSONLY ) {
                
            
        w = 1.0 / multi_aligments.size();
            
        for (size_t k =0 ; k < multi_aligments.size(); k++) {
                
                read_t r1 = multi_aligments[k].first;
                read_t r2 = multi_aligments[k].second;
                
                if (r1.chrom != "" && r2.chrom !="") {
                    //proper pairs
                    if (r1.chrom == r2.chrom && (r1.strand != r2.strand)) {
                        int s = std::min(r1.start,r2.start);
                        int e = std::max(r1.end,r2.end);
                        
                        thread_context->track->add_loc_pe(r1.chrom,s,e,w);
                    }
                    else { //bad pairs add individual segments
                        thread_context->track->add_loc_pe(r1.chrom,r1.start,r1.end,w*0.5);
                        thread_context->track->add_loc_pe(r2.chrom,r2.start,r2.end,w*0.5);
                    }
                    
                }
                else{
                    if (global_context->isPE) { //PE
                        if (r1.chrom != "") {
                            thread_context->track->add_loc_pe(r1.chrom,r1.start,r1.end,w);
                        }
                        if (r2.chrom != "") {
                            thread_context->track->add_loc_pe(r2.chrom,r2.start,r2.end,w);
                        }
                    }
                    else {//SE
                        //debug("before add loc " + std::to_string(w));
                        if (r1.strand == "-") {
                            thread_context->track->add_loc(r1.chrom,r1.end,1,w);
                        }
                        else {
                            thread_context->track->add_loc(r1.chrom,r1.start,0,w);
                        }
                    }
                }
                
            }
            
        }
        
        }
        delete read_pair;
    }
    else { // BED/BEDPE format
        read_t * cur_read1 = NULL;
      //  read_t * cur_read2 = NULL;
       // int readID = read_pair->readID;
        
        //size_t alignement_cnt = 0;
        
        double w = 1.0 ;// weight
        
        if (multi_read) {
            w = 1.0/read_pair->first_reads.size();
            
            thread_context->track->inc_multi_read();
        }
        else {
            thread_context->track->inc_uniq_read();
        }
        
        for (size_t k =0;  k < read_pair->first_reads.size() ; k++) {
            
            cur_read1 = (read_t *)read_pair->first_reads[k];
            
            if(cur_read1 != NULL ) {
                
                if (global_context->isPE) { //PE reads
                    if (!(multi_read && MULTICOUNTSONLY)) {
                        thread_context->track->add_loc_pe(cur_read1->chrom,cur_read1->start,cur_read1->end,w);
                    }
                }
                else {//SE reads
                    if(cur_read1->strand == "+" || cur_read1->strand == ".")
                    {
                        if (!(multi_read && MULTICOUNTSONLY)) {
                            thread_context->track->add_loc(cur_read1->chrom,cur_read1->start,0,w);
                        }
                    }
                    else {
                        if (!(multi_read && MULTICOUNTSONLY)) {
                            thread_context->track->add_loc(cur_read1->chrom,cur_read1->end,1,w);
                        }
                    }
                }
                if (!multi_read) { //track tag length or fragment length using uniquely mapped reads only
                    if (global_context->isPE) { //fragment length
                        thread_context->fraglength += cur_read1->end - cur_read1->start;
                        thread_context->n_frags += 1;
                    }
                    else {//tag length
                        thread_context->taglength += cur_read1->end - cur_read1->start;
                        thread_context->n_tags += 1;
                    }
                }
                
                delete cur_read1;
                
                //std::cout << "cur read " << readID  << " pid " << pid1 << std::endl;
            }
            
            
        }
        
        delete read_pair;
        
    }
}

//parse multi-reads only
template <typename T>
void process_aligned_fragment_multi(global_context_t<T> * global_context,thread_context_t<T> * thread_context,read_pair_t<T> * read_pair)
{
    //debug("in process_aligned_fragment_multi");
    if (global_context->format=="BAM" || global_context->format == "SAM") {
        
        bam1_t * cur_read1 = NULL;
        bam1_t * cur_read2 = NULL;
        int readID = read_pair->readID;
        std::string strand1 = ".";
        std::string strand2 = ".";
        std::vector<std::pair<int, int> > exons ;
        //std::string qname = "";
        std::string chrom1 = "";
        std::string chrom2 = "";
        
        size_t alignement_cnt = 0;
        bool multi_read = false;
        
        if(read_pair->first_reads.size() >1 || read_pair->second_reads.size() > 1)
        {
            multi_read = true;
            //debug("found a multi-read " + std::to_string(readID));
        }
        
        for (size_t k =0;  k < read_pair->first_reads.size() ; k++) {
            int chrom1_id = -1;
            int chrom2_id = -1;
            alignement_cnt += 1;
            
            
            cur_read1 = (bam1_t *)read_pair->first_reads[k];
            
            if (read_pair->second_reads.size() > 0) {
                cur_read2 = (bam1_t *)read_pair->second_reads[k];
            }
            
            int pid1 = -1;
            int pid2 = -1;
            int pid = -1;
            
            if(cur_read1 != NULL ) {
                if(!IS_UNMAPPED(cur_read1)){
                    uint32_t *cigar = bam_get_cigar(cur_read1);
                    //char * name = bam_get_qname(cur_read1);
                    chrom1_id  =  cur_read1->core.tid;
                    
                    //qname = std::string(name);
                    
                    if((size_t)cur_read1->core.tid < global_context->refnames.size()){
                        chrom1 = global_context->refnames[cur_read1->core.tid];
                    }
                    else {
                        std::cout << (size_t)cur_read1->core.tid  << "\t" << "missing reference sequences." << std::endl;
                        std::exit(1);
                    }
                    
                    
                    strand1 = IS_REVERSE(cur_read1) ? "-" : "+";
                    
                    exons = fetch_exon(cur_read1->core.pos, cigar,cur_read1->core.n_cigar, global_context->format);
                    
                    pid1 = global_context->peakIdx->get_ovp_peaks(chrom1,exons,global_context->shiftsize,strand1);
                    
                    bam_destroy1(cur_read1);
                    
                    
                    //std::cout << "cur read " << readID  << " pid " << pid1 << std::endl;
                }
            }
            
            if(cur_read2 != NULL){
                if(! IS_UNMAPPED(cur_read2)){
                    uint32_t *cigar = bam_get_cigar(cur_read2);
                    
                    strand2 =  IS_REVERSE(cur_read2) ? "-" : "+";
                    
                    chrom2_id = cur_read2->core.tid;
                    
                    if((size_t)cur_read2->core.tid < global_context->refnames.size()){
                        chrom2 = global_context->refnames[cur_read2->core.tid];
                    }
                    else {
                        std::cout << (size_t)cur_read2->core.tid << "\t" <<"missing reference sequences." << std::endl;
                        std::exit(1);
                    }
                    exons = fetch_exon(cur_read2->core.pos, cigar,cur_read2->core.n_cigar, global_context->format);
                    
                    pid2 = global_context->peakIdx->get_ovp_peaks(chrom2,exons,global_context->shiftsize,strand2);
                    bam_destroy1(cur_read2);
                }
            }
            
            if (chrom1_id != -1 && chrom2_id != -1){
                
                if (chrom1_id != chrom2_id ){
                    //bam_destroy1(cur_read1);
                    //bam_destroy1(cur_read2);
                    
                    continue;
                }
                
            }
            
            if(pid1 == -1 && pid2  != -1)
            {
                pid = pid2;
            }
            if(pid2 == -1 && pid1 != -1)
            {
                pid = pid1;
            }
            if (multi_read) {
                if (pid != -1) {
                    std::map<int,std::vector<int> >::iterator itr;
                    itr = thread_context->multi_read_mapTo_pids->find(readID);
                    
                    if (itr != thread_context->multi_read_mapTo_pids->end()) {
                        (* thread_context->multi_read_mapTo_pids)[readID].push_back(pid);
                    }
                    else {
                        std::vector<int> pidlist ;
                        pidlist.push_back(pid);
                        thread_context->multi_read_mapTo_pids->insert(std::pair<int,std::vector<int> > (readID,pidlist));
                    }
                }
            }
            else { //uniq read
                if (pid != -1 ) {
                    //(*thread_context->peak_reads)[pid] += 1;
                    //std::cout << "uniq peak read" << std::endl;
                    thread_context->peak_reads->at(pid) += 1.0;
                }
            }
        }
        
        if (alignement_cnt < read_pair->second_reads.size()) {
            for (size_t i = alignement_cnt; i < read_pair->second_reads.size(); i++) {
                cur_read2 = (bam1_t *)read_pair->second_reads[i];
                int pid = -1;
                //int chrom2_id = -1;
                if(! IS_UNMAPPED(cur_read2)){
                    uint32_t *cigar = bam_get_cigar(cur_read2);
                    
                    strand2 =  IS_REVERSE(cur_read2) ? "-" : "+";
                    
                    //chrom2_id = cur_read2->core.tid;
                    
                    if((size_t)cur_read2->core.tid < global_context->refnames.size()){
                        chrom2 = global_context->refnames[cur_read2->core.tid];
                    }
                    else {
                        std::cout << (size_t)cur_read2->core.tid << "\t" <<"missing reference sequences." << std::endl;
                        std::exit(1);
                    }
                    exons = fetch_exon(cur_read2->core.pos, cigar,cur_read2->core.n_cigar, global_context->format);
                    
                    pid = global_context->peakIdx->get_ovp_peaks(chrom2,exons,global_context->shiftsize,strand2);
                    bam_destroy1(cur_read2);
                }
                
                if (pid != -1) {
                    if (multi_read) {
                        std::map<int,std::vector<int> >::iterator itr;
                        itr = thread_context->multi_read_mapTo_pids->find(readID);
                        
                        if (itr != thread_context->multi_read_mapTo_pids->end()) {
                            (* thread_context->multi_read_mapTo_pids)[readID].push_back(pid);
                        }
                        else {
                            std::vector<int> pidlist ;
                            pidlist.push_back(pid);
                            thread_context->multi_read_mapTo_pids->insert(std::pair<int,std::vector<int> > (readID,pidlist));
                        }
                        
                    }
                    else { //uniq read
                        
                        //    (*thread_context->peak_reads)[pid] += 1;
                        thread_context->peak_reads->at(pid) += 1.0;
                    }
                    
                }
                
            }
        }
        
        delete read_pair;
        
    }
    else {
        read_t * cur_read1 = NULL;
        read_t * cur_read2 = NULL;
        int readID = read_pair->readID;
        
        size_t alignement_cnt = 0;
        bool multi_read = false;
        
        if(read_pair->first_reads.size() >1 || read_pair->second_reads.size() > 1)
        {
            multi_read = true;
        }
        
        for (size_t k =0;  k < read_pair->first_reads.size() ; k++) {
            //int chrom1_id = -1;
            //int chrom2_id = -1;
            alignement_cnt += 1;
            cur_read1 = (read_t *)read_pair->first_reads[k];
            
            if (read_pair->second_reads.size() > 0) {
                cur_read2 = (read_t *)read_pair->second_reads[k];
            }
            
            int pid1 = -1;
            int pid2 = -1;
            int pid = -1;
            
            if(cur_read1 != NULL ) {
                
                std::vector<std::pair<int, int> > itv_list;
                
                itv_list.push_back(std::pair<int,int> (cur_read1->start,cur_read1->end));
                
                pid1 = global_context->peakIdx->get_ovp_peaks(cur_read1->chrom,itv_list,global_context->shiftsize,cur_read1->strand);
                
                delete cur_read1;
                
                //std::cout << "cur read " << readID  << " pid " << pid1 << std::endl;
            }
            
            if(cur_read2 != NULL){
                std::vector<std::pair<int, int> > itv_list;
                
                itv_list.push_back(std::pair<int,int> (cur_read2->start,cur_read2->end));
                
                pid1 = global_context->peakIdx->get_ovp_peaks(cur_read2->chrom,itv_list,global_context->shiftsize,cur_read2->strand);
                
                delete cur_read2;
                
            }
            
            if (pid1 != -1 && pid2 != -1){
                
                if (pid1 != pid2 ){//unproper pair
                    
                    continue;
                }
                
            }
            
            if(pid1 == -1 && pid2  != -1)
            {
                pid = pid2;
            }
            if(pid2 == -1 && pid1 != -1)
            {
                pid = pid1;
            }
            if (multi_read) {
                if (pid != -1) {
                    std::map<int,std::vector<int> >::iterator itr;
                    itr = thread_context->multi_read_mapTo_pids->find(readID);
                    
                    if (itr != thread_context->multi_read_mapTo_pids->end()) {
                        (* thread_context->multi_read_mapTo_pids)[readID].push_back(pid);
                    }
                    else {
                        std::vector<int> pidlist ;
                        pidlist.push_back(pid);
                        thread_context->multi_read_mapTo_pids->insert(std::pair<int,std::vector<int> > (readID,pidlist));
                    }
                }
            }
            else { //uniq read
                if (pid != -1 ) {
                    thread_context->peak_reads->at(pid) += 1.0;
                }
            }
        }
        
        if (alignement_cnt < read_pair->second_reads.size()) {
            for (size_t i = alignement_cnt; i < read_pair->second_reads.size(); i++) {
                cur_read2 = (read_t *)read_pair->second_reads[i];
                int pid = -1;
                //int chrom2_id = -1;
                if(cur_read2 != NULL){
                    std::vector<std::pair<int, int> > itv_list;
                    
                    itv_list.push_back(std::pair<int,int> (cur_read2->start,cur_read2->end));
                    
                    pid = global_context->peakIdx->get_ovp_peaks(cur_read2->chrom,itv_list,global_context->shiftsize,cur_read2->strand);
                    
                    delete cur_read2;
                }
                
                if (pid != -1) {
                    if (multi_read) {
                        std::map<int,std::vector<int> >::iterator itr;
                        itr = thread_context->multi_read_mapTo_pids->find(readID);
                        
                        if (itr != thread_context->multi_read_mapTo_pids->end()) {
                            (* thread_context->multi_read_mapTo_pids)[readID].push_back(pid);
                        }
                        else {
                            std::vector<int> pidlist ;
                            pidlist.push_back(pid);
                            thread_context->multi_read_mapTo_pids->insert(std::pair<int,std::vector<int> > (readID,pidlist));
                        }
                        
                    }
                    else { //uniq read
                        thread_context->peak_reads->at(pid) += 1.0;
                    }
                    
                }
                
            }
        }
        
        delete read_pair;
        
    }
}



template <typename T>
void* worker(void * vargs)
{
    struct arg_struct<T> * args = (struct arg_struct<T> *) vargs;
    thread_context_t<T> * thread_context = args->arg2;
    global_context_t<T> * global_context = args->arg1;
    delete args;
    
    while (1){
        //std::pair<bam1_t *, bam1_t *> cur_read_pair;
        read_pair_t<T>  * cur_read_pair;
        while(1){
            int is_retrieved = 0;
            pthread_spin_lock(&thread_context->cur_reads_lock);
            
            if(thread_context->cur_reads->size()>0){
                cur_read_pair = thread_context->cur_reads->back();
                thread_context->cur_reads->pop_back();
                
                is_retrieved = 1;
            }
            pthread_spin_unlock(&thread_context->cur_reads_lock);
            if(global_context->all_finished && !is_retrieved) return NULL;
            
            if(is_retrieved) break;
            else usleep(tick_time);
        }
        if (global_context->flag == MULTIONLY) {
            process_aligned_fragment_multi<T>(global_context,thread_context,cur_read_pair);
        }
        else {
            process_aligned_fragment<T>(global_context,thread_context,cur_read_pair);
        }
        
        
    }
}


//start threads
template <typename T>
void init_thread(global_context_t<T> * global_context,unsigned short threadNumber)
{
    global_context->thread_number = threadNumber;
    
    int numOfpeaks = 0;
    if (global_context->peakIdx != nullptr) {
        numOfpeaks = global_context->peakIdx->get_numofpeaks();
    }
   
    
    if(threadNumber >1){
        for(int i=0;i<threadNumber;i++)
        {
            thread_context_t<T> *th_contx = new thread_context_t<T>();
            pthread_spin_init(&th_contx->cur_reads_lock, PTHREAD_PROCESS_PRIVATE);
            th_contx->thread_id = i;
            // th_contx->cur_reads = new std::vector<std::pair<bam1_t *,bam1_t *> >() ;
            
            th_contx->cur_reads = new std::vector<read_pair_t<T> *> ();
            th_contx->track = new ShortRead();
            
            if (numOfpeaks > 0) {
                //debug("init threads init peak_reads");
                th_contx->peak_reads = new std::vector<int>();
            
                for (int j=0; j < numOfpeaks; j++) {
                    th_contx->peak_reads->push_back(0);
                }
                th_contx->multi_read_mapTo_pids = new MULTI_READ_MAP ();
            }
        
            global_context->thread_contexts.push_back(th_contx);
        }
        
        for(int i=0;i<threadNumber;i++)
        {
            struct arg_struct<T> * args = new arg_struct<T>();
            args->arg1 = global_context;
            args->arg2 = global_context->thread_contexts[i];
            //int ret = pthread_create(&global_context -> thread_contexts[i].thread_object, NULL, worker, (void *)&args);
            pthread_create(&global_context -> thread_contexts[i]->thread_object, NULL, worker<T>, (void *)args);
        }
    }
    else {
        
        thread_context_t<T> *th_contx = new thread_context_t<T>();
        //th_contx->cur_reads = new std::vector<std::pair<bam1_t *,bam1_t *> >() ;
        th_contx->cur_reads = new std::vector<read_pair_t<T> *> ();
        
        if (numOfpeaks > 0) {
            //debug("init threads init peak_reads");
            th_contx->peak_reads = new std::vector<int>();
            
            for (int j=0; j < numOfpeaks; j++) {
                th_contx->peak_reads->push_back(0);
            }
            th_contx->multi_read_mapTo_pids = new MULTI_READ_MAP ();
        }
        th_contx->track = new ShortRead();
        
        global_context->thread_contexts.push_back(th_contx);
        
    }
    
    return;
}

//release memory
template <typename T>
void destroy_thread_context(global_context_t<T> * global_context)
{
    
    for(int i=0; i < global_context-> thread_number; i++)
    {
        
        if (global_context->thread_contexts[i]->multi_read_mapTo_pids != nullptr) {
            delete global_context->thread_contexts[i]->multi_read_mapTo_pids;
        }
        if (global_context->thread_contexts[i]->peak_reads != nullptr) {
            delete global_context->thread_contexts[i]->peak_reads;
        }
        
        delete global_context->thread_contexts[i]->cur_reads;
        delete global_context->thread_contexts[i]->track;
        
        if(global_context->thread_number >1){
            pthread_spin_destroy(&global_context -> thread_contexts[i]->cur_reads_lock);
        }
        delete global_context -> thread_contexts[i];
    }
    
}

//join threads
template <typename T>
void join_threads(global_context_t<T> * global_context)
{
    int xk1;
    for(xk1=0; xk1<global_context-> thread_number; xk1++)
        pthread_join(global_context -> thread_contexts[xk1]->thread_object, NULL);
}

//merge results
template <typename T>
void merge_results(global_context_t<T> * global_context, ShortRead * track)
{
    long fraglength = 0;
    long taglength = 0;
    int n_frags = 0;
    int n_tags = 0;
    
    for(size_t i=0;i < global_context->thread_number;i++)
    {
        track->append_shortreads(global_context->thread_contexts[i]->track);
        if (global_context->isPE) {
            fraglength += global_context->thread_contexts[i]->fraglength;
            n_frags += global_context->thread_contexts[i]->n_frags;
        }
        else {
            taglength += global_context->thread_contexts[i]->taglength;
            n_tags += global_context->thread_contexts[i]->n_tags;
        }
    }
    if (global_context->isPE) {
        if (n_frags == 0) {
            error("Error in estimating fragment length : not enough number of fragments.");
            std::exit(1);
        }
        fraglength = fraglength/n_frags;
        track->set_fraglength(int(fraglength));
    }
    else {
        if (n_tags == 0) {
            error("Error in estimating fragment length : not enough number of fragments.");
            std::exit(1);
        }
        taglength = taglength/n_tags;
        track->set_taglength(int(taglength));
    }
}

template <typename T>
void merge_results_redistr(global_context_t<T> * global_context, std::vector<double> * peak_reads, MULTI_READ_MAP * multi_read_mapTo_pids)
{
    
    
    for(size_t i=0;i < global_context->thread_number;i++)
    {
        std::map<int, std::vector<int>>::iterator itr ;
        for (size_t j=0;j < peak_reads->size();j++) {
            
            peak_reads->at(j) += global_context->thread_contexts[i]->peak_reads->at(j);
        }
        
        for (itr = global_context->thread_contexts[i]->multi_read_mapTo_pids->begin(); itr != global_context->thread_contexts[i]->multi_read_mapTo_pids->end(); itr ++)
        {
            multi_read_mapTo_pids->insert(std::pair<int,std::vector<int> > (itr->first,itr->second));
        }
    }
}

//read alignment files. re-distribute multi-reads among candidate peaks if the last three params are not nullptr
void parse_bam(ShortRead * track,std::string inputFile,int shiftsize,int threadNumber,int flag,Candidate_Peaks * peakIdx,std::vector<double> * peak_reads_Prime, MULTI_READ_MAP * multi_read_mapTo_pids)
{
    
    //int q_cut = mapq;
    
    global_context_t<bam1_t> global_context;
    global_context.all_finished = 0;
    global_context.flag = flag;
    
    if (flag == MULTIONLY) {
        global_context.shiftsize = shiftsize;
        global_context.peakIdx = peakIdx;
    }
    else {
        global_context.peakIdx = nullptr;
        
    }
    
    
    //global_context.shiftsize = shiftsize;
    
    //std::vector<double> peak_reads(peak_reads_Prime.size(),0.0);
    
    //open SAM/BAM file
    samFile *fp_in = NULL;
    bam1_t * aligned_read = NULL;
    bam_hdr_t *header = NULL;
    
    std::string format = "BAM";
    
    
    
    fp_in = sam_open(inputFile.c_str(), "rb");
    
    if(NULL == fp_in) {
        
        error("Could not open file " + inputFile);
        std::exit(1);
    }
    else {
        info("open file " + inputFile);
    }
    
    if(hts_get_format(fp_in)->format == sam) {format = "SAM";}
    
    global_context.format = format;
    global_context.isPE = true;
    
    aligned_read = bam_init1();
    
    //ShortRead * track = new ShortRead();
    
    //Results main_res;
    //read header
    header = sam_hdr_read(fp_in);
    
    if(NULL == header){
        //std::cout << "No header information." << std::endl;
        error("No header information in " + inputFile);
        std::exit(1);
    }
    for(int i=0;i< header->n_targets;i++)
    {
        char * p = header->target_name[i];
        int  l = header->target_len[i];
        std::string chr(p);
        
        global_context.refnames.push_back(chr);
        global_context.reflengths.insert(std::pair<std::string, int>(chr,l));
    }
    
   // debug("read in header information. number of reference sequence is " + std::to_string(global_context.refnames.size()));
    
    std::string prev_read_id = "";
    std::vector<bam1_t *> multi_read1;
    std::vector<bam1_t *> multi_read2;
    
    std::vector<bam1_t *> alignments_per_read;
    
    int lineno = 0;
    int readID = 0;
    
    //std::map<int,std::vector<int> > multi_read_mapTo_pids;
    
    
    init_thread<bam1_t>(&global_context, threadNumber);
    
    int current_thread_id = 0;
    
    
    
    while(sam_read1(fp_in, header,aligned_read) > 0)
    {
        lineno +=1 ;
        if (lineno % 1000000 == 0 ) {
           // const time_t ctt = time(0);
           // std::string cur_time = asctime(localtime(&ctt));
           // cur_time.erase(std::remove(cur_time.begin(), cur_time.end(), '\n'), cur_time.end());
           // std::cout << "INFO  @ " << cur_time << "\t" << lineno << " aligments." << std::endl;
            info(std::to_string(lineno) + " alignments.");
        }
        //std::cout << lineno << std::endl;
        
        std::string qname ="";
        char * name = bam_get_qname(aligned_read);
        if(name !=NULL)
        {
            qname = std::string(name);
        }
        //SE reads
        if (! IS_PAIRED(aligned_read))
        {
            //debug("is single end read.");
            
            global_context.isPE = false;
            //count unmapped read
            if (IS_UNMAPPED(aligned_read) ||IS_DUP(aligned_read)|| aligned_read->core.flag & BAM_FQCFAIL )
            {
                bam_destroy1(aligned_read);
                aligned_read =  bam_init1();
                continue;
            }
            
            if (qname == prev_read_id || prev_read_id == "") {
                alignments_per_read.push_back(aligned_read);
                prev_read_id = qname;
                aligned_read =  bam_init1();
                //debug("multi reads.");
                continue;
            }
            else {
                //debug("new read");
                bam1_t *cur_read = NULL;
                if (alignments_per_read.size() == 1){  //unique read
                    
                   // debug("uniq read");
                    cur_read = alignments_per_read[0];

                    if (cur_read  ){
                        readID += 1;
                        read_pair_t<bam1_t> * read_pair = new read_pair_t<bam1_t>();
                        read_pair->first_reads.push_back(cur_read);
                        read_pair->readID = readID;
                        if(global_context.thread_number >1){
                            
                            thread_context_t<bam1_t> * thread_context = global_context.thread_contexts[current_thread_id];
                            pthread_spin_lock(&thread_context->cur_reads_lock);
                            
                            //thread_context->cur_reads->push_back(std::pair<bam1_t *,bam1_t *>(cur_read,NULL));
                            thread_context->cur_reads->push_back(read_pair);
                            
                            pthread_spin_unlock(&thread_context->cur_reads_lock);
                            current_thread_id++;
                            if(current_thread_id >= global_context.thread_number) current_thread_id = 0;
                        }
                        else {
                            if (flag == MULTIONLY) {
                                process_aligned_fragment_multi<bam1_t>(&global_context,global_context.thread_contexts[0],read_pair);
                            }
                            else {
                               // debug("process_aligned_fragment");
                                process_aligned_fragment<bam1_t>(&global_context,global_context.thread_contexts[0],read_pair);
                            }
                        }
                    }
                    else {
                        bam_destroy1(aligned_read);
                    }
                    
                }
                else { //multi reads
                    //cur_read = alignments_per_read[0];
                   // debug("multi read");
                    if (flag != UNIQONLY) {
                        readID += 1;
                        read_pair_t<bam1_t> * read_pair = new read_pair_t<bam1_t>();
                    
                        for(unsigned int i=0;i<alignments_per_read.size();i++){
                            (read_pair->first_reads).push_back(alignments_per_read[i]);
                            //bam_destroy1(alignments_per_read[i]);
                        }
                        read_pair->readID = readID;
                    
                        if(global_context.thread_number >1){
                            //thread_context_t * thread_context = global_context.thread_contexts + current_thread_id;
                            thread_context_t<bam1_t> * thread_context = global_context.thread_contexts[current_thread_id];
                            pthread_spin_lock(&thread_context->cur_reads_lock);
                        
                            thread_context->cur_reads->push_back(read_pair);
                        //thread_context->cur_reads->push_back(std::pair<bam1_t *,bam1_t *>(cur_read,NULL));
                            pthread_spin_unlock(&thread_context->cur_reads_lock);
                            current_thread_id++;
                            if(current_thread_id >= global_context.thread_number) current_thread_id = 0;
                        }
                        else {
                        //process_aligned_fragment(&global_context,global_context.thread_contexts[0],std::pair<bam1_t *,bam1_t *> (cur_read,NULL));
                            if (flag == MULTIONLY) {
                                process_aligned_fragment_multi<bam1_t>(&global_context,global_context.thread_contexts[0],read_pair);
                            }
                            else {
                            //debug("process_aligned_fragment");
                            process_aligned_fragment<bam1_t>(&global_context,global_context.thread_contexts[0],read_pair);
                            }
                        }
                    }
                    else {
                        for(unsigned int i=0;i<alignments_per_read.size();i++){
                            bam_destroy1(alignments_per_read[i]);
                        }
                        
                    }
                }
                alignments_per_read.clear();
                alignments_per_read.push_back(aligned_read);
                prev_read_id = qname;
                aligned_read = bam_init1();
            }
        }
        else {//pair end read, need to keep tracking of fragment length
            size_t flag_pos = qname.find('/');
            //cur_read_id = aligned_read.qname;
            if (flag_pos != std::string::npos) {qname = qname.substr(0,flag_pos);}
            else {
                if (boost::algorithm::ends_with(qname,".1")|| boost::algorithm::ends_with(qname,".2")) {
                    flag_pos = qname.size() - 2;
                    qname = qname.substr(0,flag_pos);
                }
            }
            
            if (qname == prev_read_id ){
                if (IS_READ1(aligned_read) ){ multi_read1.push_back(aligned_read);}
                if (IS_READ2(aligned_read) ){multi_read2.push_back(aligned_read);}
                aligned_read = bam_init1();
                continue;
            }
            else {
                bam1_t * cur_read1 = NULL;
                bam1_t * cur_read2 = NULL;
                
                //multi-reads
                if (multi_read1.size() >1 || multi_read2.size() > 1){
                if (flag == UNIQONLY) {
                    for(unsigned int i=0;i< multi_read1.size();i++)
                    {
                        bam_destroy1(multi_read1[i]);
                        
                    }
                    for(unsigned int i=0;i< multi_read2.size();i++)
                    {
                        bam_destroy1(multi_read2[i]);
                    }
                    
                }
                else {
                if (multi_read1.size() >1 || multi_read2.size() > 1){
                    readID += 1;
                    read_pair_t<bam1_t> * read_pair = new read_pair_t<bam1_t>();
                    read_pair->readID = readID;
                    
                    for(unsigned int i=0;i< multi_read1.size();i++)
                    {
                        //bam_destroy1(multi_read1[i]);
                        read_pair->first_reads.push_back(multi_read1[i]);
                    }
                    for(unsigned int i=0;i< multi_read2.size();i++)
                    {
                        //bam_destroy1(multi_read2[i]);
                        read_pair->second_reads.push_back(multi_read2[i]);
                    }
                    
                    if(global_context.thread_number > 1){
                        //thread_context_t * thread_context = global_context.thread_contexts + current_thread_id;
                        thread_context_t<bam1_t> * thread_context = global_context.thread_contexts[current_thread_id];
                        
                        pthread_spin_lock(&thread_context->cur_reads_lock);
                        
                        //thread_context->cur_reads->push_back(std::pair<bam1_t *,bam1_t *>(cur_read1,cur_read2));
                        thread_context->cur_reads->push_back(read_pair);
                        
                        pthread_spin_unlock(&thread_context->cur_reads_lock);
                        current_thread_id++;
                        if(current_thread_id >= global_context.thread_number)
                            current_thread_id = 0;
                    }
                    else {
                        if (flag == MULTIONLY) {
                            process_aligned_fragment_multi<bam1_t>(&global_context,global_context.thread_contexts[0],read_pair);
                        }
                        else {
                        process_aligned_fragment<bam1_t>(&global_context,global_context.thread_contexts[0],read_pair);
                        }
                    }
                }
                }
                }
                else {
                    //uniq read
                    if (multi_read1.size() == 1) {cur_read1 = multi_read1[0];}
                    if (multi_read2.size() == 1){ cur_read2 = multi_read2[0];}
                    
                    if (cur_read1 && (IS_UNMAPPED(cur_read1) || (cur_read1->core.flag & BAM_FQCFAIL))){
                        bam_destroy1(cur_read1);
                        cur_read1 = NULL;
                    }
                    
                    if (cur_read2 && (IS_UNMAPPED(cur_read2) || (cur_read2->core.flag & BAM_FQCFAIL))){
                        bam_destroy1(cur_read2);
                        cur_read2 = NULL;
                    }
                    
                    if((cur_read1  || cur_read2) && flag != MULTIONLY){
                        readID += 1;
                        read_pair_t<bam1_t> * read_pair = new read_pair_t<bam1_t>();
                        if (cur_read1 != NULL) {
                            read_pair->first_reads.push_back(cur_read1);
                        }
                        if(cur_read2 != NULL) {
                            read_pair->second_reads.push_back(cur_read2);
                        }
                        read_pair->readID = readID;
                        
                        if(global_context.thread_number > 1){
                            //thread_context_t * thread_context = global_context.thread_contexts + current_thread_id;
                            thread_context_t<bam1_t> * thread_context = global_context.thread_contexts[current_thread_id];
                            
                            pthread_spin_lock(&thread_context->cur_reads_lock);
                            
                            //thread_context->cur_reads->push_back(std::pair<bam1_t *,bam1_t *>(cur_read1,cur_read2));
                            thread_context->cur_reads->push_back(read_pair);
                            
                            pthread_spin_unlock(&thread_context->cur_reads_lock);
                            current_thread_id++;
                            if(current_thread_id >= global_context.thread_number)
                                current_thread_id = 0;
                        }
                        else {
                            //process_aligned_fragment(&global_context,global_context.thread_contexts[0],std::pair<bam1_t *,bam1_t *> (cur_read1,cur_read2));
                            if (flag == MULTIONLY) {
                                process_aligned_fragment_multi<bam1_t>(&global_context,global_context.thread_contexts[0],read_pair);
                            }
                            else {
                            process_aligned_fragment<bam1_t>(&global_context,global_context.thread_contexts[0],read_pair);
                            }
                        }
                    }
                }
                multi_read1.clear();
                multi_read2.clear();
                prev_read_id = qname;
                if (IS_READ1(aligned_read)){
                    multi_read1.push_back(aligned_read);}
                if (IS_READ2(aligned_read)){
                    multi_read2.push_back(aligned_read);}
            }
            
            aligned_read =  bam_init1();
        }
    }
    
    
    bam_destroy1(aligned_read);
    bam_hdr_destroy(header);
    sam_close(fp_in);
    
    global_context.all_finished = 1;
    //join threads
    if(global_context.thread_number >1)
    {
        join_threads<bam1_t>(&global_context);
        
    }
    if (flag != MULTIONLY) {
        
        merge_results<bam1_t>(&global_context, track);
    }
    else {
        
        merge_results_redistr<bam1_t>(&global_context,peak_reads_Prime,multi_read_mapTo_pids);
    }
    
    
    
    destroy_thread_context<bam1_t>(&global_context);
    
    if (flag != MULTIONLY) {
        track->isPE = global_context.isPE;
        track->set_rlengths(global_context.reflengths);
    }
    
    
}


//re-distribute multi-reads in BAM format among candidate peaks.
void parse_bed(ShortRead * track,std::string inputFile,int shiftsize,int threadNumber,std::string format,int flag,Candidate_Peaks * peakIdx,std::vector<double> * peak_reads_Prime, MULTI_READ_MAP * multi_read_mapTo_pids)
{
    //'''This class provides fuctions to parsing BED files and quality controls.'''
    
    //int q_cut = mapq;
    global_context_t<read_t> global_context;
    global_context.all_finished = 0;
    global_context.shiftsize = shiftsize;
    global_context.flag = flag;
    
    if (flag == MULTIONLY) {
        global_context.shiftsize = shiftsize;
        global_context.peakIdx = peakIdx;
    }
    
    global_context.isPE = false;
    
    if (format == "BEDPE") {
        global_context.isPE = true;
    }
    
    //std::vector<double> peak_reads(peak_reads_Prime.size(),0.0);
    
    std::string prev_read_id = "";
    std::vector<read_t *> multi_read1;
    std::vector<read_t *> multi_read2;
    
    std::vector<read_t *> alignments_per_read;
    
    int lineno = 0;
    int readID = 0;
    
    //global_context.peakIdx = peakIdx;
    
    //std::map<int,std::vector<int> > multi_read_mapTo_pids;
    
    init_thread<read_t>(&global_context, threadNumber);
    
    
    int current_thread_id = 0;
    std::ifstream input;
   // boost::iostreams::filtering_istream in;
    
    try{
        //bool gzipped = isgz(inputFile);
        bool gzipped = false;
        std::string line;
        
        if(gzipped){
            
     //       input.open(inputFile, std::ios_base::in | std::ios_base::binary);
     //       try {
                
     //           in.push(boost::iostreams::gzip_decompressor());
     //           in.push(input);
     //
      //      }
      //      catch(const boost::iostreams::gzip_error& e) {
                //std::cout << e.what() << '\n';
      //          error(e.what());
      //      }
        
        }
        else {
            
            input.open (inputFile, std::ifstream::in);
        }
        
        while(! input.eof()){
            
            std::string qname;
            
            if (! gzipped && ! std::getline(input,line)){
                break;
            }
        //    if (gzipped && ! std::getline(in,line)) {
         //       break;
         //   }
            lineno +=1 ;
            if (lineno % 1000000 == 0 ) {
                //const time_t ctt = time(0);
                //std::string cur_time = asctime(localtime(&ctt));
                //cur_time.erase(std::remove(cur_time.begin(), cur_time.end(), '\n'), cur_time.end());
                //std::cout << "INFO  @ " << cur_time << "\t" << lineno << " aligments." << std::endl;
                info(std::to_string(lineno) + " alignments.");
            }
            
            if (line == "\n" || line == "" || !line.compare(0,1,"#")) {
                continue;
            }
            
            const auto strBegin = line.find_first_not_of(" ");
            const auto strEnd = line.find_last_not_of(" ");
            
            if (strBegin != std::string::npos) {
                line = line.substr(strBegin,(strEnd - strBegin + 1)) ;
            }
            else {
                continue;
            }
            
            //SE reads
            if (format == "BED")
            {
                std::string chrom,start_ss,end_ss,score_ss,strand;
                int start, end;
                
                std::stringstream ss;
                
                ss << line;
                try{
                    std::getline(ss,chrom,'\t');
                    std::getline(ss,start_ss,'\t');
                    std::getline(ss,end_ss,'\t');
                    std::getline(ss,qname,'\t');
                    std::getline(ss,score_ss,'\t');
                    std::getline(ss,strand,'\t');
                    
                    start = std::stoi(start_ss);
                    end = std::stoi(end_ss);
                    
                }
                catch (const std::invalid_argument& ia) {
                    //std::cerr << "Invalid argument: " << ia.what() << '\n';
                    //std::cerr << "Error in reading BED file : " << inputFile << '\n';
                    error("Error in reading BED file : " + inputFile);
                    
                    std::exit(1);
                }
                read_t * aligned_read = new read_t();
                aligned_read->chrom = chrom;
                aligned_read->start = start;
                aligned_read->end = end;
                aligned_read->strand = strand;
                
                if (qname == prev_read_id || prev_read_id == "") {
                    alignments_per_read.push_back(aligned_read);
                    prev_read_id = qname;
                    continue;
                }
                else {
                    read_t *cur_read = NULL;
                    if (alignments_per_read.size() == 1 && flag != MULTIONLY){  //unique read
                        
                        cur_read = alignments_per_read[0];
                        
                        readID += 1;
                        read_pair_t<read_t> * read_pair = new read_pair_t<read_t>();
                        read_pair->first_reads.push_back(cur_read);
                        read_pair->readID = readID;
                        
                        if(global_context.thread_number >1){
                            
                            thread_context_t<read_t> * thread_context = global_context.thread_contexts[current_thread_id];
                            pthread_spin_lock(&thread_context->cur_reads_lock);
                            
                            //thread_context->cur_reads->push_back(std::pair<bam1_t *,bam1_t *>(cur_read,NULL));
                            thread_context->cur_reads->push_back(read_pair);
                            
                            pthread_spin_unlock(&thread_context->cur_reads_lock);
                            current_thread_id++;
                            if(current_thread_id >= global_context.thread_number) current_thread_id = 0;
                        }
                        else {
                            if (flag == MULTIONLY) {
                                process_aligned_fragment_multi<read_t>(&global_context,global_context.thread_contexts[0],read_pair);
                            }
                            else {
                            process_aligned_fragment<read_t>(&global_context,global_context.thread_contexts[0],read_pair);
                            }
                        }
                    }
                    if (alignments_per_read.size() > 1 && flag != UNIQONLY)
                    {
                        //multi reads
                        //cur_read = alignments_per_read[0];
                        readID += 1;
                        read_pair_t<read_t> * read_pair = new read_pair_t<read_t>();
                        
                        for(unsigned int i=0;i<alignments_per_read.size();i++){
                            (read_pair->first_reads).push_back(alignments_per_read[i]);
                            //bam_destroy1(alignments_per_read[i]);
                        }
                        read_pair->readID = readID;
                        
                        if(global_context.thread_number >1){
                            
                            thread_context_t<read_t> * thread_context = global_context.thread_contexts[current_thread_id];
                            pthread_spin_lock(&thread_context->cur_reads_lock);
                            
                            thread_context->cur_reads->push_back(read_pair);
                            
                            pthread_spin_unlock(&thread_context->cur_reads_lock);
                            current_thread_id++;
                            if(current_thread_id >= global_context.thread_number) current_thread_id = 0;
                        }
                        else {
                            if (flag == MULTIONLY) {
                                process_aligned_fragment_multi<read_t>(&global_context,global_context.thread_contexts[0],read_pair);
                            }
                            else {
                            process_aligned_fragment<read_t>(&global_context,global_context.thread_contexts[0],read_pair);
                            }
                            
                        }
                    }
                    alignments_per_read.clear();
                    alignments_per_read.push_back(aligned_read);
                    prev_read_id = qname;
                    
                }
            }
            else {//pair end read BEDPE
                std::string chrom1,chrom2,start1_ss,start2_ss,end1_ss,end2_ss,score_ss,strand1,strand2;
                int start1,start2, end1,end2;
                
                std::stringstream ss;
                
                ss << line;
                try{
                    std::getline(ss,chrom1,'\t');
                    std::getline(ss,start1_ss,'\t');
                    std::getline(ss,end1_ss,'\t');
                    std::getline(ss,chrom2,'\t');
                    std::getline(ss,start2_ss,'\t');
                    std::getline(ss,end2_ss,'\t');
                    std::getline(ss,qname,'\t');
                    std::getline(ss,score_ss,'\t');
                    std::getline(ss,strand1,'\t');
                    std::getline(ss,strand2,'\t');
                    
                    start1 = std::stoi(start1_ss);
                    end1 = std::stoi(end1_ss);
                    start2 = std::stoi(start2_ss);
                    end2 = std::stoi(end2_ss);
                    
                }
                catch (const std::invalid_argument& ia) {
                    //std::cerr << "Invalid argument: " << ia.what() << '\n';
                    std::cerr << "Error in reading BED file : " << inputFile << '\n';
                    std::exit(1);
                }
                
                size_t flag_pos = qname.find('/');
                //cur_read_id = aligned_read.qname;
                if (flag_pos != std::string::npos)
                {
                    qname = qname.substr(0,flag_pos);
                }
                else {
                    if (boost::algorithm::ends_with(qname,".1") || boost::algorithm::ends_with(qname,".2")) {
                        flag_pos = qname.size() -2;
                        qname = qname.substr(0,flag_pos);
                    }
                }
                
                
                read_t * read1 = new read_t();
                read_t * read2 = new read_t();
                
                read1->chrom = chrom1;
                read2->chrom = chrom2;
                read1->start = start1;
                read2->start = start2;
                read1->end = end1;
                read2->end = end2;
                read1->strand = strand1;
                read2->strand = strand2;
                
                if (qname == prev_read_id ){
                    multi_read1.push_back(read1);
                    multi_read2.push_back(read2);
                    continue;
                }
                else {
                    read_t * cur_read1 = NULL;
                    read_t * cur_read2 = NULL;
                    
                    //multi-reads
                    if ((multi_read1.size() >1 || multi_read2.size() > 1) && flag != UNIQONLY)
                    {
                        readID += 1;
                        read_pair_t<read_t> * read_pair = new read_pair_t<read_t>();
                        read_pair->readID = readID;
                        
                        for(unsigned int i=0;i< multi_read1.size();i++)
                        {
                            //bam_destroy1(multi_read1[i]);
                            read_pair->first_reads.push_back(multi_read1[i]);
                        }
                        for(unsigned int i=0;i< multi_read2.size();i++)
                        {
                            //bam_destroy1(multi_read2[i]);
                            read_pair->second_reads.push_back(multi_read2[i]);
                        }
                        
                        if(global_context.thread_number > 1){
                            thread_context_t<read_t> * thread_context = global_context.thread_contexts[current_thread_id];
                            
                            pthread_spin_lock(&thread_context->cur_reads_lock);
                            
                            thread_context->cur_reads->push_back(read_pair);
                            
                            pthread_spin_unlock(&thread_context->cur_reads_lock);
                            current_thread_id++;
                            
                            if(current_thread_id >= global_context.thread_number)
                                current_thread_id = 0;
                        }
                        else {
                            if (flag == MULTIONLY) {
                                process_aligned_fragment_multi<read_t>(&global_context,global_context.thread_contexts[0],read_pair);
                            }
                            else {
                                process_aligned_fragment<read_t>(&global_context,global_context.thread_contexts[0],read_pair);
                            }
                            
                        }
                    }
                    if (multi_read1.size() <= 1 && multi_read2.size() <= 1 && flag != MULTIONLY)
                    {
                        //uniq read
                        if (multi_read1.size() == 1) {cur_read1 = multi_read1[0];}
                        if (multi_read2.size() == 1){ cur_read2 = multi_read2[0];}
                        
                        if (cur_read1  || cur_read2 ){
                            readID += 1;
                            read_pair_t<read_t> * read_pair = new read_pair_t<read_t>();
                            if (cur_read1 != NULL) {
                                read_pair->first_reads.push_back(cur_read1);
                            }
                            if(cur_read2 != NULL) {
                                read_pair->second_reads.push_back(cur_read2);
                            }
                            read_pair->readID = readID;
                            
                            if(global_context.thread_number > 1){
                                
                                thread_context_t<read_t> * thread_context = global_context.thread_contexts[current_thread_id];
                                
                                pthread_spin_lock(&thread_context->cur_reads_lock);
                                
                                thread_context->cur_reads->push_back(read_pair);
                                
                                pthread_spin_unlock(&thread_context->cur_reads_lock);
                                current_thread_id++;
                                if(current_thread_id >= global_context.thread_number)
                                    current_thread_id = 0;
                            }
                            else {
                                if (flag == MULTIONLY) {
                                    process_aligned_fragment_multi<read_t>(&global_context,global_context.thread_contexts[0],read_pair);
                                }
                                else {
                                    process_aligned_fragment<read_t>(&global_context,global_context.thread_contexts[0],read_pair);
                                }
                            }
                        }
                    }
                    multi_read1.clear();
                    multi_read2.clear();
                    prev_read_id = qname;
                    
                    multi_read1.push_back(read1);
                    multi_read2.push_back(read2);
                }
                
            }
        }

        input.close();
    }
    catch(std::ifstream::failure e){
        //std::cout << "error in reading file " << inputFile << std::endl;
        error("Error in reading file : " + inputFile);
    }
    
    global_context.all_finished = 1;
    //join threads
    if(global_context.thread_number >1)
    {
        join_threads<read_t>(&global_context);
        
    }
    if (flag != MULTIONLY) {
        merge_results<read_t>(&global_context, track);

    }
    else {
        merge_results_redistr(&global_context, peak_reads_Prime,multi_read_mapTo_pids);
    }
    
    destroy_thread_context<read_t>(&global_context);
    
    if (flag != MULTIONLY) {
        track->isPE = global_context.isPE;
        track->set_rlengths(global_context.reflengths);
    }
    
}

/*bool sniff ( std::string format, std::string inputfile )
{
    //Detect whether the file format is the correct format for the input file. Rule: try to find the tag size using this parser, if error occurs or tag size is too small or too big, check is failed.

    if (format == "BAM") {
        std::string magic_header = fhd.read( 3 );
        if (magic_header == "BAM")
        {
            return true;
        }
        else {
            return false;
        }
    }
    if (format == "BED") {
        
        return true;
    }
    
    return false;
}

*/

std::pair<ShortRead *, ShortRead *> read_aligmentFile( opt_t &options , int flag)
{
    info("#1 read treatment tags...");
    ShortRead *treat = new ShortRead();
    ShortRead *control = new ShortRead();
    
    
    if (options.format == "BED") {
        parse_bed(treat,options.tfile,options.shift,options.threadNum,"BED",flag,nullptr,nullptr,nullptr);
        if (options.cfile != "") {
            if (flag == BOTH) {
                parse_bed(control, options.cfile,options.shift,options.threadNum,"BED",MULTICOUNTSONLY,nullptr,nullptr,nullptr);
            }
            else {
                parse_bed(control, options.cfile,options.shift,options.threadNum,"BED",flag,nullptr,nullptr,nullptr);
            }
        }
        //SE : options.tsize is tag size
        options.tsize = treat->tsize;
        options.PE_mode = false;
    }
    
    if (options.format == "BAM" || options.format == "SAM") {
        
        parse_bam(treat,options.tfile,options.shift,options.threadNum,flag,nullptr,nullptr,nullptr);
        parse_bam(control,options.cfile,options.shift,options.threadNum,flag,nullptr,nullptr,nullptr);
        
      /*ADD  if (options.cfile != "") {
            if (flag == BOTH) {
                parse_bam(control,options.cfile,options.shift,options.threadNum,MULTICOUNTSONLY,nullptr,nullptr,nullptr);
            }
            else {
                parse_bam(control,options.cfile,options.shift,options.threadNum,flag,nullptr,nullptr,nullptr);
            }
        }*/
        if (treat->isPE) {
            options.PE_mode = true;
            options.tsize = treat->fraglength;
        }
        else {
            options.tsize = treat->tsize;
        }
    }
    if (options.format == "BEDPE")
    {
        parse_bed(treat,options.tfile,options.shift,options.threadNum,"BEDPE",flag,nullptr,nullptr,nullptr);
        if (options.cfile != "") {
            if (flag == BOTH) {
                parse_bed(control,options.cfile,options.shift,options.threadNum, "BEDPE",MULTICOUNTSONLY,nullptr,nullptr,nullptr);
            }
            else {
                parse_bed(control,options.cfile,options.shift,options.threadNum, "BEDPE",flag,nullptr,nullptr,nullptr);
            }
        //PE : options.tsize is average fragment length
        options.tsize = treat->fraglength;
        options.PE_mode = true;
        }
    }
    
    return std::pair<ShortRead *,ShortRead *>(treat, control);
}

//multi-read re-distribution among candidate peaks.
void read_distribution(ShortRead * track, std::string inputFile,Candidate_Peaks * peakIdx,opt_t options,std::vector<double> * peak_reads_Prime, MULTI_READ_MAP * multi_read_mapTo_pids)
{
    
    if (options.format == "BAM") {
        int shiftsize = options.d/2;
        parse_bam(track,inputFile,shiftsize,options.threadNum,MULTIONLY,peakIdx,peak_reads_Prime,multi_read_mapTo_pids);
    }
    if (options.format == "BED") {
        int shiftsize = options.d/2;
        parse_bed(track,inputFile,shiftsize,options.threadNum, "BED",MULTIONLY,peakIdx,peak_reads_Prime,multi_read_mapTo_pids);
    }
    
    
}


/*
 int main() {
 
 int numItr = 0;
 // std::string peakbed = "test_peak.bed0";
 char out_dir[100] = "./";
 char tfile[100] = "test_ip.bam";
 char input[100] = "test_input.bam";
 char peakbed[100] = "test_peak.bed0";
 double sf_t = 1.0;
 double sf_c = 1.0;
 int shiftsize = 50;
 bool wig = false;
 char prj_name[10] = "test";
 double pvalCutoff = 1e-5;
 int thread_num = 1;
 char format[10] = "BAM";
 
 Candidate_Peaks * peakIdx = new Candidate_Peaks(peakbed);
 
 std::vector<double> peak_reads = {12,16,10};
 std::vector<double> peak_reads_Prime = {0,0,0};
 std::vector<double> Input_peakReads = {5,4,3};
 
 
 std::vector<double> multiAlgn_weight;
 std::vector<double> multiAlgn_final_weight;
 std::vector<int> multiAlgn_To_multiRead ;
 std::vector<int> multiAlgn_To_peakdID;
 
 std::map<int,std::vector<int> > multi_read_mapTo_pids ;
 
 multi_read_mapTo_pids.insert(std::pair<int,std::vector<int> > (0,{0,1}));
 multi_read_mapTo_pids.insert(std::pair<int,std::vector<int> > (1,{1,2}));
 
 prepare_EM_structure(peakIdx,peak_reads, multi_read_mapTo_pids, multiAlgn_To_multiRead, multiAlgn_weight, multiAlgn_To_peakdID);
 
 EqualWeight_reads(multiAlgn_To_multiRead, peak_reads, peak_reads_Prime, multiAlgn_To_peakdID);
 
 for (size_t i = 0; i < multiAlgn_To_multiRead.size(); i++) {
 std::cout << "multiAlgn to multireads " << multiAlgn_To_multiRead[i] << std::endl;
 }
 for (size_t i = 0; i < multiAlgn_To_peakdID.size(); i++) {
 std::cout << "multiAlgn_To_peakdID " << multiAlgn_To_peakdID[i] << std::endl;
 }
 for (size_t i = 0; i < multiAlgn_weight.size(); i++) {
 std::cout << "multiAlgn_weight " << multiAlgn_weight[i] << std::endl;
 }
 
 for (size_t i =0; i < peak_reads_Prime.size(); i++) {
 std::cout << peakIdx->get_chrom(i) << "\t" << peakIdx->get_start(i) << "\t" << peakIdx->get_end(i) << "\t" << peak_reads_Prime[i] << std::endl;
 }
 
 //filter_peaks(peakIdx, IP_peakReads, Input_peakReads, 1.0, 1.0,1e-5,"./","test");
 
 //delete peakIdx;
 
 
 //std::vector<int> effLengths;
 
 run_EM_TEpeaks_BAM(out_dir,tfile,input,peakbed, sf_t, sf_c, shiftsize, wig, prj_name, pvalCutoff, thread_num);
 
 
 
 }*/

