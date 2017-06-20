//
//  EM_TEpeaks_BAM.cpp
//  TEToolkit_c++
//
//  Created by Ying Jin on 2/22/16.
//  Copyright (c) 2016 Ying Jin. All rights reserved.
//

#include "EM_TEpeaks.h"

//#include <malloc.h>
#include <string>
#include <stdlib.h>
#include <pthread.h>
#include <math.h>
#include <algorithm>
#include <iostream>
#include <fstream>
#include <unistd.h>
#include <ctime>
#include <sstream>

#include "htslib/sam.h"
//#include "IntervalTree.h"
#include "Candidate_Peaks.h"

#include "EMestimate_reads.h"
#include "zeroin.h"
#include "myLog.h"
//#include "Parser.h"



void prepare_EM_structure(Candidate_Peaks *peakIdx,std::vector<double> peak_reads, MULTI_READ_MAP * multi_read_mapTo_pids,std::vector<int> & multiAlgn_To_multiRead,std::vector<double> & multiAlgn_weight, std::vector<int> & multiAlgn_To_peakdID)
{
    std::map<int,std::vector<int> >::iterator itr;
    //multi-aligment -> peakID
    int multiAlgn_id = 0;
    
    //std::cout << "multi read map to peak " << multi_read_mapTo_pids.size() << std::endl;
    for (itr = multi_read_mapTo_pids->begin(); itr != multi_read_mapTo_pids->end(); itr ++) {
        //int readID = itr->first;
        //std::cout << "multi read " << readID << std::endl; 
        std::vector<int> pidlist = itr->second;
        
        double total = 0.0; // std::accumulate(counts_.begin(), counts_.end(), 0.0);;
        for (size_t i=0; i < pidlist.size(); i++) {
            //std::cout << pidlist[i] << std::endl;
            total += 1.0 * peak_reads[pidlist[i]]/peakIdx->get_length(pidlist[i]);
            multiAlgn_id += 1;
            multiAlgn_To_peakdID.push_back(pidlist[i]);
        }
        multiAlgn_To_multiRead.push_back(multiAlgn_id);
        for (size_t i=0; i < pidlist.size(); i++) {
            double w = 1.0 *peak_reads[pidlist[i]]/(peakIdx->get_length(pidlist[i])*total);
            multiAlgn_weight.push_back(w);
        }
        
    }
  //  if (multiAlgn_id > 0) {
  //      multiAlgn_To_multiRead.push_back(multiAlgn_id);
   // }
    
}

//equaly weight multi-aignments
//the input peak_reads contains unique reads only
void EqualWeight_reads(std::vector<int> multiAlgn_To_multiRead, std::vector<double> peak_reads, std::vector<double> & peak_reads_Prime,std::vector<int> multiAlgn_To_peakID)
{
    //debug("in EqualWeight_reads");
    int multiAlgn_first_idx = 0;
    //std::vector<double> tmp_peak_reads(peak_reads.size(),0.0);
    for (size_t i= 0; i < peak_reads.size(); i++) {
        peak_reads_Prime[i] = peak_reads[i];
    }
    
    
    for (size_t i=0; i< multiAlgn_To_multiRead.size(); i++) {
        //loop through all multi-alignments of each multi-read
        double w = 1.0 / (multiAlgn_To_multiRead[i] - multiAlgn_first_idx);
        
        for (int j= multiAlgn_first_idx; j < multiAlgn_To_multiRead[i]; j++) {
            int pid = multiAlgn_To_peakID[j];
            
            peak_reads_Prime[pid] += w;
        }
        multiAlgn_first_idx = multiAlgn_To_multiRead[i];
    }
}

//multi-read re-distribution among candidate peaks.
void EM_assignment(ShortRead * track, std::string inputFile,Candidate_Peaks * peakIdx,opt_t options,std::vector<double> & peak_reads_Prime)
{
    
    std::vector<double> peak_reads(peak_reads_Prime.size(),0.0);
    MULTI_READ_MAP * multi_read_mapTo_pids = new MULTI_READ_MAP();
    
    
    read_distribution(track,inputFile,peakIdx,options,&peak_reads, multi_read_mapTo_pids);


    
    //EM re-distribute multi-reads
    std::vector<double> multiAlgn_weight;
    std::vector<double> multiAlgn_final_weight;
    std::vector<int> multiAlgn_To_multiRead;
    std::vector<int> multiAlgn_To_peakID;
    
    
    /*for(size_t i = 0 ; i < peak_reads.size() ; i++)
    {
       std::cout << "peak " << i << "\t" << peak_reads[i] << std::endl;
    }*/
    prepare_EM_structure(peakIdx,peak_reads, multi_read_mapTo_pids, multiAlgn_To_multiRead, multiAlgn_weight, multiAlgn_To_peakID);
    
    /*for(size_t i = 0 ; i < peak_reads.size() ; i++)
    {
       std::cout << "peak " << i << "\t" << peak_reads[i] << std::endl;
    }
    for(size_t i = 0 ; i < multiAlgn_To_multiRead.size() ; i++)
    {
      std::cout << "multiAlignment to multireads " << multiAlgn_To_multiRead[i] << std::endl;
    }*/
    if (options.numItr == 0 ) { // eqauly distribute multi-alignment
        EqualWeight_reads(multiAlgn_To_multiRead, peak_reads, peak_reads_Prime, multiAlgn_To_peakID);
    }
    else {
        EMestimate_read(multiAlgn_To_multiRead,multiAlgn_To_peakID,multiAlgn_weight,peak_reads,peak_reads_Prime,options.numItr,peakIdx->peak_length);
    }
    /*for(size_t i = 0 ; i < peak_reads_Prime.size() ; i++)
    {
        std::cout << "peak reads prime " << peak_reads_Prime[i] << std::endl;
    }*/
    
}


void filter_peaks(std::vector<int> peaks_with_multiReads,Candidate_Peaks * peakIdx, std::vector<double> IP_peakReads, std::vector<double> Input_peakReads, opt_t options)
{
    
    try {
        //int fragsize = options.fragsize;
        std::string project_name = options.project_name;
       // double pvalCutoff = std::pow(10,options.log_pvalue);
        std::string data_outdir = options.data_outdir;
        
        //ADD int window_size = fragsize * 2;
        std::ofstream ROUT;
        
        
        ROUT.open (data_outdir +"/" +project_name+ "_multiRead_peaks.txt", std::ofstream::out);
        
        ROUT << "chrom\t" << "start\t" << "end\t" << "IP_tags\t" << "Input_tags\t" << "pileup\t" << "fold enrichment\t" << "p-value\n";
        
        for (size_t i =0; i < IP_peakReads.size(); i++) {
            
            if (peaks_with_multiReads[i] == 0) {
                continue;
            }
            double ip_tag = IP_peakReads[i];
            double input_tag = Input_peakReads[i];
            
            double fe =  input_tag >0 ? ip_tag * options.treat_scale/(input_tag * options.ctrl_scale) : (ip_tag + 0.1)* options.treat_scale/(0.1 * options.ctrl_scale);
            //normalize to 200bp window
            //ADD int peak_len = peakIdx->get_length(i);
            double pileup = 1.0 * ip_tag * options.treat_scale;
            
            //ADD int cnt_per_window = int(peakIdx->get_count(i) * window_size /peak_len);
            
            //double p = std::pow(10,peakIdx->get_pval(i) * (-10));
            
            //ADD double lam_bg = retrieve_lambda_bg(peakIdx->get_count(i),peakIdx->get_pval(i));
            
            //ADD double lambda = std::max(lam_bg,input_tag * options.ctrl_scale);
            
            //std::cout << "lam_bg = " << lam_bg << "\tlambda = " << lambda << std::endl;
            
            //ADD int cur_cnt_per_window = int(ip_tag * options.treat_scale * window_size/peak_len) ;
            
            //double pval = 1.0 - ppois(cur_cnt_per_window,lambda);
            double pval = peakIdx->get_pval(i);
            
            //inv possion to get local lambda
            if (pileup > options.pileup && fe > options.fe ) {
                ROUT << peakIdx->get_chrom(i) << "\t" << peakIdx->get_start(i) << "\t" << peakIdx->get_end(i) << "\t" << ip_tag << "\t" << input_tag << "\t" << pileup <<
                "\t" << fe << "\t" << pval << "\n";
            }
            
        }
        
        ROUT.close();
    }catch(std::ofstream::failure e ){
        std::cout << "Error in writing filtered peaks.\n" << std::endl;
    }

}



int run_EM_TEpeaks(opt_t options, ShortRead * treat, ShortRead * control,std::string peak_fname)
{

    std::vector<double> IP_peakReads;
    //std::vector<double> * IP_peakReads_ptr;
    std::vector<double> Input_peakReads;
    //std::vector<double> * Input_peakReads_ptr;
    std::vector<double> peak_reads_Prime;
    std::vector<double> ctrl_peak_reads_Prime;
    std::vector<int> peaks_with_multiReads;
    
    MULTI_READ_MAP *  multi_read_mapTo_pids = new MULTI_READ_MAP();
    //std::map<int, std::vector<int>> * multi_read_mapTo_pids_ptr;
    
    info("Read in candidate peaks ...");
    Candidate_Peaks * peakIdx = new Candidate_Peaks(peak_fname);

    
    info("Done reding candidate peaks.");
    //debug(" in run_EM_TEpeaks numofpeaks " + std::to_string(peakIdx->get_numofpeaks()));
    
    if (peakIdx == NULL) {
        std::cout << "error in read candidate peaks!" << std::endl;
        std::exit(1);
    }
    if (peakIdx->get_numofpeaks() < 10) {
        std::cout << "less than 10 candidate peak regions! Stop EM algorithm" << std::endl;
        std::exit(1);
    }
    //Candidate_Peaks peakIdx(peakbed);
    
    for(int i =0; i< peakIdx->get_numofpeaks() ; i++ )
    {
        IP_peakReads.push_back(0.0);
        Input_peakReads.push_back(0.0);
        peak_reads_Prime.push_back(0.0);
        ctrl_peak_reads_Prime.push_back(0.0);
        peaks_with_multiReads.push_back(0);
    }
    
    //parse IP /input files for multi-reads
    //IP_peakReads_ptr = &IP_peakReads;
    //multi_read_mapTo_pids_ptr = &multi_read_mapTo_pids;
    
    read_distribution(treat, options.tfile,peakIdx,options,&IP_peakReads, multi_read_mapTo_pids);
    
    if (multi_read_mapTo_pids->size()==0) {
        
        delete peakIdx;
        
        delete multi_read_mapTo_pids;
        

        
        info("no multi_reads mapped to candidate peak regions.");
        return 0;
    }
    //EM re-distribute multi-reads
    std::vector<double> multiAlgn_weight;
    std::vector<double> multiAlgn_final_weight;
    std::vector<int> multiAlgn_To_multiRead;
    std::vector<int> multiAlgn_To_peakID;
    
    //debug(" multi_read_mapTo_pids size " + std::to_string(multi_read_mapTo_pids->size()));
    
    prepare_EM_structure(peakIdx, IP_peakReads, multi_read_mapTo_pids, multiAlgn_To_multiRead, multiAlgn_weight, multiAlgn_To_peakID);
    
    //debug("after prepare EM structure numItr = " + std::to_string(options.numItr) );
    
    if (options.numItr == 0 ) { // eqauly distribute multi-alignment
        
        EqualWeight_reads(multiAlgn_To_multiRead, IP_peakReads, peak_reads_Prime, multiAlgn_To_peakID);
    }
    else {
        EMestimate_read(multiAlgn_To_multiRead,multiAlgn_To_peakID,multiAlgn_weight,IP_peakReads,peak_reads_Prime,options.numItr,peakIdx->peak_length);
        
        //debug("after EMestimate_read");
    }

    int multiAlgn_first_idx = 0;
    
    for (size_t i=0; i< multiAlgn_To_multiRead.size(); i++) {
        //loop through all multi-alignments of each multi-read
        
        for (int j= multiAlgn_first_idx; j < multiAlgn_To_multiRead[i]; j++) {
            int pid = multiAlgn_To_peakID[j];
            
            peaks_with_multiReads[pid] = 1;
        }
        multiAlgn_first_idx = multiAlgn_To_multiRead[i];
    }
    
    //info("Reding Input alignment file ..." ;
    //parse input file
    multi_read_mapTo_pids->clear();
    multiAlgn_weight.clear();
    multiAlgn_final_weight.clear();
    multiAlgn_To_multiRead.clear();
    multiAlgn_To_peakID.clear();
    //peak_reads_Prime.clear();
    
    read_distribution(control, options.cfile,peakIdx,options,&Input_peakReads,multi_read_mapTo_pids);
    
    //debug("after control read distribution");
    
    prepare_EM_structure(peakIdx, Input_peakReads, multi_read_mapTo_pids, multiAlgn_To_multiRead, multiAlgn_weight, multiAlgn_To_peakID);
    
    
    EqualWeight_reads(multiAlgn_To_multiRead, Input_peakReads, ctrl_peak_reads_Prime, multiAlgn_To_peakID);
    
    
    //filter peaks
    
    filter_peaks(peaks_with_multiReads,peakIdx, peak_reads_Prime, ctrl_peak_reads_Prime, options);
    
    delete peakIdx;
    delete multi_read_mapTo_pids;
    

    //write_filtered_peaks(peakIdx,data_outdir,project_name);
   
   // std::cout << "INFO  @ " << cur_time << "\tDone." << std::endl;
    
    return 1;
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

