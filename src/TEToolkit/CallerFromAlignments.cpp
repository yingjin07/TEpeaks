//
//  CallerFromAlignments.cpp
//  TEToolkit_c++
//
//  Created by Ying Jin on 5/31/16.
//  Copyright (c) 2016 Ying Jin. All rights reserved.
//

//Modified from MACS2

#include <stdio.h>
#include <stdlib.h>
#include <iostream>
//#include <fstream>
#include <sstream>
#include <cstdio>
#include <iomanip>
#include <map>
#include <vector>
#include <functional>
#include <math.h>
#include <boost/lexical_cast.hpp>
#include "CallerFromAlignments.h"
#include "myLog.h"
#include "zeroin.h"
#include "cStatistics.h"
#include <sys/time.h>

using boost::lexical_cast;


void apply_multiple_cutoffs ( std::vector<std::vector<double>> multiple_score_arrays, std::vector<double> multiple_cutoffs,std::vector<int> &indices)
{
    int array_length = multiple_score_arrays[0].size();
    std::vector<int> tmp;
    
    for (int i=0; i < array_length; i++) {
        tmp.push_back(0);
    }
    
    for (int i= 0; i < array_length;i++)
    {
        
        for (size_t j=0;j < multiple_cutoffs.size(); j++){
            if (multiple_score_arrays[j][i] > multiple_cutoffs[j]) {
                tmp[i] = 1;
                break;
            }
        }
        
    }
    
    for (int i =0 ; i < array_length; i++) {
        if (tmp[i] == 1) {
            indices.push_back(i);
        }
    }
}

double CallerFromAlignments::get_pscore ( int observed, double expectation )
{
    /*"""Get p-value score from Poisson test. First check existing
     table, if failed, call poisson_cdf function, then store the result
     in table.
     
     """*/
    
    double score;
    //long key_value;
    double expectation_round = round(expectation * 100000)/100000;
    std::string key_str = lexical_cast<std::string>(observed) + ":" + lexical_cast<std::string>(expectation_round);
    
    //std::hash<std::string> str_hash;
    //size_t key_value =  str_hash( key_str );
    
    if (pscore_map.find(key_str) != pscore_map.end()) {
        //debug("********found in pscore map");
        return pscore_map[key_str];
    }
    else {
        score = -1.0 * log10_poisson_cdf(observed,expectation_round,0);
        pscore_map.insert(std::pair<std::string,double>(key_str, score));
    }
    
    return score;
}

CallerFromAlignments::CallerFromAlignments(ShortRead * treat, ShortRead * ctrl,std::vector<int> ctrl_d_s,
                                           std::vector<double> ctrl_scaling_factor_s,
                                           int d ,
                                           double treat_scaling_factor ,
                                           double pseudocount ,
                                           int end_shift,
                                           double lambda_bg , bool no_lambda_flag,
                                           bool save_bedGraph,
                                           std::string  bedGraph_filename_prefix ,
                                           std::string bedGraph_treat_filename ,
                                           std::string bedGraph_control_filename ,
                                           std::string cutoff_analysis_filename ,
                                           bool save_SPMR )
{
        /*"""Initialize.

        A calculator is unique to each comparison of treat and
        control. Treat_depth and ctrl_depth should not be changed
        during calculation.

        treat and ctrl are two ShortRead objects.

        treat_depth and ctrl_depth are effective depth in million:
                                    sequencing depth in million after
                                    duplicates being filtered. If
                                    treatment is scaled down to
                                    control sample size, then this
                                    should be control sample size in
                                    million. And vice versa.

        d, sregion, lregion: d is the fragment size, sregion is the
                             small region size, lregion is the large
                             region size
                                    
        pseudocount: a pseudocount used to calculate logLR, FE or
                     logFE. Please note this value will not be changed
                     with normalization method. So if you really want
                     to set pseudocount 1 per million reads, set it
                     after you normalize treat and control by million
                     reads by `change_normalizetion_method(ord('M'))`.

        """ */
    if (ctrl_d_s.size() == 0) {
        ctrl_d_s.push_back(200);
        ctrl_d_s.push_back(1000);
        ctrl_d_s.push_back(10000);
    }
    if (ctrl_scaling_factor_s.size() == 0) {
        ctrl_scaling_factor_s.push_back(1.0);
        ctrl_scaling_factor_s.push_back(0.2);
        ctrl_scaling_factor_s.push_back(0.02);
    }
    
    
    if (treat->isPE) {
        PE_mode = true;
    }
    else {
        PE_mode = false;
    }
    this->treat = treat;
    this->ctrl = ctrl;
    
    if (this->ctrl == NULL) {
        this->ctrl = this->treat;
    }
        
    this->trackline = false;
    this->d = d;
    this->ctrl_d_s = ctrl_d_s;

    this->treat_scaling_factor = treat_scaling_factor;
    this->ctrl_scaling_factor_s= ctrl_scaling_factor_s;
    this->end_shift = end_shift;
    this->lambda_bg = lambda_bg;
    //this->pqtable = NULL;
    
    this->save_bedGraph = save_bedGraph;
    this->save_SPMR = save_SPMR;
    this->bedGraph_filename_prefix =  bedGraph_filename_prefix;
    //tmp_bytes = bedGraph_treat_filename.encode('UTF-8')
    //print bedGraph_treat_filename, tmp_bytes
    this->bedGraph_treat_filename = bedGraph_treat_filename;
    //tmp_bytes = bedGraph_control_filename.encode('UTF-8')
    //print bedGraph_control_filename, tmp_bytes
    this->bedGraph_control_filename = bedGraph_control_filename;

    this->no_lambda_flag = no_lambda_flag;

    this->pseudocount = pseudocount;

    std::vector<std::string> chr1 = this->treat->get_chrom_names();
    std::vector<std::string> chr2 = this->ctrl->get_chrom_names();
    
    for (auto ch : chr1) {
        //if (chr2.find(ch) != chr2.end())
        if (std::find(chr2.begin(), chr2.end(), ch) != chr2.end() )
        {
            this->chromosomes.push_back(ch);
        }
    }

    this->test_time = 0;

    double f = 0.3;
    // step for optimal cutoff is 0.3 in -log10pvalue, we try from pvalue 1E-10 (-10logp=10) to 0.5 (-10logp=0.3)
    while(f < 10) {
        this->pvalue_length.insert(std::pair<double,int>(f,0));
        this->pvalue_npeaks.insert(std::pair<double,int>(f,0));
        f += 0.3;
    }
    this->optimal_p_cutoff = 0;
    this->cutoff_analysis_filename = cutoff_analysis_filename;
}

CallerFromAlignments::~CallerFromAlignments(){
        /*"""Remove temparary files for pileup values of each chromosome.

        Note: This function MUST be called if the class object won't
        be used anymore.

        """*/
    
    //std::string f;

    for(auto f : pileup_data_files)
    {
        std::ifstream ff(f.second);
        if( ff.good() ){
            std::remove(f.second.c_str());
        }
    }
    return ;
}

void CallerFromAlignments::__chrom_pair_treat_ctrl ( std::pair<std::vector<int>,std::vector<double> > treat_pv, std::pair<std::vector<int>,std::vector<double> > ctrl_pv )
{
        /*"""*private* Pair treat and ctrl pileup for each region.
        
        treat_pv and ctrl_pv are [np.ndarray, np.ndarray].

        return [p, t, c] list, each element is a vector.
        """*/

    int lt = treat_pv.first.size();
    int lc = ctrl_pv.first.size();

    //std::vector<pos_tc_t> ret;

    //int pre_p = 0;
    //int index_ret = 0;
    int it = 0;
    int ic = 0;

    while (it < lt and ic < lc){
        if (treat_pv.first[it] < ctrl_pv.first[ic])
        {
            // clip a region from pre_p to p1, then set pre_p as p1.
            chr_pos_treat_ctrl.pos.push_back(treat_pv.first[it]);
            chr_pos_treat_ctrl.treat_v.push_back(treat_pv.second[it]);
            chr_pos_treat_ctrl.ctrl_v.push_back(ctrl_pv.second[ic]);

            //pre_p = treat_pv.first[it];
            //index_ret += 1;
                // call for the next p1 and v1
            it += 1;
        }
        else {
        if(treat_pv.first[it] > ctrl_pv.first[ic]){
                // clip a region from pre_p to p2, then set pre_p as p2.
            
            chr_pos_treat_ctrl.pos.push_back(ctrl_pv.first[ic]);
            chr_pos_treat_ctrl.treat_v.push_back(treat_pv.second[it]);
            chr_pos_treat_ctrl.ctrl_v.push_back(ctrl_pv.second[ic]);

                //pre_p = ctrl_pv.first[ic];
               // index_ret += 1;
                // call for the next p2 and v2
                ic += 1;
        }
        else{
                // from pre_p to p1 or p2, then set pre_p as p1 or p2.
            chr_pos_treat_ctrl.pos.push_back(treat_pv.first[it]);
            chr_pos_treat_ctrl.treat_v.push_back(treat_pv.second[it]);
            chr_pos_treat_ctrl.ctrl_v.push_back(ctrl_pv.second[ic]);

           // pre_p = treat_pv.first[it];
            //index_ret += 1;
                // call for the next p1, v1, p2, v2.
            it += 1;
            ic += 1;

        }
        }
    }
    
    //return ret;

}


void CallerFromAlignments::__pileup_treat_ctrl_a_chromosome ( std::string chrom , bool uniqOnly )
{
        /*"""After this function is called, chr_pos_treat_ctrl will
        be reset and assigned to the pileup values of the given
        chromosome.
        
        """*/
    
    std::pair<std::vector<int>, std::vector<double> > treat_pv, ctrl_pv;
    //long i;
    //float t;
    std::ifstream f;

    if (std::find(chromosomes.begin(),chromosomes.end(),chrom) == chromosomes.end()) {
        error("chromosome is not valid : " + chrom);
        std::exit(1);
    }
    //assert chrom in self.chromosomes, "chromosome %s is not valid." % chrom;

    // check backup file of pileup values. If not exists, create
    // it. Otherwise, load them instead of calculating new pileup values.
    if (uniqOnly) {
    
    if (pileup_data_files.find( chrom ) != pileup_data_files.end())
    {
        
     try{
         f.open( pileup_data_files[ chrom ],std::ifstream::in );
         chr_pos_treat_ctrl.pos.clear();
         chr_pos_treat_ctrl.treat_v.clear();
         chr_pos_treat_ctrl.ctrl_v.clear();
         
        while(! f.eof()){
            std::string line, ss_pos,ss_treat,ss_ctrl;
            std::stringstream ss;
            
            if (! std::getline(f,line)){
                break;
            }
            
            ss << line;
            std::getline(ss,ss_pos,'\t');
            std::getline(ss,ss_treat,'\t');
            std::getline(ss,ss_ctrl,'\t');
            
            chr_pos_treat_ctrl.pos.push_back(std::stol(ss_pos));
            chr_pos_treat_ctrl.treat_v.push_back(std::stod(ss_treat));
            chr_pos_treat_ctrl.ctrl_v.push_back(std::stod(ss_ctrl));
            
         }
            
         f.close();
            return;
        }
        catch(std::ifstream::failure e ){
            struct timeval start;
            gettimeofday(&start, NULL);
            std::string temp_filename = lexical_cast<std::string>(start.tv_sec) + "." + lexical_cast<std::string> ( start.tv_usec);
            //temp_fd, temp_filename = mkstemp();
            //os.close(temp_fd);
            pileup_data_files[ chrom ] = temp_filename;
        }
    }
    else{
        //temp_fd, temp_filename = mkstemp();
        //os.close(temp_fd);
        //pileup_data_files[ chrom ] = temp_filename;
        struct timeval start;
        gettimeofday(&start, NULL);
        std::string temp_filename = lexical_cast<std::string>(start.tv_sec) + "." + lexical_cast<std::string> ( start.tv_usec);
        
        pileup_data_files[ chrom ] = temp_filename;
    }
    }


    // reset or clean existing chr_pos_treat_ctrl
    chr_pos_treat_ctrl.pos.clear();
    chr_pos_treat_ctrl.treat_v.clear();
    chr_pos_treat_ctrl.ctrl_v.clear();

    std::vector<double> scale_factors;
    std::vector<int> ds;
    ds.push_back(this->d);
    scale_factors.push_back(treat_scaling_factor);
    
    if (PE_mode){
        
            treat_pv = treat->pileup_a_chromosome_pe ( chrom, scale_factors, 0.0 );
    }
    else{
        
            treat_pv = treat->pileup_a_chromosome( chrom, ds, scale_factors,uniqOnly, 0.0,true, this->end_shift );
        
        
    }
    
    
        
    if (not no_lambda_flag){
            if (PE_mode){
                // note, we pileup up PE control as SE control because
                // we assume the bias only can be captured at the
                // surrounding regions of cutting sites from control experiments.
                ctrl_pv = ctrl->pileup_a_chromosome_c( chrom, this->ctrl_d_s, this->ctrl_scaling_factor_s, this->lambda_bg );
            }
            else{
                
                
                ctrl_pv = ctrl->pileup_a_chromosome( chrom, this->ctrl_d_s, this->ctrl_scaling_factor_s,uniqOnly,
                                                         this->lambda_bg,
                                                        false );
                
            }
    }
    else{
        std::vector<int> pos;
        pos.push_back(treat_pv.first[treat_pv.first.size()-1]);
        std::vector<double> val;
        val.push_back(this->lambda_bg);
        
        ctrl_pv = std::pair<std::vector<int>,std::vector<double> > (pos,val);
         //   ctrl_pv = [treat_pv.first[][0][-1:], np.array([self.lambda_bg,], dtype="float32")]; // set a global lambda
    }
    

    __chrom_pair_treat_ctrl( treat_pv, ctrl_pv);
    
    // clean treat_pv and ctrl_pv
    treat_pv.first.clear(); // = [];
    treat_pv.second.clear();
    ctrl_pv.first.clear(); //  = [];
    ctrl_pv.second.clear();
    
    // save data to temporary file
    if (uniqOnly) {
    
    try{
        std::ofstream f(this->pileup_data_files[ chrom ],std::ofstream::out);
        //    cPickle.dump( chr_pos_treat_ctrl, f , protocol=2 );
        
        for (size_t i=0 ; i < chr_pos_treat_ctrl.pos.size(); i ++ ) {
            
            double tr_v = round(chr_pos_treat_ctrl.treat_v[i] * 100000)/100000;
            double ct_v = round(chr_pos_treat_ctrl.ctrl_v[i] * 100000)/100000;
            
            
            f  << chr_pos_treat_ctrl.pos[i] << "\t" << tr_v << "\t" << ct_v << "\n";
        }
        f.close();
    }
    catch(std::ofstream::failure e){
        // fail to write then remove the key in pileup_data_files
        this->pileup_data_files.erase(chrom);
    }
    }
    
    return ;
}



std::vector<double> CallerFromAlignments::__cal_pscore ( std::vector<double> array1, std::vector<double> array2 )
{
    std::vector<double> s;
    
    if (array1.size() != array2.size()) {
        error("Two arraies shouble have the same length.");
        std::exit(1);
    }
    
    for (size_t i=0 ; i< array1.size();i++)
    {
        s.push_back( get_pscore( int(array1[i]), array2[i] ));
    }
    return s;
}



void CallerFromAlignments::__cal_pvalue_qvalue_table ()
{
       /* """After this function is called, pqtable is built. All
        chromosomes will be iterated. So it will take some time.
        
        """*/

    std::map<double,int> pvalue_stat ;

    info( "#4 Start to calculate pvalue stat..." );

    long N = 0;
    std::vector<double> keys;
    
    for(size_t i =0; i< chromosomes.size(); i++ )
    {
        std::string chrom = chromosomes[ i ];
        int pre_p = 0;

        //debug("In __cal_pvalue_qvalue_table i = "+ std::to_string(i));
        __pileup_treat_ctrl_a_chromosome( chrom , true); //unique reads only
        //debug("after call pileup_treat_ctrl_a_chromosome");
        
       

        for (size_t j=0; j< chr_pos_treat_ctrl.pos.size() ; j++)
        {
            double this_v = get_pscore( int(chr_pos_treat_ctrl.treat_v[j]), chr_pos_treat_ctrl.ctrl_v[j] );
            
            int this_l = chr_pos_treat_ctrl.pos[j] - pre_p;
            

            N += this_l;
            
            if (pvalue_stat.find( this_v ) != pvalue_stat.end()){
                pvalue_stat[ this_v ] += this_l;
            }
            else{
                pvalue_stat.insert(std::pair<double,int>(this_v,this_l));
                keys.push_back(this_v);
            }
            pre_p = chr_pos_treat_ctrl.pos[j];
        }
        //nhcal += chr_pos_treat_ctrl.pos.size();
        
    }
    
    //logging.debug ( "calculate pvalue/access hash for %d times" % nhcal )
    //logging.debug ( "access hash for %d times" % nhcal )
    //nhval = 0;

    //N = sum(pvalue_stat.values());// # total length
    int k = 1;//                           # rank
    double f = -1.0 * log10(N);
    //double pre_v = -2147483647;
    //int pre_l = 0;
    double pre_q = 2147483647; //             # save the previous q-value

    // pqtable = new Float64HashTable();
    
    //unique_values = sorted(pvalue_stat.keys(), reverse=True) #sorted(unique_values,reverse=True);
    std::sort(keys.begin(), keys.end(), std::greater<double>());
    
    
    for (auto key : keys)
    {
        double v = key;

        //double rounded_up = (ceil)(v * 100000) / 100000;
        int rounded_up = (int)(ceil)(v * 100000);
        //std::string key_str = lexical_cast<std::string>(v) ;
        
        //std::hash<std::string> str_hash;
        //size_t key_value =  str_hash( key_str );
        
        int l = pvalue_stat[key];
        
        double q = v + (log10(k) + f);
     //   debug("N = " + std::to_string(N) + " k " + std::to_string(k) + " f = " + std::to_string(f) +" pval " + std::to_string(v) + " qval " + std::to_string(q));
        q = std::max(0.0,std::min(pre_q,q)); //           # make q-score monotonic
        if (pqtable.find(rounded_up) == pqtable.end()) {
            pqtable.insert(std::pair<int,double>(rounded_up,q));
        }
        else {
            pqtable[rounded_up] = q;
        }
        pre_q = q;
        k +=l;
        //nhcal += 1;
    }
    
    try{
        std::ofstream f("test_pqtable.txt",std::ofstream::out);
        //    cPickle.dump( chr_pos_treat_ctrl, f , protocol=2 );
        for(auto k : pqtable){
            
            f << k.first << "\t" << k.second << "\n";
            
        }
        f.close();
    }
    catch(std::ofstream::failure e){
        // fail to write then remove the key in pileup_data_files
        this->pileup_data_files.erase(chrom);
    }
}

void CallerFromAlignments::call_peaks ( PeakIO *peaks, std::vector<std::string> scoring_function_symbols, std::vector<double> score_cutoff_s, bool uniqOnly, int min_length , int max_gap , bool call_summits , bool auto_cutoff )
{
        /*"""Call peaks for all chromosomes. Return a PeakIO object.
        
        scoring_function_s: symbols of functions to calculate score. 'p' for pscore, 'q' for qscore, 'f' for fold change, 's' for subtraction. for example: ['p', 'q']
        score_cutoff_s    : cutoff values corresponding to scoring functions
        min_length        : minimum length of peak
        max_gap           : maximum gap of 'insignificant' regions within a peak. Note, for PE_mode, max_gap and max_length are both set as fragment length.
        call_summits      : boolean. Whether or not call sub-peaks.
        save_bedGraph     : whether or not to save pileup and control into a bedGraph file
        """*/
    
    
    //std::string chrom;
    std::string s;
    //bytes tmp_bytes;

    //prepare p-q table
    if (pqtable.size() == 0 && uniqOnly){
        info("#4 Pre-compute pvalue-qvalue table...");
        if (auto_cutoff){
            info("#4 Cutoff will be automatically decided!");
            //ADD __pre_computes( max_gap, min_length );
        }
        else{
            __cal_pvalue_qvalue_table();
        }
    }

    // prepare bedGraph file
 /*ADD   if (save_bedGraph){
        bedGraph_treat_f.open( bedGraph_treat_filename, std::ofstream::out);
        bedGraph_ctrl_f.open( bedGraph_control_filename, std::ofstream::out );

        info ("#3 In the peak calling step, the following will be performed simultaneously:");
        info ("#3   Write bedGraph files for treatment pileup (after scaling if necessary)... " + bedGraph_filename_prefix + "_treat_pileup.bdg");
        info ("#3   Write bedGraph files for control lambda (after scaling if necessary)... " +bedGraph_filename_prefix + "_control_lambda.bdg");

        if (save_SPMR){
            info ( "#3   --SPMR is requested, so pileup will be normalized by sequencing depth in million reads." );
        }
        else {
            if (treat_scaling_factor == 1){
                info ( "#3   Pileup will be based on sequencing depth in treatment." );
            }
            else {
                info ( "#3   Pileup will be based on sequencing depth in control." );
            }
        }

        if (trackline){
                // this line is REQUIRED by the wiggle format for UCSC browser
            std::string tmp_bytes = ("track type=bedGraph name=\"treatment pileup\" description=\"treatment pileup after possible scaling for \'" + bedGraph_filename_prefix + " \'\"\n");
            bedGraph_treat_f << tmp_bytes ;
            tmp_bytes = ("track type=bedGraph name=\"control lambda\" description=\"control lambda after possible scaling for " + "\'" + bedGraph_filename_prefix + "\'\"\n");
            bedGraph_ctrl_f << tmp_bytes ;
        }
    }
*/
    info("#4 Call peaks for each chromosome...");
    for ( auto chrom : chromosomes){
            // treat/control bedGraph will be saved if requested by user.
        __chrom_call_peak_using_certain_criteria ( peaks, chrom, scoring_function_symbols, score_cutoff_s, min_length, max_gap, uniqOnly,call_summits, save_bedGraph );
    }

    // close bedGraph file
    if (save_bedGraph){
        bedGraph_treat_f.close();
        bedGraph_ctrl_f.close();
        save_bedGraph = false;
    }
}

bool CallerFromAlignments::__close_peak_wo_subpeaks (bool uniqOnly,std::vector<peak_content_t> peak_content, PeakIO *peaks, int min_length, std::string chrom, int smoothlen, std::vector<std::vector<double> > score_array_s, std::vector<double> score_cutoff_s)
{
       /* """Close the peak region, output peak boundaries, peak summit
        and scores, then add the peak to peakIO object. peak_content contains [start, end, treat_p, ctrl_p, index_in_score_array], peaks: a PeakIO object

        """*/

    int peak_length = peak_content[peak_content.size()-1].e - peak_content[0].s;
    
    if (peak_length >= min_length)
    { // # if the peak is too small, reject it
        //debug("peak length greater than min length");
        std::vector<int> tsummit;
        std::vector<int> tsummit_index;
        int summit_pos   = 0;
        double summit_value = 0.0;
        
       
       // debug("in __close_peak_wo_subpeaks");
        for(size_t i=0; i< peak_content.size();i++)
        {
            int tstart = peak_content[i].s;
            int tend = peak_content[i].e;
            double ttreat_p = peak_content[i].t_val;
           // double tctrl_p = peak_content[i].c_val;
            //int tlist_scores_p = peak_content[i].index;
            
           // debug("tstart = " + std::to_string(tstart) + "\t" + std::to_string(tend));
            double tscore = ttreat_p; // # use qscore as general score to find summit
            if(summit_value == 0 or summit_value < tscore){
                tsummit.clear();
                tsummit_index.clear();
                tsummit.push_back((tend + tstart) / 2);
                tsummit_index.push_back(i);
                summit_value = tscore;
            }
            else {
                if (summit_value == tscore){
                    //# remember continuous summit values
                    tsummit.push_back(int((tend + tstart) / 2));
                    tsummit_index.push_back( i );
                }
            }
            
        }
            // the middle of all highest points in peak region is defined as summit
        int midindex = int((tsummit.size() + 1) / 2) - 1;
        summit_pos    = tsummit[ midindex ];
        int summit_index  = tsummit_index[ midindex ];
        
        
        double summit_treat = peak_content[ summit_index ].t_val;
        double summit_ctrl = peak_content[ summit_index ].c_val;
        
        //# this is a double-check to see if the summit can pass cutoff values.
        for (size_t i=0 ; i < score_cutoff_s.size();i++)
        {
            
            if (score_cutoff_s[i] > score_array_s[ i ][ peak_content[ summit_index ].index ])
            {
                return false; // # not passed, then disgard this peak.
            }
        }
        
        
        double summit_p_score = get_pscore( int(summit_treat), summit_ctrl );
        //double summit_q_score ;
       // std::string key_str = lexical_cast<std::string>(summit_p_score) ;
       // std::hash<std::string> str_hash;
       // size_t key_value =  str_hash( key_str );
        //double rounded_up = (ceil)(summit_p_score * 100000) / 100000;
        
        int rounded_up = (int)(ceil)(summit_p_score * 100000);
        
        double summit_q_score = 0;
        if (uniqOnly) {
            
            summit_q_score = pqtable[ rounded_up ];

        }
        
        
        
        peaks->add( chrom,           // chromosome
                       peak_content[0].s, // start
                       peak_content[peak_content.size()-1].e, // end
                       summit_pos, // summit position
                       summit_q_score, // score at summit
                       summit_treat, // pileup
                       summit_p_score, // pvalue
                       1.0* ( summit_treat + pseudocount ) / (summit_ctrl + pseudocount ), // fold change
                       summit_q_score // qvalue
                      );
        //# start a new peak
        return true;
    }
    else {
        
        return false;
    }
}

bool CallerFromAlignments::__close_peak_with_subpeaks (std::vector<peak_content_t> peak_content, PeakIO *peaks, int min_length, std::string chrom, int smoothlen, std::vector<std::vector<double> > score_array_s, std::vector<double> score_cutoff_s,float min_valley )
{
       /* """Algorithm implemented by Ben, to profile the pileup signals
        within a peak region then find subpeak summits. This method is
        highly recommended for TFBS or DNAase I sites.
        
        """*/

    int peak_length = peak_content[ peak_content.size()-1].e - peak_content[ 0 ].s;
    
    if (peak_length < min_length){
        return false;
    }// if the region is too small, reject it

    //Add 10 bp padding to peak region so that we can get true minima
    int end = peak_content[peak_content.size()-1].e + 10;
    int start = peak_content[0].s - 10;
    int start_boundary = 10; //// # this is the offset of original peak boundary in peakdata list.

    
    if (start < 0){
        start_boundary = 10 + start;// # this is the offset of original peak boundary in peakdata list.
        start = 0;
    }

    std::vector<double> peakdata(end-start,0); // # save the scores (qscore) for each position in this region
    std::vector<int> peakindices(end-start,0); // # save the indices for each position in this region
    for (size_t i=0;i < peak_content.size();i++)
    {
        //(tstart, tend, ttreat_p, tctrl_p, tlist_scores_p) = peak_content[i];
        //tscore = self.pqtable[ get_pscore(int(ttreat_p), tctrl_p) ] # use qscore as general score to find summit
        int tscore = peak_content[i].t_val;// # use pileup as general score to find summit
        int m = peak_content[i].s - start + start_boundary;
        int n = peak_content[i].e - start + start_boundary;
        for (int k = m; k< n; k++) {
            peakdata[k] = tscore;
            peakindices[k] = i;
        }
    }
    std::vector<int> summit_offsets;
    maxima(peakdata,summit_offsets,smoothlen); // # offsets are the indices for summits in peakdata/peakindices array.
        //print "maxima:",summit_offsets
    if (summit_offsets.size() == 0){
            // **failsafe** if no summits, fall back on old approach #
        return __close_peak_wo_subpeaks(true,peak_content, peaks, min_length, chrom, smoothlen, score_array_s, score_cutoff_s);
    }
    else{
        // remove maxima that occurred in padding
        //m = np.searchsorted(summit_offsets, start_boundary);
        std::sort(summit_offsets.begin(),summit_offsets.end());
        std::vector<int>::iterator m = std::lower_bound(summit_offsets.begin(),summit_offsets.end(),start_boundary);
        std::vector<int>::iterator n = std::upper_bound(summit_offsets.begin(),summit_offsets.end(), peak_length + start_boundary);
        
        
        summit_offsets.erase (n+1,summit_offsets.end());
        summit_offsets.erase(summit_offsets.begin(),m-1);
        //summit_offsets = summit_offsets[m:n];
    }
    
    summit_offsets = enforce_peakyness(peakdata, summit_offsets);
    
        //print "enforced:",summit_offsets
    if (summit_offsets.size() == 0){
            // **failsafe** if no summits, fall back on old approach #
        return __close_peak_wo_subpeaks(true,peak_content, peaks, min_length, chrom, smoothlen, score_array_s, score_cutoff_s);
    }
    
    for (size_t k=0; k < summit_offsets.size(); k++)
    {
        int summit_index = peakindices[summit_offsets[k]];
        int summit_offset = summit_offsets[k] - start_boundary;

        double summit_treat = peak_content[ summit_index ].t_val;
        double summit_ctrl = peak_content[ summit_index ].c_val;

        double summit_p_score = get_pscore( int(summit_treat), summit_ctrl );
        //std::string key_str = lexical_cast<std::string>(summit_p_score) ;
        
        //std::hash<std::string> str_hash;
        //size_t key_value =  str_hash( key_str );
        
        
        //double summit_q_score = 0; // pqtable[ key_value ];
        
        //ADD double rounded_up = (ceil)(summit_p_score * 100000) / 100000;
        int rounded_up = (int)(ceil)(summit_p_score * 100000);
        double summit_q_score= pqtable[ rounded_up ];
        
        for(size_t i=0; i < score_cutoff_s.size();i++)
        {
            if (score_cutoff_s[i] > score_array_s[ i ][ peak_content[ summit_index ].index ]){
                return false; // # not passed, then disgard this summit.
            }
        }

        peaks->add( chrom,
                       peak_content[ 0 ].s,
                       peak_content[peak_content.size()-1].e,
                       start + summit_offset,
                       summit_q_score,
                       summit_treat,
                       summit_p_score,
                       1.0* ( summit_treat + pseudocount ) / ( summit_ctrl + pseudocount ), //# fold change
                       summit_q_score
                      );
    }
    //start a new peak
    return true;
}

std::vector<double> CallerFromAlignments::__cal_qscore ( std::vector<double> array1, std::vector<double> array2 )
{
    std::vector<double> s;

    if(array1.size() != array2.size())
    {
        error("tow lists should have the same length.");
        std::exit(1);
    }
    
    //debug("pqtable size " + std::to_string(pqtable.size()));
    for (size_t i=0; i< array1.size(); i++)
    {
        double p = get_pscore( int(array1[i]), array2[i] );
        
        
        //std::string key_str = lexical_cast<std::string>(p) ;
        
        //std::hash<std::string> str_hash;
        //size_t key_value =  str_hash( key_str );
        
       //ADD double rounded_up = (ceil)(p * 100000) / 100000;
        int  rounded_up = (int)(p * 100000);
        
        double q = 0;

       // stringstream(q_str) >> q;
        if (pqtable.find(rounded_up) == pqtable.end()) {
            int k1 = rounded_up + 1;
            
            if (pqtable.find(k1) == pqtable.end()) {
                int k2 = rounded_up - 1;
                if (pqtable.find(k2) == pqtable.end()) {
                   // debug("didn't find matching qvalue in pqtable " + std::to_string(rounded_up) + "\t" + std::to_string(k1) + "\t" + std::to_string(k2));
                    q = 0;
                }
                else {
                    q = pqtable[k2];
                }
            }
            else {
                q = pqtable[k1];
            }
            
            
        }
        else {
           // q = pqtable[rounded_up];
            q = pqtable[rounded_up];
        }
        
       
        
       // debug("p = " + std::to_string(p) + " q = " + std::to_string(q));
        s.push_back(q);

    }
        
    
   
    return s;
}


void CallerFromAlignments::__chrom_call_peak_using_certain_criteria ( PeakIO *peaks, std::string chrom, std::vector<std::string> scoring_function_s, std::vector<double> score_cutoff_s, int min_length,int max_gap, bool uniqOnly, bool call_summits, bool save_bedGraph )
{
        /*""" Call peaks for a chromosome.

        Combination of criteria is allowed here.

        peaks: a PeakIO object, the return value of this function
        scoring_function_s: symbols of functions to calculate score as score=f(x, y) where x is treatment pileup, and y is control pileup
        save_bedGraph     : whether or not to save pileup and control into a bedGraph file
        """*/

    //list score_array_s ;         //# list to keep different types of scores
    std::vector<peak_content_t> peak_content_list  ;         //#  to store information for a
                                       /* #  chunk in a peak region, it
                                        #  contains lists of: 1. left
                                        #  position; 2. right
                                        #  position; 3. treatment
                                        #  value; 4. control value;
                                        #  5. list of scores at this
                                        #  chunk*/
    //long  lastp;
    //float tp, cp;
    
    if (scoring_function_s.size() != score_cutoff_s.size()) {
        error("number of functions and cutoffs should be the same!");
        std::exit(1);
    }
        
    //peak_content = [] ; //          # to store points above cutoff

    //# first, build pileup, chr_pos_treat_ctrl
    //# this step will be speeped up if pqtable is pre-computed.
    //debug("In __chrom_call_peak_using_certain_criteria ");
    //debug("score_cutoff size " + std::to_string(score_cutoff_s.size()));
    //debug("score cutoff value : " + std::to_string(score_cutoff_s[0]));
    
    __pileup_treat_ctrl_a_chromosome( chrom,uniqOnly );
    
    //debug("In __chrom_call_peak_using_certain_criteria  size of chr_pos_treat_ctrl " + std::to_string(chr_pos_treat_ctrl.pos.size()));
    

    //# while save_bedGraph is true, invoke __write_bedGraph_for_a_chromosome
    if (save_bedGraph){
       // __write_bedGraph_for_a_chromosome ( chrom );
        info("write bedgraph for a chromosome.");
    }
    
    
    //keep all types of scores needed
    //t0 = ttime()
    std::vector<std::vector<double> > score_array_s ;
    for( size_t i=0;i< scoring_function_s.size();i++)
    {
        if (scoring_function_s[i] == "p"){
            //debug("scoring function p");
            score_array_s.push_back(__cal_pscore( chr_pos_treat_ctrl.treat_v, chr_pos_treat_ctrl.ctrl_v ));
        }
        if(scoring_function_s[i] == "q")
        {
            //debug("scoring function q");
            score_array_s.push_back(__cal_qscore( chr_pos_treat_ctrl.treat_v, chr_pos_treat_ctrl.ctrl_v )) ;
            
        }
     //ADD   if(scoring_function_s[i] == "f")
     //ADD   {
     //ADD       score_array_s.push_back(__cal_FE( chr_pos_treat_ctrl.treat_v, chr_pos_treat_ctrl.ctrl_v ) );
     //ADD   }
    //ADD    if(scoring_function_s[i] == "s")
    //ADD    {
    //ADD        score_array_s.push_back(__cal_subtraction( chr_pos_treat_ctrl.treat_v, chr_pos_treat_ctrl.ctrl_v ));
    //ADD    }
    }
    
    // get the regions with scores above cutoffs
    //above_cutoff = np.nonzero( apply_multiple_cutoffs(score_array_s,score_cutoff_s) )[0];// # this is not an optimized method. It would be better to store score array in a 2-D ndarray?
    //above_cutoff_index_array = np.arange(pos_array.shape[0],dtype="int32")[above_cutoff];// # indices
    //above_cutoff_endpos = pos_array[above_cutoff];// # end positions of regions where score is above cutoff
    //above_cutoff_startpos = pos_array[above_cutoff-1];// # start positions of regions where score is above cutoff

    std::vector<int> above_cutoff_startpos, above_cutoff_endpos;
    std::vector<int> above_cutoff_index_array;
    
    //print "time to build peak regions: %f" % self.test_time
    
    //return peaks;

    
    apply_multiple_cutoffs(score_array_s,score_cutoff_s,above_cutoff_index_array);
    
    if (above_cutoff_index_array.size() == 0){
            //# nothing above cutoff
        //debug("nothing above cutoff");
        return ;
    }
    
    for (size_t j =0; j < above_cutoff_index_array.size(); j++) {
            above_cutoff_endpos.push_back(chr_pos_treat_ctrl.pos[above_cutoff_index_array[j]]);
            above_cutoff_startpos.push_back(chr_pos_treat_ctrl.pos[above_cutoff_index_array[j]-1]);
    }

    if (above_cutoff_index_array[0] == 0){
        //# first element > cutoff, fix the first point as 0. otherwise it would be the last item in data[chrom]['pos']
        above_cutoff_startpos[0] = 0;
        
    }
   

    //#print "apply cutoff -- chrom:",chrom,"  time:", ttime() - t0
    //# start to build peak regions
    //#t0 = ttime()

    //# first bit of region above cutoff

    peak_content_t pt ;
    
    pt.s =above_cutoff_startpos[0];
    pt.e =above_cutoff_endpos[0];
    pt.index = above_cutoff_index_array[0];
    pt.t_val = chr_pos_treat_ctrl.treat_v[pt.index];
    pt.c_val = chr_pos_treat_ctrl.ctrl_v[pt.index];
    
    peak_content_list.push_back(pt);
    
    long lastp = pt.e;
    
    
    for (size_t i =1; i < above_cutoff_startpos.size(); i++ )
    {
        peak_content_t pt2;
        
        pt2.s = above_cutoff_startpos[ i ];
        pt2.e = above_cutoff_endpos[ i ];
        pt2.index = above_cutoff_index_array[ i ];
       // acs_ptr += 1;
        //ace_ptr += 1;
        //acia_ptr+= 1;
        pt2.t_val = chr_pos_treat_ctrl.treat_v[ pt2.index ];
        pt2.c_val = chr_pos_treat_ctrl.ctrl_v[ pt2.index ];
        long tl = pt2.s - lastp;
        
        if (tl <= max_gap){
                // append.
            peak_content_list.push_back(pt2 );
            lastp = pt2.e; // #above_cutoff_endpos[i]
        }
        else{
                //close
            if (call_summits){
                __close_peak_with_subpeaks (peak_content_list, peaks, min_length, chrom, min_length, score_array_s, score_cutoff_s ); // # smooth length is min_length, i.e. fragment size 'd'
            }
            else{
                
                __close_peak_wo_subpeaks (uniqOnly,peak_content_list, peaks, min_length, chrom, min_length, score_array_s, score_cutoff_s ); //# smooth length is min_length, i.e. fragment size 'd'
             }//end else
            peak_content_list.clear();
            peak_content_list.push_back(pt2);
            lastp = pt2.e;// #above_cutoff_endpos[i];
        }//end else
    } //end for
    
    //save the last peak
             
    if (peak_content_list.size()==0){
        return ;
    }
    else{
        if (call_summits){
            __close_peak_with_subpeaks (peak_content_list, peaks, min_length, chrom, min_length, score_array_s, score_cutoff_s ); // # smooth length is min_length, i.e. fragment size 'd'
        }
        else{
            __close_peak_wo_subpeaks (uniqOnly,peak_content_list, peaks, min_length, chrom, min_length, score_array_s, score_cutoff_s ); //# smooth length is min_length, i.e. fragment size 'd'
        }
    }

    //#print "close peaks -- chrom:",chrom,"  time:", ttime() - t0
   // return peaks;
}
