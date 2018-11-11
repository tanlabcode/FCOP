///////////////////////////////////////////////////////////////////////////////
//  Written by  : Li Teng, The University of Iowa, li-teng@uiowa.edu 
//  Date    : Aug. 2013
///////////////////////////////////////////////////////////////////////////////

#include <time.h>
#include <stdio.h>
#include <unistd.h>
#include <fstream>
#include <assert.h>
#include <iostream>
#include <math.h>
#include <cstdlib>
#include <fcntl.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <stdio.h>
#include <algorithm>
#include <set>
#include <vector>
#include <map>
#include <ctime>
#include <cstring>

using namespace std;

#define DEBUG
#undef DEBUG

#define MAX_LINE	2000
#define NAME_LEN	2000 

int NUMBER_GENE, NUMBER_TF;
char CUTFILE[100]; 

int S, G, E, M; // input parameters 
float P; // input parameters 
float *data; 

//============ Search Scheme 3_2 =======================================
void recurSearch3_2(vector<int> head, vector<int> tail, vector< pair< vector<int>, pair<vector<int>, float> > >* gsSet, int minS, int minE, float minP, vector< vector<pair<int,int> > >* invList, vector<pair<int,int> > ILintsec, float *EnhancerP);

vector< pair< vector<int>, pair<vector<int>, float> > > SampleGeneSearch3_2(int size, int minS, int minE, float minP, vector< vector<pair<int,int> > >* invList, float *EnhancerP)
{
    vector< pair< vector<int>, pair<vector<int>, float> > > gsSet;
    for (int n=0; n<= size-minS; n++)
    {
        vector<int> head (1, n);
//        cout<<"Processing head["<<n<<"]"<<endl;
                
        //get tail
        vector<int> tail;
        for (int m = n+1; m < size; m++)
        {
            if ((*invList)[m].size() >= (int) minE)

                tail.push_back(m);
        }
        recurSearch3_2(head, tail, &gsSet, minS, minE, minP, invList, (*invList)[n], EnhancerP);
    }
    return gsSet;
}
    
void recurSearch3_2(vector<int> head, vector<int> tail, vector< pair< vector<int>, pair< vector<int>, float > > >* gsSet, int minS, int minE, float minP, vector< vector<pair<int,int> > >* invList, vector<pair<int,int> > ILintsec, float *EnhancerP)
{
    vector<pair<int,int> > g ((*invList)[head[0]]);
    
    //get the intersection
    int n = head.size()-1; // get the tail of the head, the newest element that was added to the head
    //vector<pair<int,int> > tempG (g);
    //g.resize(min(tempG.size(), (*invList)[head[n]].size()));
	
    g.resize(set_intersection(ILintsec.begin(), ILintsec.end(), (*invList)[head[n]].begin(), (*invList)[head[n]].end(), g.begin())-g.begin());
    
    vector<int> gSet;
    for(int n=0; n<g.size(); n++)
    {
        gSet.push_back(g[n].first);
    }
    gSet.resize(unique(gSet.begin(), gSet.end())-gSet.begin());

    //check the total length
    if(head.size()+tail.size() < minS)
        return;
 
    if(gSet.size()<minE)
        return;
    
    while(tail.size() != 0)
    {
        vector<int> headTemp = head;
        headTemp.push_back(tail[0]);
        tail.erase(tail.begin());
        recurSearch3_2(headTemp, tail, gsSet, minS, minE, minP, invList, g, EnhancerP);
    }

    //check the number of possible supports	
    if((head.size()<=minS) && (head.size() >=2 )) // find Cuts for all the patterns
//    if(head.size()==minS) // only for tunning
    {
        int T = gSet.size();
        float pij;	

        float *p1, *p2, *pi;
        p1 = new float [T + 1 - minE];
        p2 = new float [T + 1 - minE];
        pi = new float [T];

        for(int n = 0; n< T; n++){
            float ei= 1;
    	    for(int m =0; m <head.size(); m++){
	        ei = ei * data[ NUMBER_TF * gSet.at(n) + head.at(m)];		  
		//cout<<data[ NUMBER_TF * gSet.at(n) + head.at(m)]<<"\t";
  	    }	   
	    pi[n] = ei * EnhancerP[gSet.at(n)];  // probability of current transaction support itemset
	}

        int maxx;

        maxx = T;
        int flag = 1;
        int kk = minE;
        while(flag && (kk<=maxx)){
  	    for(int j=0; j < T +1 -kk; j++)
	        p1[j] = 1;

	    for(int i =0; i<kk; i++){
	        for(int j =0; j< T+1-kk ; j++){
	            if(j==0){
		        p2[0] = pi[i + j ] * p1[j] + (1 - pi[i + j]) * 0;
		    }				
		    else{
		        p2[j] = pi[i + j] * p1[j] + (1-pi[i + j]) * p2[j-1];
		    }
	        }		
	        for (int j =0; j< T+1-kk; j++)
		    p1[j] = p2[j];	
	    }
            pij = p2[T-kk] ;

            if(flag && (pij < minP)){   // the probability of finding the pattern in random data is less than a threshold 
               flag = 0;
               pij = kk;
            }
            else{
               kk++;   
            }

        }
        cout << kk <<"\t";
        delete[] p1;
        delete[] p2;
        delete[] pi;
 
        for(int i =0 ; i< head.size();i++)
           cout<<head.at(i)<<"\t";
        cout<<endl;

        pair<vector<int>, float> head2( head, pij);
	pair< vector<int>, pair<vector<int>, float> >  gs(gSet, head2); 
	(*gsSet).push_back(gs);  
    }

}

//============ Search Scheme 4 =======================================
void recurSearch4(vector<int> head, vector<int> tail, vector< pair< vector<int>, pair<vector<int>, float> > >* gsSet, int minS, int minE, float minP, vector< vector<pair<int,int> > >* invList, vector<pair<int,int> > ILintsec, float *EnhancerP, vector< pair< int , pair< int, vector<int> > > >* Cutoffs);

vector< pair< vector<int>, pair<vector<int>, float> > > SampleGeneSearch4(int size, int minS, int minE, float minP, vector< vector<pair<int,int> > >* invList, float *EnhancerP)
{
    // Read Cutoff for patterns from file Cuts_HepG2
    vector< pair< int , pair< int, vector<int> > > >Cutoffs;   // < cutoff, <length, <pattern> > >
	
    char buf[MAX_LINE];
	int len;
        ifstream infile;
	infile.open(CUTFILE);
	
	if(!infile) {
        cout << "Cannot open file CUTFILE.\n";
        exit (-1);
    }
    int flag;
    while(infile.getline(buf, MAX_LINE, '\n')){  // delim defaults to '\n'
        len = strlen(buf);
				
	int patterncutoff = 0;
	int valuenum = 0;
	int currentvalue = 0;
	vector<int> pattern;        
	flag = 1;	
//        cout<<len<<";"<<buf<<";"<<endl;

	for(int i =0; i<len ; i++){
	    if( buf[i] > 57 || buf[i] <48){
	        if(valuenum == 0){ //first value read
                    if (currentvalue > minE){
                        patterncutoff = currentvalue;  //only read cutoff when it's larger than minE
                    }
                    else{
                        flag = 0;
                        i = len;
                    }                    
		}
		else{
		    pattern.push_back(currentvalue);	
		}
		valuenum = valuenum + 1;
		currentvalue = 0;
	    }
	    else{
		currentvalue = currentvalue * 10 + buf[i] - 48;
	    }
	}
        if( buf[len-1] != '\n'){  // based on different platform, change row symbol would be different
           valuenum = valuenum + 1;
           pattern.push_back(currentvalue);
        }

        if(flag){
	    pair< int, vector<int> > lpattern(valuenum - 1, pattern);
  	    pair< int, pair<int, vector<int> > > clpattern(patterncutoff, lpattern) ; 
	    Cutoffs.push_back(clpattern);
        }
    }
    infile.close();
//    cout<<"Read Cutoff for "<<Cutoffs.size()<<" patterns."<<endl;
		
    vector< pair< vector<int>, pair<vector<int>, float> > > gsSet;
    for (int n=0; n<= size-minS; n++)
    {
        vector<int> head (1, n);
//        cout<<"Processing head["<<n<<"]"<<endl;
                
        //get tail
        vector<int> tail;
        for (int m = n+1; m < size; m++)
        {
            if ((*invList)[m].size() >= (int) minE)

                tail.push_back(m);
        }
        recurSearch4(head, tail, &gsSet, minS, minE, minP, invList, (*invList)[n], EnhancerP, &Cutoffs);
    }
    return gsSet;
}
    
void recurSearch4(vector<int> head, vector<int> tail, vector< pair< vector<int>, pair< vector<int>, float > > >* gsSet, int minS, int minE, float minP, vector< vector<pair<int,int> > >* invList, vector<pair<int,int> > ILintsec, float *EnhancerP, vector< pair< int , pair< int, vector<int> > > >* Cutoffs)
{
    vector<pair<int,int> > g ((*invList)[head[0]]);
    
    //get the intersection
    int n = head.size()-1; // get the tail of the head!! the newest element that was added to the head
    //vector<pair<int,int> > tempG (g);
    //g.resize(min(tempG.size(), (*invList)[head[n]].size()));
	
    g.resize(set_intersection(ILintsec.begin(), ILintsec.end(), (*invList)[head[n]].begin(), (*invList)[head[n]].end(), g.begin())-g.begin());
    
    vector<int> gSet;
    for(int n=0; n<g.size(); n++)
    {
        gSet.push_back(g[n].first);
    }
    gSet.resize(unique(gSet.begin(), gSet.end())-gSet.begin());

    //check the total length
    if(head.size()+tail.size() < minS)
        return;
 
    if(gSet.size()<minE)
        return;
    
    while(tail.size() != 0)
    {
        vector<int> headTemp = head;
        headTemp.push_back(tail[0]);
        tail.erase(tail.begin());
        recurSearch4(headTemp, tail, gsSet, minS, minE, minP, invList, g, EnhancerP, Cutoffs);
    }

    //check the number of possible supports	
    if(head.size()>=minS)
//    if(head.size()==minS) // only for tunning
    {
	    // in Cutoffs find the corresponding cutoff, if not exist using minE
		//cout<<"Check in Cutoffs for "<<(*Cutoffs).size()<<" patterns..."<<endl;
		int cutoff = minE;
		for(int i = 0; i<(*Cutoffs).size(); i++){
		   //cout<< "cutoff = "<<(*Cutoffs)[i].first << ", length = "<<(*Cutoffs)[i].second.first <<", pattern[0] = "<<(*Cutoffs)[i].second.second.at(0)<<endl;
		    if (head == (*Cutoffs)[i].second.second){
			   //cout<< head.at(0)<<","<<head.at(1)<<"vs"<<(*Cutoffs)[i].second.second.at(0)<<(*Cutoffs)[i].second.second.at(1)<<endl;
			    if((*Cutoffs)[i].first > minE)
			        cutoff = (*Cutoffs)[i].first;
		    }
		}
		
        int T = gSet.size();
        float pij;	

        float *p1, *p2, *pi;
        p1 = new float [T + 1 - minE];
        p2 = new float [T + 1 - minE];
        pi = new float [T];

        for(int n = 0; n< T; n++){
            float ei= 1;
    	    for(int m =0; m <head.size(); m++){
	        ei = ei * data[ NUMBER_TF * gSet.at(n) + head.at(m)];		  
		//cout<<data[ NUMBER_TF * gSet.at(n) + head.at(m)]<<"\t";
  	        }	   
	        pi[n] = ei * EnhancerP[gSet.at(n)];  // probability of current transaction support itemset
	    }

        for(int j=0; j < T +1 -cutoff; j++)
	        p1[j] = 1;

	    for(int i =0; i<cutoff; i++){
	        for(int j =0; j< T+1-cutoff ; j++){
	            if(j==0){
	                p2[0] = pi[i + j ] * p1[j] + (1 - pi[i + j]) * 0;
		        }				
		        else{
		            p2[j] = pi[i + j] * p1[j] + (1-pi[i + j]) * p2[j-1];
		        }
	        }		
	        for (int j =0; j< T+1-cutoff; j++)
		        p1[j] = p2[j];	
	    }
  	    //cout << p2[T - minE]<<endl;
        pij = p2[T-cutoff] ;
   
        delete[] p1;
        delete[] p2;
        delete[] pi;

	if(pij >= minP){

	    for(int n=0; n<(*gsSet).size(); n++){ //*******************trim
                if ((*gsSet)[n].second.first.size()>=head.size()){
		    if(includes((*gsSet)[n].second.first.begin(), (*gsSet)[n].second.first.end(), head.begin(), head.end()))
                    {
                        return; // current head is a found frequent itemset or is a subset of it
                    }  
		}
            }
	      
            pair<vector<int>, float> head2( head, pij);
            pair< vector<int>, pair<vector<int>, float> >  gs(gSet, head2); 
	    (*gsSet).push_back(gs);  
	}
    }
}



int main(int argc, char *argv[])
{
    char   filename[100];
    char   filenameTF[100];
    char   fileEnhancerP[100];
    
    int    op1_1=0, op1_2=0, op1_3=0, op1_4=0, op1_5=0, op1_6 =0, op1_7=0;
    int    option1 = 0;
    int    option2 = 0;
    int    option3 = 0;
    ofstream out;
	
    if(argc < 6) {
        cout << endl << "Usage: " << argv[0] << " [other options]" << endl;
	cout << "==================== Inputes ==========================" << endl;
	cout << "-f     -ffilename      File for input data" << endl;
	cout << "-X     -X5000         No. of rows of input data (genes/enhancers)" << endl;
	cout << "-Y     -Y62           No. of columns of input data (TFs)" << endl<<endl;
	cout << "-C     -CCuts_GM12878 Pattern specific cutoff file" << endl;
	cout << "-s     -s5            Pattern length" << endl;
	cout << "-m     -m80           Least support"<<endl;
	cout << "-p     -p0.8          Pvalue cutoff in Step1; Minimal frequentness probability cutoff in Step2" << endl;
	cout << "-------------------- Optional Inputes ------------------"<<endl;
	cout << "-t     -tTFlist    File for the TFs" << endl;
	cout << "-h     -hEnhancerProbabilities  Files for the enhancer confidence levels" <<endl; 
	cout << "==================== Examples ==========================" << endl;
        cout << "Step1: Generate Pattern Specific Cutoff at pvalue = 0.001" <<endl;
	cout << " -fmapping_GM12878 -X5000 -Y62 -s5 -p0.001  -hEnhancerP  > Cuts_GM12878" << endl<<endl;
        cout << "Step2: Search with Frequentness Probability on Minimum Support and Pattern Specific Cutoff" <<endl;
	cout << " -fmapping_GM12878 -X5000 -Y62 -CCuts_GM12878 -m80 -p0.8 " << endl;
	cout << " -fmapping_GM12878 -X5000 -Y62 -CCuts_GM12878 -m80 -p0.8 -tTFfile -hEhancerP" << endl;
        cout << "========================================================" << endl << endl;
	exit(1);
    }	
    for(int i=1; i<argc; i++){
        if( argv[i][0] == '-' )
	    switch( argv[i][1] ){
	        case 'f':  	
                    strcpy( filename, argv[i]+2 ); op1_1=1;       break;
	        case 'X':
	            NUMBER_GENE = (int) atoi(argv[i]+2); op1_2=1; break;
                case 'Y':
	            NUMBER_TF = (int) atoi(argv[i]+2); op1_3 =1; break;
	        case 'C':  	
                    strcpy( CUTFILE, argv[i]+2 ); op1_4=1;       break;
		case 'm':
		   M = (int) atoi(argv[i]+2); op1_5 = 1; break;
		case 'p':
	           P = (float) atof(argv[i]+2); op1_6 = 1; break;
                case 's':
                   S = (int) atoi(argv[i]+2); op1_7 = 1; break;
  	        case 't':  	
                   strcpy( filenameTF, argv[i]+2 ); option2=1;       break;  //TF list available
		case 'h':
		   strcpy(fileEnhancerP, argv[i]+2); option3 =1; break;

	        default: cout << "Error Input Options: " << argv[i] << endl; exit(0); break;
	    }
        else{	cout << "Error Input Options: " << argv[i] << endl;  exit(0); }
    }	
	
    //cout <<op1_1<< ", "<<op1_2<< ", "<<op1_3<< ", "<<op1_4<< ", "<<op1_5<< ", "<<endl;
    if( op1_1 & op1_2 & op1_3 & !op1_4 & ! op1_5 & op1_6 & op1_7 ){
        option1 = 1;
    }
    else if (op1_1 & op1_2 & op1_3 & op1_4 & op1_5 & op1_6 & !op1_7 ) {
        option1 = 4;
    }
    else{
        cout <<"Inputs should be such combinations:"<<endl;
        cout <<"-f -X -Y -s -p  Search with Frequentness Probability on Least Support and Pattern Specific Cutoff"<<endl;
 	cout <<"-f -X -Y -C -m -p    Search with Frequentness Probability on Least Support and Pattern Specific Cutoff"<<endl;
        cout <<"-t and -h are optional"<<endl;
	exit(0); 
    }
	
    // start clock
    time_t start,end;
    double dif;
    time (&start);

    char    TFlist[100][20];
    float   EnhancerP[NUMBER_GENE];
    
    
    vector< vector<int> > maxCliq[NUMBER_GENE];
    //=====================================================
    //Read maxCliq from input file !! 
    ifstream infile;
    infile.open(filename);
	
    vector<int> head;
    data = new float[NUMBER_GENE*NUMBER_TF];
	
    for(int i = 0; i < NUMBER_GENE; i++)
    {	
        char buf[200000],sub[2000];
	int k, p, columnID, len;
	float ii;
	infile.getline(buf, 200000, '\n');
//	printf("%s", buf);  getchar();	
	
        len = strlen( buf );
        assert(len <= 200000);    
        p=0;  //pointer index for sub
	k= 0; //pointer index for buf
	columnID = 0;  //index of current column
	sub[0] = '\0';
    
        while( buf[k] != '\n' && buf[k]!=0xd && k <= len){
//	cout<<buf[k]<<".";
 	    if (buf[k] == '\t' || buf[k] == ' ')
            {
	        sub[p] = '\0';
		ii = atof(sub);
		data[i * NUMBER_TF + columnID] = ii;  //load raw data in data
		if(ii>0){
		    head.push_back(columnID);
		    //cout<<"read "<<ii<<",";
		}				   
                sub[0] = '\0';
                k++;
		p=0;
		columnID++;
            }
            else
            {
		sub[p] = buf[k];
		p++;
		k++;				
            }
	}		
	sub[p] = '\0';
	ii = atof(sub);
	data[i*NUMBER_TF + columnID] = ii; //load raw data in data
	if(ii>0)
	    head.push_back(columnID);
	/*cout<<"Push row "<<i<< " into maxCliq:";
	for(int m =0; m <head.size(); m++)
	    cout<<head.at(m)<<",";
        cout<<endl;*/		
		
	if(head.size() >= S){
	    maxCliq[i].push_back(head);
	}
	/*else{  // trim the genes that does not support any itemset longer than S
	    //vector<int> emptyhead;
	    //maxCliq[i].push_back(emptyhead);
	}*/
	head.clear(); 
	#ifdef DEBUG
	    if(maxCliq[i].size()>0)
		cout<<"Read row "<< i+1 <<" for "<< maxCliq[i].size()<< " maxCliq of size" << maxCliq[i].at(0).size()<<endl;	
	#endif
    }
        
		
    infile.close();
    //cout<<"Done reading input data!"<<endl;
	
    #ifdef DEBUG
        cout<<"Display data:"<<endl;
	for(int i=0;i<NUMBER_TF;i++)
	    cout<<data[i]<<"\t";
	cout<<endl<<endl;;
	for(int i = 0;i<NUMBER_TF;i++)
	    cout<<data[ (NUMBER_GENE-1)*NUMBER_TF + i]<<"\t";
	cout<<endl;
	getchar();
    #endif
    
    //========================================================= 
    vector< vector<pair<int,int> > > invList;
    for(int eachS = 0; eachS < NUMBER_TF; eachS++)
    {
        vector<pair<int,int> > temp;
        for(int i = 0; i < NUMBER_GENE; i++)
        {
            for(int n = 0; n < maxCliq[i].size(); n++)
            {
                if (binary_search ((maxCliq[i]).at(n).begin(), (maxCliq[i]).at(n).end(), eachS))
                {
                    temp.push_back(make_pair (i,n));
                }
            }
        }
        invList.push_back(temp);
    }    
    
    // read the TF list for output ------------------ 
    int no_TF;
    if(option2) {
        FILE *fp;
        char TFname[20];
        fp = fopen(filenameTF,"r");
        if(fp == NULL) perror(filenameTF);
        else {  
            no_TF = 0;
            while( fscanf(fp, "%s", TFname) == 1){
                strcpy(TFlist[ no_TF ],TFname);
                no_TF++;  
            }
        }      
        fclose(fp);
        /*cout<<"Done reading "<< no_TF << " TFs from "<< filenameTF <<endl;
        cout<<"Check:: first one: "<< TFlist[0] <<endl;
        cout<<"Check:: last one: "<< TFlist[no_TF-1]<<endl<<endl; */
                      
        if(no_TF != NUMBER_TF){
             cout<<"Input file wrong! Check the no. of rows in "<< filenameTF<< " and no. of columns in "<< filename << "." <<endl;
             return -1;
        }
    }  
		
    for (int i = 0; i<NUMBER_GENE;i++) EnhancerP[i] = 1;
    // read the Enhancer Probability list ------------------ 
    float ep;
    int no_Gene;
    if(option3) {
        FILE *fp;
        fp = fopen(fileEnhancerP,"r");
        if(fp == NULL) {
	    perror(fileEnhancerP);
	    no_Gene=0;
		//return -1;
	}
        else {  
            no_Gene =0;
            while( fscanf(fp, "%f", &ep) == 1){
                EnhancerP[no_Gene] = ep;
		// cout<<no_Gene<<" , " <<ep<<endl;
                no_Gene++;  
            }
        }      
        fclose(fp);          
		   
	if(no_Gene != NUMBER_GENE){
            cout<<no_Gene<<", "<<NUMBER_GENE<<endl;
	    cout<<"Input Wrong! Check the no. of rows in "<< fileEnhancerP<<endl;
	    cout<<"Run without assigning probability for enhancers."<<endl<<endl;
	    option3 = 0;
            //return -1;
        }
    } 

    //------------------------------------------------- 
    for(int i=0; i<NUMBER_GENE;i++){
	maxCliq[i].clear();
    }
			
    switch( option1 ){
        case 1: {
	    vector< pair< vector<int>, pair< vector<int>, float> > > gsSet;
            gsSet = SampleGeneSearch3_2(NUMBER_TF,S,1,P,&invList,EnhancerP);
            break;
        }
	case 4: {
	    vector< pair< vector<int>, pair< vector<int>, float> > > gsSet;
  	    gsSet = SampleGeneSearch4(NUMBER_TF,2,M,P,&invList,EnhancerP); // find all significant patterns longer than 2 

	    for (int n=0; n<gsSet.size(); n++){
	        cout <<gsSet[n].first.size()<<"\t"<<gsSet[n].second.second<<"\t";
		for(int i=0; i<gsSet[n].second.first.size();i++){
		    if(i == gsSet[n].second.first.size() -1 ){
                        if(option2)
                            cout<<TFlist[gsSet[n].second.first.at(i)]<<endl;
                        else
                            cout<< gsSet[n].second.first.at(i)+1<<endl; // output TF ids start with 1
                    }
                    else{
                        if(option2)
                            cout << TFlist[gsSet[n].second.first.at(i)]<<"\t";
                        else
                            cout << gsSet[n].second.first.at(i)+1<<"\t";  // output TF ids start with 1                                       
                    }      
                }                
            }

            break;
        }

    } // end switch
  time (&end);
  dif = (double)(difftime (end,start)/ (double) 60);
  dif = (double)(difftime (end,start));
//  printf ("Elapsed time for generating GS sets: %.2lf seconds.\n", dif );  

}


