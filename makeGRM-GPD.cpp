#include <string.h>
#include <math.h>
#include <cstdio>
#include <string>
#include <iostream>
#include <sstream>
#include <fstream>
#include <cstdlib>

#define PACK_DENSITY 4
#define MASK0   3 // 3 << 2 * 0
#define MASK1  12 // 3 << 2 * 1
#define MASK2  48 // 3 << 2 * 2
#define MASK3 192 // 3 << 2 * 3

using namespace std;
#include "Eigen/Dense"
using namespace Eigen;

void decode_plink(char *output, const char *input, const int lengthInput){
  int i, k;
  char tmp, geno;
  int a1, a2;
  
  for(i=0;i<lengthInput;++i){
    tmp = input[i];
    k   = PACK_DENSITY * i;
    geno      = (tmp & MASK0);
    a1        = !(geno & 1);
    a2        = !(geno >> 1);
    output[k] = (geno == 1) ? 3 : a1 + a2;
    k++;
    
    geno      = (tmp & MASK1) >> 2; 
    a1        = !(geno & 1);
    a2        = !(geno >> 1);
    output[k] = (geno == 1) ? 3 : a1 + a2;
    k++;
    
    geno      = (tmp & MASK2) >> 4; 
    a1        = !(geno & 1);
    a2        = !(geno >> 1);
    output[k] = (geno == 1) ? 3 : a1 + a2;
    k++;
    
    geno      = (tmp & MASK3) >> 6; 
    a1        = !(geno & 1);
    a2        = !(geno >> 1);
    output[k] = (geno == 1) ? 3 : a1 + a2;
  }
}

// calcGRM(z,zz_k,G,W,B,bedfile,n,p,nChr,packed,unpacked,numBytes,j_start,j_stop);
void calcGRM( MatrixXf &z, MatrixXf &zz_k, MatrixXf &G, MatrixXf &W, MatrixXf &B,
              string bedfile, int n, int p, int nChr, char *packed, char *unpacked, 
              int numBytes, int *j_start, int *j_stop, int *nchrs){
  ifstream influx;
  influx.open(bedfile.c_str(), std::ios::in | std::ios::binary);
  if(!influx){
    cerr << "[readGenotypes] Error reading file "<<bedfile<<endl;
    exit(1);
  }
  influx.seekg(0, ifstream::end);
  influx.seekg(3, ifstream::beg);
  
  int i,j,k;
  int nEff;
  float x;
  float dp_j, h_j, s_j;

  float pi_k;
  float C = 0.;
    
  int Meff = 0;
  for(int iChr=0;iChr<nChr;iChr++){
    int chr  = 1 + iChr;
    int jsp1 = j_start[iChr]+1;
    cout<<"# => Processing chromosome "<<chr<<" (SNP_"<<jsp1<<" to SNP_"<<j_stop[iChr]<<"). <= \n";
    k = 0;
    for(j=j_start[iChr];j<j_stop[iChr];j++){
      influx.read((char*)packed, sizeof(char) * numBytes);
      decode_plink(unpacked, packed, numBytes);
      dp_j = 0.0;
      nEff = 0;
      
      for(i=0;i<n;i++){
        x      = (float)((int) unpacked[i]);
        z(i,k) = x;
        if(x!=3.0){
          nEff++;
          dp_j += x;
        }
      }
      dp_j = dp_j / nEff;
      h_j  = dp_j * (1 - 0.5*dp_j);
      s_j  = sqrt(h_j);
      
      if(s_j>0.){
        Meff++;
        for(i=0;i<n;i++){
          if(z(i,k)==3.){
            z(i,k) = 0.; // imputed to the mean
          }else{
            z(i,k) = (z(i,k) - dp_j) / s_j;
          }
        }
        k++;
      }else{
        cout<<"*** SNP"<<j<<" is monomorphic! - please remove such SNPs.\n";
        exit(1);
      }
    }// j
    
    if(k!=nchrs[iChr]) cout<<"Mismatch!\n";
    pi_k = ((float) nchrs[iChr]) / p; 
    C   += pi_k * pi_k;
    
    // cout<<"\t\t*** Last value of k = "<<k<<endl;
    zz_k.noalias()   = z.block(0,0,n,k) * z.block(0,0,n,k).transpose();
    G.noalias()     += zz_k;
    W.noalias()     += zz_k * zz_k;
    z.setZero();
  }// iChr
  
  cout<<"# C constant is "<<C<<endl;
  cout<<"# GRM calculated over "<<Meff<<" / "<<p<<" SNPs.\n";
  float Mf  = (float) Meff;
  float nf  = (float) n;
  float nm2 = nf * Mf * Mf;
  
  B.noalias()  = ( G * G - W );
  W.noalias() -= nf * G;
  G.noalias()  = (1. / Mf) * G;
  
  B.noalias()  = (1./nm2) * B;
  W.noalias()  = (1./nm2) * W;
  
  // Add normaliszing constants + mean trace etc....
  float G_bar     = G.trace() / n;
  float W_bar     = W.trace() / n;
  float B_bar     = B.trace() / n;
  
  cout<<"# Tr(G)/n = "<<G_bar<<" - Tr(W)/n = "<<W_bar<<" - Tr(B)/n = "<<B_bar<<".\n";
  cout<<"#\n";
  
  B.noalias()  = (1. - C) * (1./B_bar) * B;
  W.noalias()  = C * (1./W_bar) * W;
  
  influx.close();
}

// Main funtion used on cov_files.txt with a sequential access
int main(int argc, char *argv[]){
  // Input arguments
  string bfile       = "none";
  string grmPrefix   = "none";
  bool verbose       =  true;
  // Indices
  string sw;
  int i,j;

  if(argc==1){
    cerr<<"\tArguments must be specified. Type --help for more details."<<endl;
    exit(1);
  }
  
  // Read arguments
  sw = argv[1];
  if (sw == "--help"){
    cerr<<"\t--bfile      : Binary PLINK format for genotypes."<<endl;
    cerr<<"\t--out        : Prefix for output file: [prefix].grm.bin. Default is [none]."<<endl;
    exit(1);
  }else{
    if (argc == 1) {
      cerr<<"\tArguments must be specified. Type --help for more details."<<endl;
      exit(1);
    }
  }
  
  for(i = 1; i<argc;i++){
    sw = argv[i];
    if (sw == "--bfile"){
      bfile = argv[i + 1];
    }
    if (sw == "--out"){
      grmPrefix = argv[i + 1];
    }
  }
 

  if(bfile=="none"){
    cerr<<"\tA prefix must be specified for files [prefix].bed, [prefix].bim and [prefix].fam"<<endl;
    exit(1);
  }

  string bedfile = bfile+".bed";
  string bimfile = bfile+".bim";
  string famfile = bfile+".fam";
  
  initParallel();
  int nThreads = nbThreads( );

  clock_t tic = clock();
  time_t t = time(0);   // get time now
  struct tm * now = localtime( & t );
  
  if(verbose){
    cout <<"# >>> CALCULATING GRMs <<<"<<endl;
    cout <<"# Analysis starts : ";
    cout << (now->tm_year + 1900) << '-'
         << (now->tm_mon + 1) << '-'
         <<  now->tm_mday << " at "
         <<  now->tm_hour <<":"<<now->tm_min<<":"<<now->tm_sec
         <<  ".\n";
    cout<<"# Input files: "<<endl;
    cout<<"# BED: "<<bfile<<".bed.\n";
    cout<<"# BIM: "<<bfile<<".bim.\n";
    cout<<"# FAM: "<<bfile<<".fam.\n";
    cout<<"# Output file names: "<<grmPrefix<<".grm.[bin|id]\n";
    cout<<"# Using "<<nThreads<<" threads.\n";
  }
  
  // Few tools
  string line = "";
  string tok  = ""; 
  ifstream tmpStream;
  
  // Get number of SNPs
  if(verbose){
    //cout<<"## 1-Counting the number of SNPs..."<<endl;
  }  
  int p = -1;
  tmpStream.open(bimfile.c_str());
  while(tmpStream){
    getline(tmpStream,line);
    p++;
  }
  tmpStream.close();
    
  int nChr  = 0;
  int *chroms = new int[p];
  tmpStream.open(bimfile.c_str());
  for(j=0;j<p;j++){
    getline(tmpStream,line);
    stringstream ss;
    ss << line;
    ss >> tok;
    chroms[j] = atoi(tok.c_str());
    if(chroms[j]>nChr) nChr = chroms[j];
  }
  tmpStream.close();
  
  int *nchrs   = new int[nChr];
  int *j_start = new int[nChr]; //  j >= j_start
  int *j_stop  = new int[nChr]; //  j  < j_stop
  
  int k;
  for(k=0;k<nChr;k++){
    nchrs[k] = 0;
    j_start[k] = 0;
    j_stop[k] = 0;
  }
  
  j_start[0]     = 0;
  j_stop[nChr-1] = p;
  
  int jp1 = 0;
  int kp1 = 0;
  for(j=0;j<p;j++){
    k = chroms[j]-1;
    kp1 = k + 1;
    jp1 = j + 1;
    nchrs[k]++;
    if(jp1<p){
      if(chroms[jp1] > chroms[j]){
        j_start[kp1] = jp1;
        j_stop[k]    = jp1;
        //cout<<"\tTransition from chromosome "<<chroms[j]<<" to chromosome "<<chroms[jp1]<<" happens for j="<<j<<".\n";
      }
    }
  }
  
  int m = 0;
  for(k=0;k<nChr;k++){
    if(nchrs[k]>m) m = nchrs[k];
  }
  if(verbose){
    cout<<"# Found "<<p<<" SNPs across "<<nChr<<" chromosomes."<<endl;
    cout<<"# Maximum block size (max number of SNPs on a given chromosome) = "<<m<<" SNPs."<<endl;
  }
  
  // Get sample size
  if(verbose){
    //cout<<"2-Counting the number of samples..."<<endl;
  }    
  int n = -1; 
  tmpStream.open(famfile.c_str());
  while(tmpStream){
    getline(tmpStream,line);
    n++;
  }
  tmpStream.close();  
  if(verbose){
    cout<<"# Found "<<n<<" samples."<<endl;
  }
  
  //exit(1);
  
  string *FID = new string[n];
  string *IID = new string[n];
  tmpStream.open(famfile.c_str());
  for(i=0;i<n;i++){
    getline(tmpStream,line);
    stringstream ss;
    ss << line;
    ss >> FID[i];
    ss >> IID[i];
  }
  tmpStream.close();
  
  
  
  int numBytes   = (int)ceil((double)n / PACK_DENSITY);
  char* packed   = new char[numBytes];
  char* unpacked = new char[numBytes * PACK_DENSITY];
  
  MatrixXf G    = MatrixXf::Zero(n,n); // standard GRM
  MatrixXf W    = MatrixXf::Zero(n,n); // G_within - W matrix
  MatrixXf B    = MatrixXf::Zero(n,n); // G_between - B matrix
  MatrixXf zz_k = MatrixXf::Zero(n,n); // G_k = zz_k / M_k;
  MatrixXf    z = MatrixXf::Zero(n,m);

  calcGRM(z,zz_k,G,W,B,bedfile,n,p,nChr,packed,unpacked,numBytes,j_start,j_stop,nchrs);
  
  z.resize(0,0); // free - memory
  int ncolX = 1; // Just intercept
  int ncolA = 1 + ncolX + 3 * n;
 
  // Output grm
  string grmOutIdFile;
  string grmOutBinFile;
  
  for(int ip=0;ip<3;ip++){
    if(ip==0){
      grmOutIdFile  = grmPrefix+".G.grm.id";
      grmOutBinFile = grmPrefix+".G.grm.bin";
      cout<<"# Writing standard GRM...\n";
    }
    if(ip==1){
      grmOutIdFile  = grmPrefix+".W.grm.id";
      grmOutBinFile = grmPrefix+".W.grm.bin";
      cout<<"# Writing within-chromosome GRM...\n";
    }
    if(ip==2){
      grmOutIdFile  = grmPrefix+".B.grm.id";
      grmOutBinFile = grmPrefix+".B.grm.bin";
      cout<<"# Writing between-chromosome GRM...\n";
    }
    
    ofstream fileId(grmOutIdFile.c_str());
    for(i=0;i<n;i++){
      fileId<<FID[i]<<"\t"<<IID[i]<<endl;
    }
    fileId.close();
    
    float f_buf = 0.0;
    int size = sizeof (float);
    
    float sGd  = 0.;
    float sG2d = 0.;
    
    float sGo  = 0.;
    float sG2o = 0.;
    
    int npd = 0;
    int npo = 0;
 
    fstream A_Bin(grmOutBinFile.c_str(), ios::out | ios::binary);
    if (!A_Bin) throw (("Error: can not open the file [" + grmOutBinFile + "] to write.").c_str());
    for (i = 0; i < n; i++) {
      for (j = 0; j <= i; j++) {
        if(ip==0) f_buf = G(i, j);
        if(ip==1) f_buf = W(i, j);
        if(ip==2) f_buf = B(i, j);
        if(i==j){
          sGd  += f_buf;
          sG2d += f_buf * f_buf;
          npd++;
        }else{
          sGo  += f_buf;
          sG2o += f_buf * f_buf;
          npo++;
        }
        A_Bin.write((char*) &f_buf, size);
      }
    }
    A_Bin.close();
    cout << "# GRM of " << n << " individuals has been saved in the file [" + grmOutBinFile + "] (in binary format)." << endl;
    
    if(npd != n) cout<<"# **** we might have a problem here!\n";
    
    float meanGd = sGd / npd;
    float meanGo = sGo / npo;
    float varGd  = sG2d / (npd-1) - meanGd * meanGd * (1. + 1. / (npd - 1));
    float varGo  = sG2o / (npo-1) - meanGo * meanGo * (1. + 1. / (npo - 1));
    
    cout<<"# mean(G_diag) = "<<meanGd<<" - var(G_diag) = "<<varGd<<endl;
    cout<<"# mean(G_Offdiag) = "<<meanGo<<" - var(G_Offdiag) = "<<varGo<<" (over "<<npo<<" pairs of individuals)."<<endl;
    cout << "#\n";
  }
  
  time_t t2 = time(0);   // get time now
  struct tm * now2 = localtime( & t2 );
  cout <<"# Analysis ends: ";
  cout << (now2->tm_year + 1900) << '-'
       << (now2->tm_mon + 1) << '-'
       <<  now2->tm_mday << " at "
       <<  now2->tm_hour <<":"<<now2->tm_min<<":"<<now2->tm_sec
       << ".\n";
  clock_t toc = clock();
  float time_elapsed = (float)(toc - tic) / CLOCKS_PER_SEC;
  printf("# Time elapsed: %f seconds.\n\n", time_elapsed / nThreads);
  
  
  delete [] packed;
  delete [] unpacked;    
  delete [] FID;
  delete [] IID;
  delete [] chroms;
  delete [] nchrs;
  delete [] j_start;
  delete [] j_stop;

  return EXIT_SUCCESS;
}


