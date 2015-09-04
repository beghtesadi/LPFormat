#include <stdio.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>
#include <math.h>

/* scip includes */
#include "objscip/objscip.h"
#include "objscip/objscipdefplugins.h"

/* namespace usage */
using namespace std;
using namespace scip;

void GenerateDemand (vector<int> &x , vector<double> &y);
double CalProb (int n, double mean);
int CalMin (int  a, int b);
double mean=5;
double K=2;


int main ()
{
  vector<double> purchcost = {1,4,7};//{1,4,7};
  vector<int >vialsize ={1,5,10};
  double ConfRate=0.8;
  double Shortbeta=1;
  int I=3;
  double Wastbeta=1;
  double M=100;
  int Nper_to_spoil=1;
  int T=5; //T is the number of stages
  SCIP* scip = NULL;
  SCIP_CALL( SCIPcreate(&scip) );

  //FILE* fps = fopen("/home/beghtesadi/SCIPproject/lpformat/PSize.xls", "w");
  //SCIPinfoMessage(scip, fps, "NStages\tK\tNperspoil\tNxvar\tNyvar\tNzvar\tNsubcon\tNmascon\n")

   
  int NScen=pow(K, T-1);
  int Nnodes=0;
  for (int i=0; i<T ; i++)
  {
	Nnodes=Nnodes+pow(K,i);
  }
  vector < double > Dprob;
  vector < int > Demand;
  vector < double > Nodedemand;
  vector < double > Scenprob;
  vector < double > Scenprob2;
  vector < double > Nodeprobe;


  //GenerateDemand(Demand, Dprob);
  for (int k=0 ; k<K ; k++)
  {
	  Dprob.push_back(1/K);
	  //choose a random number between 1 and 10 for demand
	  //Demand.push_back( rand() % 10 + 1 );
  }


  //generating the nodes matrix
  vector< vector<int> > nodemat;
  vector< vector<int> > jmat;

	int n=1;
	for (int t=1 ; t<T+1 ; t++)
	{
		vector<int>row; //creat an empty vector
        vector<int>jrow;
		int len=NScen/pow(K,t-1);
		for (int i=0 ; i<pow(K,t-1) ; i++)
		{
			for (int j=0 ; j<len ; j++)
			{
				row.push_back(n);
                jrow.push_back(i);
			}
		n++;
		}
		nodemat.push_back(row);
        jmat.push_back(jrow);
	}

	for (int i=0 ; i<pow(K,T-1) ; i++)
	{
		Scenprob2.push_back(0);
	    Scenprob.push_back(0);
	}


	for (int j=0 ; j<K ; j++){
		//Nodedemand.push_back(Demand[j]);
		Nodedemand.push_back(rand() % 6 + 2);
		Nodeprobe.push_back(Dprob[j]);
		Scenprob[j]=Dprob[j];}

	for (int t=3; t<T+1 ; t++){
		int n=0;
		for (int i=0 ; i<pow(K,t-2) ; i++){
			for (int j=0 ; j < K ; j++){
				//Nodedemand.push_back(Demand[j]);
				Nodedemand.push_back(rand() % 6 + 2);
				Nodeprobe.push_back(Dprob[j]);
				Scenprob2[n]=Scenprob[i]* Dprob[j];
				n++;}}
		int ssize= (int)Scenprob2.size();
		for(int ii=0 ; ii<ssize ; ii++){
			Scenprob[ii]=Scenprob2[ii];}
	}

  
   /***********************
    * Version information *
    ***********************/

   SCIPprintVersion(scip, NULL);
   SCIPinfoMessage(scip, NULL, "\n");

   /* include default plugins */
   SCIP_CALL( SCIPincludeDefaultPlugins(scip) );

   /* set verbosity parameter */
   SCIP_CALL( SCIPsetIntParam(scip, "display/verblevel", 5) );
   SCIP_CALL( SCIPsetBoolParam(scip, "display/lpinfo", TRUE) ); 

   SCIP_CALL( SCIPsetSeparating(scip, SCIP_PARAMSETTING_DEFAULT, TRUE) );

   /* create empty problem */
   SCIP_CALL( SCIPcreateProb(scip, "vaccine", 0, 0, 0, 0, 0, 0, 0) );

  //vector < SCIP_Real > xub = {50 , 10 , 5};
  //Add x_{im} variables
  char var_name[255];
  vector< vector< vector< SCIP_VAR* > > > x_var(I);
  for (int i = 0; i < I; ++i)
   { 
    x_var[i].resize(T-1);
    for (int t=2; t<T+1 ; t++)
     {
      int len = pow(K,T-t);
      x_var[i][t-2].resize(pow(K,t-1)); 
      for (int j=0; j<pow(K,t-1) ; j++)
      {
        SCIP_VAR* xvar;
        SCIPsnprintf(var_name, 255, "x%d_%d", i , nodemat[t-1][j*len] );
        SCIP_CALL( SCIPcreateVar(scip,
                     &xvar,                                                      // returns new index
                     var_name,                                                   // name
                     0.0,                                                        // lower bound
                     SCIPinfinity(scip),                                                         // upper bound
                     Nodeprobe[nodemat[t-1][j*len]-2] *purchcost[i],             // objective
                     SCIP_VARTYPE_INTEGER,                                       // variable type
                     true,                                                       // initial
                     false,                                                      // forget the rest ...
                     0, 0, 0, 0, 0) );
        SCIP_CALL( SCIPaddVar(scip, xvar) );
        cout << "Nodeprobe" << Nodeprobe[nodemat[t-1][j*len]-2] << endl;
        cout << "purchcost" << purchcost[i] << endl;
        cout << Nodeprobe[nodemat[t-1][j*len]-2] *purchcost[i] << endl;
        x_var[i][t-2][j] = xvar;
      } 
    }
   } 
  /*Add y_{imn} variables*/
  vector< vector< vector< vector< vector< SCIP_VAR* > > > > > y_var (I);
  SCIP_Real ycoef =0;
  for (int i =0 ; i<I; i++)
  {
   y_var[i].resize(T-1);
   for (int t=2; t<T+1 ; t++)
   {
    int len = pow(K,T-t);
    y_var[i][t-2].resize(pow(K,t-1));
    for (int j=0 ; j<pow(K,t-1) ; j++)
    {
     int NNN= CalMin(Nper_to_spoil, T-t);     
     y_var[i][t-2][j].resize (NNN+1); 
     for (int kk=0; kk <NNN+1; kk++)
      {
       y_var[i][t-2][j][kk].resize(pow(K,kk));
       for(int jj=0; jj < pow(K,kk) ; jj++)
       { 
        SCIP_VAR* yvar;
        SCIPsnprintf(var_name, 255, "y%d_%d_%d_%d", i , nodemat[t-1][j*len] , kk , jj );
        SCIP_CALL( SCIPcreateVar(scip,
                     &yvar,                                                      // returns new index
                     var_name,                                                   // name
                     0.0,                                                        // lower bound
                     SCIPinfinity(scip),                                      // upper bound
                     ycoef,                                                   // objective
                     SCIP_VARTYPE_CONTINUOUS,                                       // variable type
                     true,                                                       // initial
                     false,                                                      // forget the rest ...
                     0, 0, 0, 0, 0) );
         SCIP_CALL( SCIPaddVar(scip, yvar) );
         y_var[i][t-2][j][kk][jj] = yvar;
         //SCIPinfoMessage(scip, NULL, "New Y variable is Y_%d_%d_%d_%d_%d", i, t, j, t+kk, j * (pow(_K, t+kk-1)/pow(_K, t-1))+jj);
       } // jj
      }//kk
     }//j
    }//t
   }   
  
  /*Add z_{m} variables */
  vector< vector< SCIP_VAR* > > z_var(T-1);
  for (int t=2 ; t<T+1; t++)
  {
     z_var[t-2].resize(pow(K,t-1));
     int len=pow(K,T-t);
     for( int j=0; j<pow(K,t-1) ; j++)
     {
       SCIP_VAR* var;
       SCIPsnprintf(var_name, 255, "Z%d", nodemat[t-1][j*len] );
       SCIP_CALL( SCIPcreateVar(scip,
                     &var,                   // returns new index
                     var_name,               // name
                     0.0,                    // lower bound
                     1.0,                    // upper bound
                     0,                      // objective
                     SCIP_VARTYPE_BINARY,   // variable type
                     true,                   // initial
                     false,                  // forget the rest ...
                     0, 0, 0, 0, 0) );
      SCIP_CALL( SCIPaddVar(scip, var) );
      z_var[t-2][j] = var;
    }
  }

  FILE* fp = fopen("/home/beghtesadi/SCIPproject/lpformat/original.dec", "w");
  SCIPinfoMessage(scip, fp, "PRESOLVED\n0\nNBLOCKS\n%d\n" , Nnodes-1);

   
  
 // SCIPinfoMessage(scip, fp, "0");

  /*Add constraint 2b, wastage constraint*/
  char con_name[255];
  vector< vector < vector <SCIP_CONS* > > > w_con(T-1);
  for (int t=2; t<T+1 ; t++)
   {
     int NN= CalMin(Nper_to_spoil,T-t);
     w_con[t-2].resize(pow(K,t-1));
     for (int j=0 ; j<pow(K,t-1); j++)
     {  
      w_con[t-2][j].resize(pow(K,NN));
      for (int n=0; n<pow(K,NN); n++)
      {
       SCIP_CONS* ww_con = NULL;
       SCIPsnprintf(con_name, 255, "w_con%d_%d_%d", t,j, n);
       SCIP_CALL( SCIPcreateConsLinear(scip, & ww_con, con_name, 0, NULL, NULL,
                     -SCIPinfinity(scip),    /* lhs */
                     Wastbeta,               /* rhs */
                     true,                   /* initial */
                     false,                  /* separate */
                     true,                   /* enforce */
                     true,                   /* check */
                     true,                   /* propagate */
                     false,                  /* local */
                     false,                   /* modifiable */
                     false,                  /* dynamic */
                     false,                  /* removable */
                     false) );               /* stickingatnode*/
      w_con[t-2][j][n]=ww_con;
      for (int ii=0; ii<I; ii++)
     {
       SCIP_CALL( SCIPaddCoefLinear(scip, w_con[t-2][j][n], x_var[ii][t-2][j], vialsize[ii]) );
       //SCIPinfoMessage(subscip, NULL, "X at this constraint X_%d", ii);
       for(int kk=0; kk < NN+1 ; kk++)
       {
         int leng = NScen/pow(K,t-1);
         int jjmat = jmat[t+kk-1][leng*j+n* NScen/pow(K,t+NN-1)];
         //assert(y_var[ii][t+kk-2][_jjmat] != NULL);
         SCIP_CALL( SCIPaddCoefLinear(scip, w_con[t-2][j][n], y_var[ii][t-2][j][kk][jjmat-j*pow(K,kk)], -1) );
         //SCIPinfoMessage(subscip, NULL, "Y at this constraint Y_%d_%d_%d", ii, t+kk, _jjmat);
        }
     }//ii
    assert (w_con[t-2][j][n] != NULL);
    SCIP_CALL( SCIPaddCoefLinear(scip, w_con[t-2][j][n], z_var[t-2][j], -1000) );
    SCIPinfoMessage(scip, NULL, "z at this constraint Z" );
    SCIP_CALL( SCIPaddCons(scip, w_con[t-2][j][n]) );
    }
   } 
   }
  /*Add constraint (2c) */ 
  vector< vector < vector < vector < SCIP_CONS* > > > > sec_con(I);
  for (int i=0; i<I; i++)
  {
   sec_con[i].resize(T-1);
   for (int t=2; t<T+1 ; t++)
   {
    int NN = CalMin (T-t,Nper_to_spoil);
    sec_con[i][t-2].resize(pow(K,t-1));
    for (int j=0 ; j<pow(K,t-1); j++)
    {
     sec_con[i][t-2][j].resize(pow(K,NN));
      for (int n=0; n<pow(K,NN); n++)
      {
        SCIP_CONS* s_con = NULL;
        SCIPsnprintf(con_name, 255, "sec_con%d_%d_%d_%d", i ,t,j, n);
        SCIP_CALL( SCIPcreateConsLinear(scip, & s_con, con_name, 0, NULL, NULL,
                     0,                      
                     SCIPinfinity(scip),    
                     true,                   
                     false,                  
                     true,                  
                     true,                   
                     true,                  
                     false,                 
                     false,
                     false,                  
                     false,                  
                     false) );               
       sec_con[i][t-2][j][n] = s_con;
       SCIP_CALL( SCIPaddCoefLinear(scip, sec_con[i][t-2][j][n], x_var[i][t-2][j], vialsize[i]) ); 
       for(int kk=0; kk < NN+1 ; kk++)
        {
          int leng = NScen/pow(K,t-1);
          int jjmat = jmat[t+kk-1][leng*j+n* NScen/pow(K,t+NN-1)];
          //assert (y_var[ii][t+kk-2][_jjmat] != NULL);
          SCIP_CALL( SCIPaddCoefLinear(scip, sec_con[i][t-2][j][n], y_var[i][t-2][j][kk][jjmat-j*pow(K,kk)], -1) );
        }
      SCIP_CALL( SCIPaddCons(scip,sec_con[i][t-2][j][n]) );
      }
    }
   }
  }
   /* Print the decomposition structure*/
   for (int t=2 ; t<T+1 ; t++)
   {
    int NN= CalMin(Nper_to_spoil,T-t);
    int len = pow(K,T-t);
    for (int j=0 ; j<pow(K,t-1); j++)
    {
     SCIPinfoMessage(scip , fp, "BLOCK\t%d\n", nodemat[t-1][j*len]-1);
     for (int n=0; n<pow(K,NN); n++)
      {
        SCIPinfoMessage(scip , fp, "w_con%d_%d_%d\n", t ,j, n);
       }
     for (int i=0 ; i<I ; i++)
      {
        for (int n=0; n<pow(K,NN); n++)
         {
           SCIPinfoMessage(scip , fp, "sec_con%d_%d_%d_%d\n", i, t ,j, n);
          }
       }
     }
    }
  SCIPinfoMessage(scip , fp, "MASTERCONSS\n");
  
  /*Add shortage constraint */
  vector< vector < vector <SCIP_CONS* > > > sh_con(T-1);
  SCIP_Real mwj = M ;
  for(int t=2; t<T+1; t++)
  {
    sh_con[t-2].resize(pow(K,t-1));
    int lenght=NScen/pow(K,t-1); 
    for (int j=0 ; j<pow(K,t-1) ; j++)
      {
         sh_con[t-2][j].resize(pow(K,T-t));
         for (int q=0 ; q<pow(K,T-t) ; q++)
            {
                SCIP_CONS* con;
                SCIPsnprintf(con_name, 255, "sh_con%d_%d_%d", t, j, q);
                SCIP_CALL( SCIPcreateConsLinear(scip, &con, con_name, 0, NULL, NULL,
                     Nodedemand[nodemat[t-1][j*lenght]-2]-Shortbeta,    /* lhs */
                     SCIPinfinity(scip),     /* rhs */
                     true,                   /* initial */
                     false,                  /* separate */
                     true,                   /* enforce */
                     true,                   /* check */
                     true,                   /* propagate */
                     false,                  /* local */
                     false,                   /* modifiable */
                     false,                  /* dynamic */
                     false,                  /* removable */
                     false) );               /* stickingatnode */
                sh_con[t-2][j][q] = con;
                SCIPinfoMessage(scip , fp, "sh_con%d_%d_%d\n", t, j, q);
                SCIP_CALL( SCIPaddCoefLinear(scip, sh_con[t-2][j][q], z_var[T-2][j * pow(K,T-t) + q], Nodedemand[nodemat[t-1][j*lenght]-2]) );
                for (int i=0 ;i<I ; i++)
                {
                 int NNN = CalMin (Nper_to_spoil ,t-2);
                 for (int kk=0 ; kk < NNN + 1 ; kk++)
                 {
                   int lenn = NScen/pow(K, t-1);
                   int jjmat = jmat[t-kk-1][lenn*j];
                   SCIP_CALL( SCIPaddCoefLinear(scip, sh_con[t-2][j][q], y_var[i][t-kk-2][jjmat][kk][j - jjmat*pow(K,kk)], 1) );
                   }//kk
                 }//i
                SCIP_CALL( SCIPaddCons(scip, sh_con[t-2][j][q]) );
             }
       }
   } //Done to here

   /*Add constraint (2f) */
   SCIP_Real zcoef = -1;
   vector<vector<SCIP_CONS*>>z_con (NScen);  
   for (int i=0; i<NScen; i++)
   {  
      z_con[i].resize(T-1);
      for (int t=2; t<T+1 ; t++)
       {
          SCIP_CONS* con;
          SCIPsnprintf(con_name, 255, "z_con%d_%d", i, t);
          SCIP_VAR* index = z_var[T-2][i];
          SCIP_Real coeff = 1;
          SCIP_CALL( SCIPcreateConsLinear(scip, &con, con_name, 1, &index, &coeff,
                     0,                      /* lhs */
                     SCIPinfinity(scip),     /* rhs */
                     true,                   /* initial */
                     false,                  /* separate */
                     true,                   /* enforce */
                     true,                   /* check */
                     true,                   /* propagate */
                     false,                  /* local */
                     false,                   /* modifiable */
                     false,                  /* dynamic */
                     false,                  /* removable */
                     false) );               /* stickingatnode */
                SCIP_CALL( SCIPaddCons(scip, con) );
                z_con[i][t-2] = con;
                SCIPinfoMessage(scip , fp, "z_con%d_%d\n", i, t);
                int jjjmat= jmat [t-1][i];
                SCIP_CALL( SCIPaddCoefLinear(scip, z_con[i][t-2], z_var[t-2][jjjmat], zcoef) );
        }
    }

   /*Add constraint 2e */
  SCIP_CONS*  alpha_con = NULL;
   SCIP_CALL( SCIPcreateConsLinear(scip, & alpha_con,"AlphaCon" , 0, NULL, NULL,
                     -SCIPinfinity(scip),    /* lhs */
                     1-ConfRate,             /* rhs */
                     true,                   /* initial */
                     false,                  /* separate */
                     true,                   /* enforce */
                     true,                   /* check */
                     true,                   /* propagate */
                     false,                  /* local */
                     false,                   /* modifiable */
                     false,                  /* dynamic */
                     false,                  /* removable */
                     false) );               /* stickingatnode */

   for (int i = 0; i < NScen; ++i)
	 SCIP_CALL( SCIPaddCoefLinear(scip, alpha_con, z_var[T-2][i], Scenprob[i]) );
    SCIP_CALL( SCIPaddCons(scip, alpha_con) );
    SCIPinfoMessage(scip , fp, "AlphaCon");

   //SCIPinfoMessage(scip, fps, "T\tK\tNper_to_spoil\tNxvar\tNyvar\tNBinary\tNsubcon\tNmascon\n");}}}

         
   SCIP_CALL( SCIPwriteOrigProblem(scip, "original.lp", "lp", FALSE) );


   /*************
    *  Solve    *
    *************/

  SCIP_CALL( SCIPsolve(scip) );


   /**************
    * Statistics *
    *************/
  SCIP_CALL( SCIPprintStatistics(scip, NULL) );

   SCIP_CALL( SCIPprintBestSol(scip, NULL, FALSE) );



   /********************
    * Deinitialization *
    ********************/

   SCIP_CALL( SCIPfree(&scip) );

   //BMScheckEmptyMemory();

   return 0;
}//for main 

void GenerateDemand (vector<int> &x , vector<double> &y){
        int a, b;
        a = CalMin(mean,ceil((K-1)/2));
        if (mean > ceil ((K-1)/2))
           b= floor((K-1)/2);
        else
           b= K-mean-1;
        double prob = 0;
        for (int i= 0; i< mean-a +1 ; i++){
          prob = prob +CalProb (i , mean);
          }
        double probb = prob;
        x.push_back(mean-a);
        y.push_back(prob);
        //y.push_back(CalProb(mean-a,mean));
	for (int i=1 ; i<a+b ; i++){
           y.push_back(CalProb(mean-a+i , mean)); 
	   x.push_back(mean-a+i);
           probb = probb +y[i];}
        x.push_back (mean+b);
        y.push_back(1-probb); 
        //y.push_back(1-probb); 
        
}

double CalProb (int n, double mean) 
 {
  double prob=exp(-1*mean);
  if (n > 0)
  {
    for (int i=1 ; i<n+1 ; i++)
    {
     prob = prob * mean / i ;
    }
   } 
  return (prob);
 }

int CalMin (int  a, int b)
{
   int min;
   if (a < b)
     min=a;
   else
     min=b;
   return (min);
}

