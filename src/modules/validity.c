/*-------------------------------------------------------------------*/
/* Cluster Validity index        Qinpei Zhao                         */
/*                                                                   */
/* Different types of cluster validity indices                       */
/*                                                                   */    
/*                                                                   */
/* Naming conventions used in the code                               */
/*                                                                   */
/*    TS        training set (data objects)                          */
/*    CB        codebook (cluster representatives, centroids)        */
/*    P         partitioning (pointing from TS to CB)                */
/*                                                                   */
/*    p-prefix  pointer, e.g. pTS is pointer to the training set TS  */
/*                                                                   */
/*-------------------------------------------------------------------*/

#define ProgName       "ClusterValidity"
#define VersionNumber  "Version 0.0" /* qp */
#define LastUpdated    "08.03.2008"

/*-------------------------------------------------------------------*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "cb.h"
#include "math.h"
#include "validity.h"
#include "interfc.h"
#include "memctrl.h"



#define DBL_MAX   10000000000
#define pi  3.1415
#define max(a, b)  (((a) > (b)) ? (a) : (b))

/*-------------------------------------------------------------------
VECTORTYPE MeanOfTs(TRAININGSET *pTS, CODEBOOK *pCB);
double ValidityMSE(TRAININGSET *pTS, CODEBOOK *pCB, PARTITIONING *pP);
double ValiditySSB(TRAININGSET *pTS, CODEBOOK *pCB, PARTITIONING *pP);
double ValidityDunn(TRAININGSET *pTS, CODEBOOK *pCB, PARTITIONING *pP);
double ValidityDBI(TRAININGSET *pTS, CODEBOOK *pCB, PARTITIONING *pP);
double ValidityFratio(TRAININGSET *pTS, CODEBOOK *pCB, PARTITIONING *pP);
//another implementation of F-ratio i.e. f-test from ISMO
double FTest(TRAININGSET *pTS, CODEBOOK *pCB, PARTITIONING *pP);
double Silhouette(TRAININGSET *pTS, CODEBOOK *pCB, PARTITIONING *pP);
double SilhouetteCI(TRAININGSET *pTS, CODEBOOK *pCB, PARTITIONING *pP);
double Xie_Beni(TRAININGSET *pTS, CODEBOOK *pCB, PARTITIONING *pP);
double Xie_Beni_Fuzzy(TRAININGSET *pTS, CODEBOOK *pCB, PARTITIONING *pP);
double ValidityBIC(TRAININGSET *pTS, CODEBOOK *pCB, PARTITIONING *pP);
double RandIndex(TRAININGSET *pTS, PARTITIONING *pP, PARTITIONING *pPtruth);
double CorrectedRI(TRAININGSET *pTS, PARTITIONING *pP, PARTITIONING *pPtruth);
double NMI(TRAININGSET *pTS, CODEBOOK *pCB, PARTITIONING *pP, PARTITIONING *pPtruth);

double S_Dbw(TRAININGSET *pTS, CODEBOOK *pCB, PARTITIONING *pP);
-------------------------------------------------------------------*/



/* =================== PUBLIC FUNCTIONS ============================= */
/*-------------------------------------------------------------------*/
// function for getting mean vector of total dataset
VECTORTYPE MeanOfTs(TRAININGSET *pTS, CODEBOOK *pCB)
{
	int          i = 0;
	int          k = 0;
	double       vector_sum = 0;
	double       vector_mean = 0;
	VECTORTYPE   MeanVector;
	MeanVector = CreateEmptyVector(VectorSize(pCB));
	// mean vector of all training set/
	// VectorSize => dimension
	for( k = 0; k < VectorSize(pCB); k++ )
	{
		for( i = 0; i < BookSize(pTS); i++)	
		{
			vector_sum += (llong) VectorScalar(pTS,i,k);
		}
		vector_mean = vector_sum/BookSize(pTS);
		MeanVector[k] = vector_mean;
	}
	return MeanVector;
}
// end of MeanOfTs


// MSE:
double ValidityMSE(TRAININGSET *pTS, CODEBOOK *pCB, PARTITIONING *pP)
{
	int		   i = 0;
	int        j = 0;
	double     sum = 0;
	double     MSE = 0;
	double     dist = 0;
	for( i = 0; i < BookSize(pTS); i++)
	{
		j = Map(pP,i);
		dist = VectorDistance(Vector(pTS, i), Vector(pCB, j), 
			VectorSize(pTS), MAXLLONG, EUCLIDEANSQ);
		sum += dist * VectorFreq(pTS, i);	
		//fprintf(stdout, "VectorFreq=%f\n", VectorFreq(pTS, i));
	}
	// double SSW = sum;
	MSE = sum /(TotalFreq(pTS) * VectorSize(pTS));
	return MSE;
}
//end of MSE

//SSB: sum-of-squares of variance between clusters
double ValiditySSB(TRAININGSET *pTS, CODEBOOK *pCB, PARTITIONING *pP)
{
	int        i = 0;
	int        j = 0;
	double     SSB = 0;

	for (i = 0; i < BookSize(pCB); i++){
		for(j = i+1; j < BookSize(pCB); j++)
		{
			double temp = VectorDistance(Vector(pCB, i), Vector(pCB, j), 
			VectorSize(pCB), MAXLLONG, EUCLIDEANSQ);
			SSB = SSB + temp;
		}
	}
	//average SSB
	//SSB = SSB*2/(M*(M-1));
	return SSB;

}
//end SSB

double ValidityDunn(TRAININGSET *pTS, CODEBOOK *pCB, PARTITIONING *pP)
{
	// min distance between clusters	
	double		min_dist = DBL_MAX;
	double		max_dist = 0;
	double		m_max = 0; //max in k clusters
	double		dist = 0;
	int         i = 0;
	int         j = 0;
	int         k = 0;
	int         t = 0;
	
	for( i = 0; i < PartitionCount(pP); i++)
		for(j = i+1; j < PartitionCount(pP); j++)
			for( k = FirstVector(pP, i); !EndOfPartition(k); k = NextVector(pP, k))
				for( t = FirstVector(pP, j); !EndOfPartition(t); t = NextVector(pP, t))
				{
					dist = VectorDistance(Vector(pTS, k), Vector(pTS, t),
									VectorSize(pTS), MAXLLONG, EUCLIDEANSQ);
					if (dist < min_dist)
						min_dist = dist;			
				}

		
	// max distance in clusters
	for( k = 0; k < PartitionCount(pP); k++)
	{
		for( i = 0; i < BookSize(pTS); i++)
			for( j = i+1; j < BookSize(pTS); j++)
			{
				if(Map(pP,i)== k && Map(pP, j) == k)
					dist = VectorDistance(Vector(pTS, i), Vector(pTS, j),
						VectorSize(pTS), MAXLLONG, EUCLIDEANSQ);
				if(dist > max_dist)
					max_dist = dist;
			}
		if(max_dist > m_max)
			m_max = max_dist;
	}
	return min_dist/m_max;

}/*ValidityDunn( )*/

// DBI implementaion:
/* 2k7-12-7
DBI: the ratio of the sum of within-cluster scatter to between-cluster
separation. the lower the better
*/
double ValidityDBI(TRAININGSET *pTS, CODEBOOK *pCB, PARTITIONING *pP)
{
	int		   i = 0;
	int        j = 0;
	int        t = 0;
	int        NumEachCluster = 0;
	llong      sum = 0;
	llong      DistClusters;
	double     DBIndex = 0;
	double     SumRatio = 0;
	double     CurrRatio = 0;
	llong      Scatter[PartitionCount(pP)];	
		
	//get the Si,Sj - within cluster scatter
	for( i = 0; i < PartitionCount(pP); i++)
	{ 
		NumEachCluster = 0;
		sum = 0;
		for( t = 0; t < BookSize(pTS); t++)
		{
			if(Map(pP,t) == i)
			{	
				NumEachCluster++;
				llong temp = sqrt(VectorDistance(Vector(pTS, t), Vector(pCB, i),
					VectorSize(pTS), MAXLLONG, EUCLIDEANSQ));
				sum += temp*temp;
			}
		}
		Scatter[i] = sum/NumEachCluster;
	}
	//get the ratio here	
	for(i = 0; i < PartitionCount(pP); i++)
	{
		double MaxRatio = 0;
		llong  temp = 0;
		for(j = i+1; j < PartitionCount(pP); j++)
		{
			CurrRatio = 0;
			DistClusters = 0;
			//distance between 2 centroids
			DistClusters = sqrt(VectorDistance(Vector(pCB, i), Vector(pCB, j),
					VectorSize(pCB), MAXLLONG, EUCLIDEANSQ));
			DistClusters *= DistClusters;
			//ratio of the sum of within cluster scatter to between cluster separation
			temp = Scatter[i]+Scatter[j];
			
			CurrRatio = (double)temp/(double)DistClusters;
			if(CurrRatio > MaxRatio)
				MaxRatio = CurrRatio;
		}
		SumRatio += MaxRatio;
	}
	DBIndex = SumRatio/PartitionCount(pP);	
	return DBIndex;	
	
}
//end of DBI



// f-ratio F-ratio is defined as the ratio of the within-class variance 
// against the between-class variance. the smaller, the more seperated 
// the clusters are.
double ValidityFratio(TRAININGSET *pTS, CODEBOOK *pCB, PARTITIONING *pP)
{
	int		   i = 0;
	int        j = 0;
	int        t = 0;
	double     InterGroup_sum = 0;
	double     temp_dist = 0;
	double     F_ratio = 0;
	double	   sum = 0;
	double     IntraGroup = 0;
	VECTORTYPE MeanVector;
	// mean vector of all training set/
	MeanVector = MeanOfTs(pTS,pCB);
	// calculate between-group variance
	for( j = 0; j < PartitionCount(pP); j++)
	{
		int cluster_size = 0;
		for( t = FirstVector(pP, j); !EndOfPartition(t); t = NextVector(pP, t))
		{
			cluster_size ++;
		} // actually cluster_size = VectorFreq(pCB, i);
		temp_dist = sqrt(VectorDistance(Vector(pCB, j), MeanVector,
			VectorSize(pCB), MAXLLONG, EUCLIDEANSQ));
		temp_dist *= temp_dist;
		temp_dist /= PartitionCount(pP);
		InterGroup_sum += cluster_size * temp_dist;				
	}
//	InterGroup_sum /= (VectorSize(pCB)-1);
	// calculate within-group variance
	for( i = 0; i < BookSize(pTS); i++)
	{
		j = Map(pP, i);
		temp_dist = sqrt(VectorDistance(Vector(pTS, i), Vector(pCB, j), 
			VectorSize(pTS), MAXLLONG, EUCLIDEANSQ) * VectorFreq(pTS, i));
		sum += temp_dist*temp_dist;
	}
	IntraGroup = sum; // / (BookSize(pTS)-PartitionCount(pP))
//	IntraGroup /= (BookSize(pTS)-VectorSize(pCB));
//	InterGroup_sum = InterGroup_sum / PartitionCount(pP);
	F_ratio =  IntraGroup/InterGroup_sum ;
	return F_ratio;
}
//end of F_ratio

//F-test: refer to implementation in ftest.c 
//function float CalculateFTest(DistanceInfo* DI, TRAININGSET* TS,
//    CODEBOOK* CB, PARTITIONING* PA, void* AuxiliaryData)
double FTest(TRAININGSET *pTS, CODEBOOK *pCB, PARTITIONING *pP)
{
	int          j = 0;
	double       QB = 0; //distance between clusters
	double       NMSE = 0;
	double       MSE = 0;
	double       temp_dist = 0;
	VECTORTYPE   MeanVector;
	// mean vector of all training set/
	MeanVector = MeanOfTs(pTS,pCB);
	// calculate between-group variance		
	for( j = 0; j < PartitionCount(pP); j++)
	{		
		// EUCLIDEANSQ (without square)
		temp_dist = sqrt(VectorDistance(Vector(pCB, j), MeanVector,
			VectorSize(pCB), MAXLLONG, EUCLIDEANSQ));
		temp_dist *= temp_dist;	
		QB += VectorFreq(pCB,j)*temp_dist;
	}
	// Inter-Group
	QB /= VectorSize(pTS) * (BookSize(pCB) - 1);
	
	//Intra-Group
	//calculte MSE distance intra-cluters	
	MSE = ValidityMSE(pTS,pCB,pP); 
	// end MSE
	NMSE = TotalFreq(pTS)* MSE;		
	NMSE /= BookSize(pTS) - BookSize(pCB);
	return NMSE*1000/QB;
}
//end F-test

/*
1. for the ith object, calculate its average distance to all other objects in its cluster. --- Ai
2. for the ith object and any other clusters not containing the object, calculate the
object's average distance to all the objects in the given cluster. Find the minimum such value
with respect to all clusters.  --- Bi
3. for the ith object, the silhouette coefficient is SCi = (Bi-Ai)/max(Ai, Bi)
4. overall measure on the clustering can obtained by taking the average of all the points.
ref: cluster analysis: basic concepts and algorithms, charpter 8.
*/
double Silhouette(TRAININGSET *pTS, CODEBOOK *pCB, PARTITIONING *pP)
{
	int		i = 0;
	int		j = 0;
	int     s = 0;
	int		t = 0;
	double  SC = 0; // finnal silhouette coefficient
	//calculate the Ai, Bi, and SCi
	for(i = 0; i < BookSize(pTS); i++)
	{
		double SCi = 0;		
		double Ai = 0;
		double Bi = DBL_MAX;
		j = Map(pP,i);
		//get Ai here;
		//first one
		for( t = FirstVector(pP, j); !EndOfPartition(t); t = NextVector(pP, t))
			Ai += sqrt(VectorDistance(Vector(pTS, i),Vector(pTS,t),
					VectorSize(pTS), MAXLLONG, EUCLIDEANSQ))/VectorFreq(pCB,j);
		//second one, can delete
/*		for(s = 0; s < BookSize(pTS); s++)
		{
			t = Map(pP,s);
			if(t == j && s!=i)
				Ai += sqrt(VectorDistance(Vector(pTS, i), Vector(pTS, s),
					VectorSize(pTS), MAXLLONG, EUCLIDEANSQ))/VectorFreq(pCB,j);					
		}*/
//		PrintMessage("Ai %f\n", Ai);

		//get Bi here
		for(s = 0; s < PartitionCount(pP); s++)
		{
			double dis = 0;
			if(s != j)
			{
				for(t = FirstVector(pP,s);!EndOfPartition(t); t = NextVector(pP, t))
					dis+= sqrt(VectorDistance(Vector(pTS, i),Vector(pTS,t),
						VectorSize(pTS), MAXLLONG, EUCLIDEANSQ))/VectorFreq(pCB,s);
			}
			if(dis != 0 && dis < Bi)
				Bi = dis;
		}
		SCi = (Bi-Ai)/max(Ai,Bi);
		SC += SCi;
	}
	SC /= BookSize(pTS);
	return SC;
}
//end for Silhouette
/*
An alternative to silhouette coefficient. use statistics of clusters instead of 
individual points.
smallest distance between two clusters;
cluster's standard deviation;
ratio of them;
ref: H. Zhu, H.L. Zhu. Clustering analysis using data range aware seeding and 
                       agglomerative expectation maximization, 2008
*/
double SilhouetteCI(TRAININGSET *pTS, CODEBOOK *pCB, PARTITIONING *pP)
{
	int		i = 0;
	int		j = 0;
	int		t = 0;
	double  SCI = 0;
	
	
	for( i = 0; i < BookSize(pCB); i++)
	{
		double  sd = 0; //standard deviation
		double  diff = DBL_MAX;
				
		for(t = FirstVector(pP,i);!EndOfPartition(t); t = NextVector(pP, t))
			sd += sqrt(VectorDistance(Vector(pTS, t),Vector(pCB,i),
					VectorSize(pTS), MAXLLONG, EUCLIDEANSQ))/VectorFreq(pCB,i);
		//min distance between clusters
		for(j = 0; j < BookSize(pCB); j++)
		{
			double temp = 0;
			if(j != i)
				temp = abs(sqrt(VectorDistance(Vector(pCB, j),Vector(pCB,i),
					VectorSize(pTS), MAXLLONG, EUCLIDEANSQ)));
			if(temp != 0 && temp < diff)
				diff = temp;
		}
		SCI += diff/sd *VectorFreq(pCB,i)/BookSize(pTS);
	}
	return SCI;

}
//end for SCI

double Xie_Beni(TRAININGSET *pTS, CODEBOOK *pCB, PARTITIONING *pP)
{
	int		   i = 0;
	int        j = 0;
	double     SE = 0;
	double     MSE = 0;
	double     CentroidMSE = 0;
	double     Xie_Beni_Index = 0;
	double     temp_dist = 0;
	double     dist = 0;
	double     min_dist = 0;
	VECTORTYPE MeanVector;
	MSE = ValidityMSE(pTS,pCB,pP);
	SE = TotalFreq(pTS) * MSE;
	// mean vector of dataset
	MeanVector = MeanOfTs(pTS,pCB);		
	for( i = 0; i < BookSize(pCB); i++)
	{		
		// EUCLIDEANSQ (without square)
		temp_dist = VectorDistance(Vector(pCB, i), MeanVector,
			VectorSize(pCB), MAXLLONG, EUCLIDEANSQ);
		CentroidMSE += temp_dist;
		CentroidMSE /= BookSize(pCB);
	}
	// min distance between clusters	
	min_dist = VectorDistance(MeanVector, Vector(pCB,1),
						VectorSize(pCB), MAXLLONG, EUCLIDEANSQ);
	for( i = 0; i < BookSize(pCB); i++)
	{
		for(j = i+1; j < BookSize(pCB); j++)
		{
			dist = VectorDistance(Vector(pCB, i), Vector(pCB, j),
					VectorSize(pCB), MAXLLONG, EUCLIDEANSQ);
			if (dist < min_dist)
				min_dist = dist;
		}
	}
	// Xie-Beni value
	Xie_Beni_Index = (VectorSize(pCB)*SE + CentroidMSE) / min_dist;
	return Xie_Beni_Index;
}
// end for Xie-Beni

/*
// XieBeni for fuzzy clustering, fuzzy partitions, different with hard clustering
double Xie_Beni_Fuzzy(TRAININGSET *pTS, CODEBOOK *pCB, FuzzyPartitioning *pP)
{


}
//end XieBeni_Fuzzy
*/

// BIC
double ValidityBIC(TRAININGSET *pTS, CODEBOOK *pCB, PARTITIONING *pP)
{
	int		   i = 0;
	int        j = 0;
	double     variance = 0;
	double     LogLikelihood = 0;
	double     BIC = 0;
	for( i = 0; i < BookSize(pCB); i++)
	{
		int   Num_Each = 0;
		int   Num = 0;
		int   dim = 0;
		for(j = 0; j < BookSize(pTS); j++)
		{
			// calculate variance of each cluster
			if(Map(pP,j) == i)
				variance += VectorDistance(Vector(pTS,j), Vector(pCB,i),
					VectorSize(pTS), MAXLLONG, EUCLIDEANSQ);
		}
		variance /= BookSize(pTS)-1;	//-BookSize(pCB)
		Num_Each = VectorFreq(pCB, i);
		Num = TotalFreq(pTS);
		dim = VectorSize(pCB);
		LogLikelihood += Num_Each*log(Num_Each)-Num_Each*log(Num)-dim*Num_Each*log(2*pi)/2
			-Num_Each*log(variance)/2-(Num_Each-BookSize(pCB))/2;
	}
	BIC = LogLikelihood - BookSize(pCB)*log(TotalFreq(pTS))/2;
	return BIC;
}
// end for BIC

double Scatter(TRAININGSET *pTS, CODEBOOK *pCB, PARTITIONING *pP)
{
	int		   i = 0;
	double     scatt = 0;

	double     varianceofcluster = 0;
	VECTORTYPE   MeanVector;
	double     VarianceofData = 0;

	MeanVector = MeanOfTs(pTS,pCB);
	
    //variance of data set
	for(i=0; i< BookSize(pTS); i++)
	{
		VarianceofData = VectorDistance(Vector(pTS, i), MeanVector,
			VectorSize(pCB), MAXLLONG, EUCLIDEANSQ);
	}

	// average of variance of clusters

	varianceofcluster = variance_cluster(pTS, pCB, pP);

	varianceofcluster /= BookSize(pCB);
	scatt = varianceofcluster/ VarianceofData;
	return scatt;
}

// sum of Variance of clusters
double variance_cluster(TRAININGSET *pTS, CODEBOOK *pCB, PARTITIONING *pP)
{
	int		   i = 0;
	int        j = 0;
	double     varianceofcluster = 0;

	// standard deviation
	for( j = 0; j < BookSize(pCB); j++)
	{
		for( i = 0; i < BookSize(pTS); i++)
		{
			if (j == Map(pP, i));
			{
				varianceofcluster += sqrt(VectorDistance(Vector(pTS, i), Vector(pCB, j), 
					VectorSize(pTS), MAXLLONG, EUCLIDEANSQ));				
			}
		}
	}
	return varianceofcluster;
}

// S_Dbw
// Understanding of Internal Clustering Validation Measures, ICDM 2010, "S_Dbw is the best."
double S_Dbw(TRAININGSET *pTS, CODEBOOK *pCB, PARTITIONING *pP)
{ 
	int		   i = 0;
	int        j = 0;
	int        k = 0;
	int        d = 0;
    double     S_Dbw = 0;
	double     Den_bw = 0;
	double     Scatt = 0;
	double     stdev = 0;
	VECTORTYPE   VectorU;
	VectorU = CreateEmptyVector(VectorSize(pCB));

	// scatter 
	Scatt = Scatter(pTS, pCB, pP);
	// standard deviation
	double temp = variance_cluster(pTS, pCB, pP);
	stdev = sqrt(temp)/BookSize(pCB);


	double fvalue = 0;
	for(i = 0; i < BookSize(pCB); i++)
	{
		for(j = 0; j < BookSize(pCB); j++)
		{
			if(j != i)
			{
				double fi = 0;
				double fj = 0;
				double fij = 0;
				double temp1 = 0;
				double temp2 = 0;
				double temp3 = 0;
				for( k = 0; k < BookSize(pTS); k++)
				{
					if(i == Map(pP, k))
					{
						temp1 = sqrt(VectorDistance(Vector(pTS, k), Vector(pCB, i),
							VectorSize(pCB), MAXLLONG, EUCLIDEANSQ));
						if(temp1 <= stdev)
						{
							fi += 1;
							//fprintf(stdout, "temp1=%f; stdev=%f,,,,", temp1, stdev);
						}
						
					}
					if(j == Map(pP, k))
					{
						temp2 = sqrt(VectorDistance(Vector(pTS, k), Vector(pCB, j),
							VectorSize(pCB), MAXLLONG, EUCLIDEANSQ));
						if(temp2 <= stdev)
						{
							fj += 1;
							//fprintf(stdout, "temp2=%f; stdev=%f,,,,", temp2, stdev);
						}
					}
				
					if(i == Map(pP, k) || j == Map(pP, k))
					{
						//take the middle point of Ci and Cj
						for( d = 0; d < VectorSize(pTS); d++)	
						{
							VectorU[d] = ((llong) VectorScalar(pCB,i,d)+ (llong) VectorScalar(pCB,j,d))/2;
							
						}
						temp3 = sqrt(VectorDistance(Vector(pTS, k), VectorU,
							VectorSize(pCB), MAXLLONG, EUCLIDEANSQ));
						if(temp3 <= stdev)
						{
							fij += 1;
							//fprintf(stdout, "temp3=%f; stdev=%f\n ", temp3, stdev);
						}
					}
				}
				if((fij != 0) & (max(fi,fi) != 0)){
					fvalue += fij/max(fi, fj);
					fprintf(stdout, "fi=%f; fj=%f; fij=%fl fvalue = %f\n ", fi, fj, fij, fij/max(fi, fj));
				}
				
			} // i!=j
		} // for j
	} // for i
	Den_bw = fvalue/((BookSize(pCB)* (BookSize(pCB)-1)));
	fprintf(stdout, "Scatt=%f; Den_bw=%f;\n", Scatt, Den_bw);
	S_Dbw = Scatt + Den_bw;
	return S_Dbw;
}
//end S_Dbw



// ----------External index
/* Need two partitions, which are pP and pPtruth here. If there is no
ground truth partitions, then we can use MonteCarlo method.
get pairs of points and
ss - they both belong to pP and pPtruth
sd - they both belong to pP but not pPtruth
ds - they don't belong to pP but belong to pPtruth
dd - not belong to pP and pPtruth
count the numbers of each
*/
double RandIndex(TRAININGSET *pTS, PARTITIONING *pP, PARTITIONING *pPtruth)
{
	double ss = 0;  
	double sd = 0;
	double ds = 0;
	double dd = 0;
	int i = 0;
	int j = 0;
	double indexvalue = 0;
	for(i = 1; i < BookSize(pTS)+1; i++)
	{
		for(j = i+1; j < BookSize(pTS)+1; j++)
		{
			
			if(Map(pP, j) == Map(pP, i) && Map(pPtruth, j) == Map(pPtruth, i))
				ss++;
			if(Map(pP, j) == Map(pP, i) && Map(pPtruth, j) != Map(pPtruth, i))
				sd++;
			if(Map(pP, j) != Map(pP, i) && Map(pPtruth, j) == Map(pPtruth, i))
				ds++;
			if(Map(pP, j) != Map(pP, i) && Map(pPtruth, j) != Map(pPtruth, i))
				dd++;
			
		}
	}

	indexvalue = (ss+dd)/(ss+sd+ds+dd);
	return indexvalue;
}
// corrected rand index = [(a+d)-((a+b)(a+c)+(c+d)(b+d))*1/m ] /[m-((a+b)(a+c)+(c+d)(b+d))*1/m]
double CorrectedRI(TRAININGSET *pTS, PARTITIONING *pP, PARTITIONING *pPtruth)
{
	double ss = 0;  
	double sd = 0;
	double ds = 0;
	double dd = 0;
	double m = 0;
	double part1 = 0;
	double part2 = 0;
	int i = 0;
	int j = 0;
	double indexvalue = 0;
	for(i = 1; i < BookSize(pTS)+1; i++)
	{
		for(j = i+1; j < BookSize(pTS)+1; j++)
		{
			
			if(Map(pP, j) == Map(pP, i) && Map(pPtruth, j) == Map(pPtruth, i))
				ss++;
			if(Map(pP, j) == Map(pP, i) && Map(pPtruth, j) != Map(pPtruth, i))
				sd++;
			if(Map(pP, j) != Map(pP, i) && Map(pPtruth, j) == Map(pPtruth, i))
				ds++;
			if(Map(pP, j) != Map(pP, i) && Map(pPtruth, j) != Map(pPtruth, i))
				dd++;
			
		}
	}
	m = ss+sd+ds+dd;
	part1 = (ss+dd) -((ss+sd)*(ss+ds)+(ds+dd)*(sd+dd))/m;
	part2 = m - ((ss+sd)*(ss+ds)+(ds+dd)*(sd+dd))/m;
	indexvalue = part1/part2;
	return indexvalue;
}

// function for normalized mutual information
double NMI(TRAININGSET *pTS, CODEBOOK *pCB, PARTITIONING *pP, PARTITIONING *pPtruth)
{
	double entropy1 = 0;
	double entropy2 = 0;
	double MI = 0;
	double indexvalue = 0;
	int k = 0;
	int l = 0;
	int m = 0;
	int n = 0;
	int j = 0;
	
// calculate mutual information	
	for (k = 0; k < BookSize(pCB); k++)
		for( l = 0; l < BookSize(pCB); l++)
		{
			double c1 = 0;
			double c2 = 0;
			double c3 = 0;
			for(j = 0; j< BookSize(pTS); j++){
				if(Map(pP, j) ==k)
					c1++;
				if(Map(pPtruth, j) == l)
					c2++;	
				if(Map(pP, j) ==k && Map(pPtruth, j) == l)
					c3++;
			}
			c1 /= BookSize(pTS);
			c2 /= BookSize(pTS);
			c3 /= BookSize(pTS);
			if(c3 != 0)
				MI += c3*log(c3/(c1*c2));
		}
// calculate the entropy sum for normalization.
	for( m = 0; m < BookSize(pCB); m++)
	{
		double counts1 = 0;
		double counts2 = 0;

		for( n = 0; n < BookSize(pTS); n++)
		{
			if (Map(pP, n) == m)
				counts1++;
			if (Map(pPtruth, n) == m)
				counts2++;
		}
		counts1 /= BookSize(pTS);
		counts2 /= BookSize(pTS);
		entropy1 += counts1 * log(counts1);
		entropy2 += counts2 * log(counts2);
	}

	//normalized mutual information (NMI): 1 represents perfert agreement
//	indexvalue = MI/sqrt(entropy1*entropy2);
	//variation of information (VI): 0 represents perfert agreement
	indexvalue = 0 - entropy1 - entropy2 - 2* MI;
	return indexvalue;
}

/*-------------------------------------------------------------------*/
/*-------------------------------------------------------------------*/
/*---------------------end clustering validity-----------------------*/


