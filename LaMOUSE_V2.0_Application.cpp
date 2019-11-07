
#include <iostream>
#include <math.h>
#include <stdlib.h>
#include <time.h>
#include <fstream>
#include <string>
#include <ctime>
#include <random>  

#include<omp.h>



//La'MOUSE
//Version 2.0
//This version integrates different habitats

//Created December 26, 2018
//Last revision 6/19/2019


using namespace std;

time_t seconds;

//Base variables and parameters
int t,Test,i,j,DataCounter, bWindow,bTau,AlphaWindow,AlphaTau,MJAWindow,MJATau;
double CurrentTime,TimeToRun,TotalReactions,TimeOfReaction,UniSeed,Sum,Picker;
double b[10],bMax, dN,dJ,dA,MNJ,MJA, Ro,MJARo, Alpha[10],AlphaSlope,AlphaIntercept, Mu,AlphaMu,MJAMu,mj,ma,TrapRateN,TrapRateJ,TrapRateA,Lambda;
double AlphaMin, bMin[10], MJAMin, MJAMax,AlphaMax,bThresh,AlphaThresh,MJAThresh;
double Reaction[1000], ReactionProbability[1000];
int Sites,WhichReaction, NumReacts, StartTime;
int SN[10], SJ[10], SA[10],IN[10], IJ[10], IA[10], RN[10], RJ[10], RA[10];
double BA[10], BM, BN, g, M,V;

//ABC variables follow
int T,Trial,MaxTrials,Hits,CurrentDay,LastRecordDay,Extinct,DaysElapsed, TrapNights[10][100000];
double Rain[10][100000],LST[10][100000],NDVI[10][100000], TotalSJ[10][100000], TotalIJ[10][100000], TotalRJ[10][100000], TotalSA[10][100000], TotalIA[10][100000], TotalRA[10][100000],  TotalSimSN[10][100000], TotalSimIN[10][100000], TotalSimRN[10][100000], TotalSimSJ[10][100000], TotalSimIJ[10][100000],TotalSimRJ[10][100000],TotalSimSA[10][100000], TotalSimIA[10][100000], TotalSimRA[10][100000],TotalRats[10][100000],TotalSimRats[10][100000],TotalJuv[10][100000],TotalAdu[10][100000],TotalSimJuv[10][100000],TotalSimAdu[10][100000],JPCR[10][100000],APCR[10][100000],JSero[10][100000],ASero[10][100000];
double Data[10][10000][200], AverageRain,Day[10][10000];
double R0NG[10][10000], R0OS[10][10000], MeanR0NG[10], MeanR0OS[10];
int SNCaptured[10][10000], INCaptured[10][10000], RNCaptured[10][10000], SJCaptured[10][10000], IJCaptured[10][10000], RJCaptured[10][10000], SACaptured[10][10000], IACaptured[10][10000], RACaptured[10][10000], TotalJuvCap[10][10000], TotalAduCap[10][10000];

//debugging vars
int AddUpRats,AddUpSN,AddUpIN,AddUpRN,AddUpSJ,AddUpIJ,AddUpRJ,AddUpSA,AddUpIA,AddUpRA, NumAdds;
double MeanRats,MeanSN,MeanIN,MeanRN,MeanSJ,MeanIJ,MeanRJ,MeanSA,MeanIA,MeanRA;

//parallel vars
int Files, NumFiles;
//char StringNum[50];

ifstream in_Data[10];
//ofstream out_Densities[10];
//ofstream out_Pars[10];
ofstream out_Posterior;
//ofstream out_DSCompare[10];
//ofstream out_Testing[10];
//ofstream out_ErrorLog[10];

random_device rd;   // non-deterministic generator  
mt19937 gen(rd());  // to seed mersenne twister.  

//DEFINE FUNCTIONS
double Big_Rand(double MIN, double MAX)
{
	uniform_real_distribution<> dist(MIN, MAX); // distribute results between Min and Max
	return dist(gen);
}

int Rand_Int(int MIN, int MAX)
{
	uniform_int_distribution<> dist(MIN, MAX); // distribute results between Min and Max
	return dist(gen);
}

int Fish_Rand(double p)
{
	poisson_distribution<> distr(p); // distribute results between 1 and 6 inclusive.  
	return distr(gen);
}

int Bin_Rand(int Tries, double Srate)
{
	binomial_distribution<> distr(Tries,Srate); // distribute results between 1 and 6 inclusive.  
	return distr(gen);
}

double Norm_Rand(double Mu, double SD)
{ 
	normal_distribution<> distr(Mu, SD); // 
	return distr(gen);
}

double LogNorm_Rand(double p, double q)
{ 
	lognormal_distribution<> distr(p, q); // p is Alpha (shape) and q is Beta (rate)
	return distr(gen);
}

double Gam_Rand(double p, double q)
{
	gamma_distribution<> distr(p, q); // p is Alpha (shape) and q is Beta (rate)
	return distr(gen);
}

//DONE DEFINING FUNCTIONS

void CalculateRates();
void PickReaction();
void ImplementReaction();
void SummaryStats();
void RecordData();
void Prior();
void IngestData();
void Initialize();
void IntegrateBinaryEnvironment();

//#define _CRT_SECURE_NO_WARNINGS


int main(int argc, char *argv[])
{
	
	
	//cout << "Enter the number of sites\n";
	//cin >> Sites;
	Sites = 2;
	StartTime = 365;//set to 365 to allow space for back averaging of rainfall

	//string DataName[10];


	//Use this for application to Faranah data
	in_Data[0].open("BantouData.tsv");
	in_Data[1].open("TanganyaData.tsv");
	//string OutIndex;
	//cout << "Enter output file index #\n";//This should be a number only. i.e., first run is 1, second run is 2, etc...
	//cin >> OutIndex;
	string OutIndex = argv[1];

	//out_Densities[0].open("BantouDensities"+OutIndex+ ".csv");
	//out_Densities[1].open("TanganyaDensities" + OutIndex + ".csv");

	//out_Pars[0].open("BantouPars" + OutIndex + ".csv");
	//out_Pars[1].open("TanganyaPars" + OutIndex + ".csv");

	//out_DSCompare[0].open("BantouDSCompare" + OutIndex + ".csv");
	//out_DSCompare[1].open("TanganyaDSCompare" + OutIndex + ".csv");

	//out_Testing[0].open("BantouTesting" + OutIndex + ".csv");
	//out_Testing[1].open("TanganyaTesting" + OutIndex + ".csv");

	//out_ErrorLog[0].open("BantouError" + OutIndex + ".csv");
	//out_ErrorLog[1].open("TanganyaError" + OutIndex + ".csv");

	out_Posterior.open("FaranahPosterior"+ OutIndex + ".csv");

	//Ingest rainfall data **********Must be tab delimited .tsv************
	IngestData();

	//out_ErrorLog[0] << "TotAbundDist,TotalFreqDist\n";
	//out_ErrorLog[1] << "DifMeanJuv,DifMeanAdu,CorJuvJuv,CorAduAdu,DifFreqIJ,DifFreqRJ,DifFreqIA,DifFreqRA\n";
	//out_Densities[0] << "CurrentTime,CurrentDay,Rain,LST,NDVI,AverageRain,b,SN,SJ,SA,IN,IJ,IA,RN,RJ,RA\n";
	//out_Densities[1] << "CurrentTime,CurrentDay,Rain,LST,NDVI,AverageRain,b,SN,SJ,SA,IN,IJ,IA,RN,RJ,RA\n";
	//out_Pars[0] << "bWindow, bTau , bThresh, bMin, bMax, MNJ, MJA, Alpha, dN, dJ, dA,BA,BM,BN,g,V,M, mj,ma, TrapRateJ,TrapRateA,Lambda\n";
	//out_Pars[1] << "bWindow, bTau , bThresh, bMin, bMax, MNJ, MJA, Alpha, dN, dJ, dA,BA,BM,BN,g,V,M, mj,ma, TrapRateJ,TrapRateA,Lambda\n";
	out_Posterior << "bWindow, bTau , bThresh, bMin(0), bMin(1), bMax, MNJ, MJA,Alpha(0),Alpha(1),dN, dJ, dA,BA(0),BA(1),BM,BN,g,V,M, mj,ma, TrapRateJ,TrapRateA,Lambda,R01_NG,R02_NG\n";

	//out_DSCompare <<" Day,  Rain , LST[k] , NDVI, TotalSJ ,TotalSimSJ, TotalIJ, TotalSimIJ ,TotalRJ ,TotalSimRJ,TotalSA ,TotalSimSA, TotalIA , TotalSimIA ,TotalRA,TotalSimRA \n";
	//out_DSCompare << "Day, Rain,LST,NDVI ,TotalRats,TotalSimRats\n";
	//out_DSCompare << "Day, Rain,LST,NDVI ,Juv,SimJuv,Adu,SimAdu\n";
	//out_DSCompare << "Day, Rain,LST,NDVI ,Juv,SimJuv,Adu,SimAdu,Data_FIJ,Sim_FIJ,Data_FRJ,Sim_FRJ,Data_FIA,Sim_FIA,Data_FRA,Sim_FRA\n";
	//out_DSCompare[0] << "Day,Rain,LST,NDVI,TotalJuv,TotalSimJuv,SimJuvCap,TotalAdu,TotalSimAdu,SimAduCap,TotalSJ,TotalSimSJ,SimSJCap,TotalIJ,TotalSimIJ,SimIJCap,TotalRJ,TotalSimRJ,SimRJCap,TotalSA,TotalSimSA,SimSACap,TotalIA, TotalSimIA,SimIACap,TotalRA,TotalSimRA,SimRACap\n";
	//out_DSCompare[1] << "Day,Rain,LST,NDVI,TotalJuv,TotalSimJuv,SimJuvCap,TotalAdu,TotalSimAdu,SimAduCap,TotalSJ,TotalSimSJ,SimSJCap,TotalIJ,TotalSimIJ,SimIJCap,TotalRJ,TotalSimRJ,SimRJCap,TotalSA,TotalSimSA,SimSACap,TotalIA, TotalSimIA,SimIACap,TotalRA,TotalSimRA,SimRACap\n";

	//START ABC LOOP HERE

	Hits = 0;//Hit counter for posterior
	Trial = 0;

	while (Hits < 50)
	{
		cout << "trials:" << Trial << " Hits:" << Hits << "\n";

		//Draw parameters from priors
		Prior();

		//Initialize
		Initialize();//Sets initial densities and the "CurrentTime"


		T = 0;//This is just a counter for data recording and later use in ABC
	
		LastRecordDay = CurrentTime-1;//An indicator variable that helps record data only on new integer days. Prevents overflow of record arrays. Subtract 1 so that the first day is recorded (initialize sets current time to 365)
		
		Extinct = 0;//An indicator variable used to identify cases where mice go extinct prior to end of run. These are then excluded from the posterior.

		AddUpSN = 0;//THIS IS A DEBUGGING TOOL ONLY
		AddUpIN = 0;//THIS IS A DEBUGGING TOOL ONLY
		AddUpRN = 0;//THIS IS A DEBUGGING TOOL ONLY
		AddUpSJ = 0;//THIS IS A DEBUGGING TOOL ONLY
		AddUpIJ = 0;//THIS IS A DEBUGGING TOOL ONLY
		AddUpRJ = 0;//THIS IS A DEBUGGING TOOL ONLY
		AddUpSA = 0;//THIS IS A DEBUGGING TOOL ONLY
		AddUpIA = 0;//THIS IS A DEBUGGING TOOL ONLY
		AddUpRA = 0;//THIS IS A DEBUGGING TOOL ONLY
		AddUpRats = 0;
		NumAdds = 0;//THIS IS A DEBUGGING TOOL ONLY
		do//start time loop
		{

			CurrentDay = int(CurrentTime);//Discretize time so that daily rainfall can be used

			//cout << CurrentTime << "\n";
			//cout << CurrentDay << "\n";
			//cin >> Test;

			//Calculate summary values for environmental variables and use them to calculate current values of environmentally dependent model parameters
			//IntegrateSimpleEnvironment();
			IntegrateBinaryEnvironment();

			/*cout << bMax << "," << Ro << "," << AverageRain << "," << Mu << "\n";
			cout << b << "\n";
			cin >> Test;
			*/

			//Execute Reactions
			CalculateRates();
			PickReaction();
			ImplementReaction();
			//Reactions complete

			//THE FOLLOWING CHUNK OF CODE IS FOR DEBUGGING ONLY
			/*
			AddUpSN = AddUpSN + SN;
			AddUpSJ = AddUpSJ + SJ;
			AddUpSA = AddUpSA + SA;
			AddUpIN = AddUpIN + IN;
			AddUpIJ = AddUpIJ + IJ;
			AddUpIA = AddUpIA + IA;
			AddUpRN = AddUpRN + RN;
			AddUpRJ = AddUpRJ + RJ;
			AddUpRA = AddUpRA + RA;
			AddUpRats = AddUpRats + SN + SJ + SA + IN + IJ + IA + RN + RJ + RA;
			NumAdds++;

			if (SN < 0 || SJ < 0 || SA < 0 || IN < 0 || IJ < 0 || IA < 0 || RN < 0 || RJ < 0 || RA < 0)
			{
				cout << SN << "," << SJ << "," << SA << "," << IN << "," << IJ << "," << IA << "\n";
				cin >> Test;
			}
			//END OF DEBUGGING CODE*/

			//cout << SN << "," << SJ << "," << SA << "," << IN << "," << IJ << "," << IA << "\n";

		
			//Output densities (generally shut down to increase speed)
			/*
			if (CurrentDay % 200 == 0 && CurrentDay > LastRecordDay)
			{
				for (i = 0; i < Sites; i++)
				{
					out_Densities[i] << CurrentTime << "," << CurrentDay << "," << Rain[i][CurrentDay] << "," << LST[i][CurrentDay] << "," << NDVI[i][CurrentDay] << "," << AverageRain << "," << b[i] << "," << SN[i] << "," << SJ[i] << "," << SA[i] << "," << IN[i] << "," << IJ[i] << "," << IA[i] << "," << RN[i] << "," << RJ[i] << "," << RA[i] << "\n";
					out_Densities[i].flush();
				}
			}*/

			//Record simulated and real data on an (integer) daily basis
			DaysElapsed = CurrentDay - LastRecordDay;//Count how many days past between reactions. We need to force data recording on a daily basis even if multiple days pass before a new reaction occurs
			if (CurrentDay > LastRecordDay)//Record simulated data no more than once per day
			{
				RecordData();
				LastRecordDay = int(CurrentTime);
			}
			//Done with ABC things

		
			//cout << CurrentTime << "\n";
			//cout << CurrentDay << "\n";
			//cin >> Test;

			//Update time to next reaction
			UniSeed = 0;
			while (UniSeed == 0)
			{
				UniSeed = Big_Rand(0, 1);
			}
			TimeOfReaction = -log(UniSeed) / TotalReactions;

			CurrentTime = CurrentTime + TimeOfReaction;

		} while (CurrentTime < TimeToRun);//End time loop

		//THE BELOW IS ALL JUST DEBUGGING CODE
		MeanSN = AddUpSN / (1.0*NumAdds);
		MeanIN = AddUpIN / (1.0*NumAdds);
		MeanRN = AddUpRN / (1.0*NumAdds);
		MeanSJ = AddUpSJ / (1.0*NumAdds);
		MeanIJ = AddUpIJ / (1.0*NumAdds);
		MeanRJ = AddUpRJ / (1.0*NumAdds);
		MeanSA = AddUpSA / (1.0*NumAdds);
		MeanIA = AddUpIA / (1.0*NumAdds);
		MeanRA = AddUpRA / (1.0*NumAdds);
		MeanRats = AddUpRats / (1.0*NumAdds);
		//out_Testing[0] << BN << "," << BM << "," << BA << "," << MeanRats << "," << MeanSN << "," << MeanIN << "," << MeanRN << "," << MeanSJ << "," << MeanIJ << "," << MeanRJ << "," << MeanSA << "," << MeanIA << "," << MeanRA << "\n";//DEBUGGING CODE ONLY. TO BE REMOVED



		if (CurrentDay==TimeToRun-1&&Extinct != 1)//calculate summary stats and test acceptance only if the simulation ran long enough to actually compare simulated and real data.
		{

			//cout << "Calling summaries next\n";
			//cin >> Test;
			SummaryStats();

		}

		Trial++;
	}//END of ABC trials loop


}



void CalculateRates()
{
	int TotRat[10],TotInfect[10],CR;

	CR = 0;//CR is just a reaction counter
	for (i = 0; i < Sites; i++)//for the number of habitats/sites run these reactions
	{
		TotRat[i] = SN[i] + IN[i] + RN[i] + SJ[i] + SA[i] + IJ[i] + IA[i] + RJ[i] + RA[i];//Sum up the density of **ALL** rats (changed from just free-living rats on May 27, 19). This is used for density dependence (assuming only free-living forms compete for resources) and for frequency dependent infection (assuming you only run into other free living individuals)
		TotInfect[i] = IN[i] + IJ[i] + IA[i];//Sum up the number of infected rats

		if (SN[i] < 0 || IN[i] < 0 || RN[i] < 0 || SJ[i] < 0 || SA[i] < 0 || IJ[i] < 0 || IA[i] < 0 || RJ[i] < 0 || RA[i] < 0)//If any density goes negative, throw an error and await acknowledgement
		{
			cout << "WARNING: Negative densities have occured\n";
			cout << SN[i] << "," << IN[i] << "," << RN[i] << "," << SJ[i] << "," << SA[i] << "," << IJ[i] << "," << IA[i] << "," << RJ[i] << "," << RA[i] << "\n";
			cin >> Test;
		}


		//cout << "TotalRat:" << TotRat[i];
		//cin >> Test;
		//Base MOUSE
		Reaction[CR] = b[i] * SA[i];//A baby rat is born
		CR++;
		Reaction[CR] = dN * SN[i];//A baby rat dies
		CR++;
		Reaction[CR] = Alpha[i] * SN[i] * (TotRat[i]);//Competition kills a baby rat
		CR++;
		Reaction[CR] = MNJ * SN[i];//A baby rat matures
		CR++;
		Reaction[CR] = dJ * SJ[i];//A juvenile rat dies
		CR++;
		Reaction[CR] = Alpha[i] * SJ[i] * (TotRat[i]);//Competition kills a juvenile rat
		CR++;
		Reaction[CR] = MJA * SJ[i];//A juvenile rat matures
		CR++;
		Reaction[CR] = dA * SA[i];//An adult rat dies
		CR++;
		Reaction[CR] = Alpha[i] * SA[i] * (TotRat[i]);//Competition kills an adult rat.
		CR++;


		//Additions for LaMOUSE
		Reaction[CR] = (BA[i] / (1.0 + Lambda * (TotRat[i] - 1.0))) * SA[i] * IA[i];//An adult is infected by an infectious adult
		CR++;
		Reaction[CR] = (BA[i] / (1.0 + Lambda * (TotRat[i] - 1.0))) * SA[i] * IJ[i];//An adult is infected by an infectious juvenile
		CR++;
		Reaction[CR] = (BA[i] / (1.0 + Lambda * (TotRat[i] - 1.0))) * SJ[i] * IA[i];//A juvenile is infected by an infectious adult
		CR++;
		Reaction[CR] = (BA[i] / (1.0 + Lambda * (TotRat[i] - 1.0)))* SJ[i] * IJ[i];//A juvenile is infected by an infectious juvenile
		CR++;
		Reaction[CR] = BM * SN[i] * IA[i];//A newborn is infected by an infectious adult (presumably its mother)
		CR++;
		Reaction[CR] = BN * SN[i] * IN[i];//A newborn is infected by an infectious newborn (presumably its nest mate)
		CR++;

		Reaction[CR] = g * IN[i];//An infected newborn recovers
		CR++;
		Reaction[CR] = g * IJ[i];//An infected juvenile recovers
		CR++;
		Reaction[CR] = g * IA[i];//An infected adult recovers
		CR++;

		Reaction[CR] = dN * IN[i];//An infected newborn dies
		CR++;
		Reaction[CR] = dJ * IJ[i];//An infected juvenile dies
		CR++;
		Reaction[CR] = dA * IA[i];//And infected adult dies
		CR++;
		Reaction[CR] = dN * RN[i];//A recovered newborn dies
		CR++;
		Reaction[CR] = dJ * RJ[i];//A recovered juvenile dies
		CR++;
		Reaction[CR] = dA * RA[i];//A recovered adult dies
		CR++;

		Reaction[CR] = b[i] * V* IA[i];//An infected adult gives birth to an infected baby (vertical transmission)
		CR++;
		Reaction[CR] = b[i] * (1.0 - V)* IA[i];//An infected adult gives birth to an UNinfected baby (no vertical transmission)
		CR++;
		Reaction[CR] = b[i] * M * RA[i];//A recovered adult gives birth to a recovered baby (maternal antibody transfer)
		CR++;
		Reaction[CR] = b[i] * (1.0 - M)* RA[i];//A recovered adult gives birth to a susceptible baby (no maternal antibody transfer)
		CR++;

		Reaction[CR] = Alpha[i] * IN[i] * (TotRat[i]);//Competition kills an infected Newborn
		CR++;
		Reaction[CR] = Alpha[i] * RN[i] * (TotRat[i]);//Competition kills a recovered Newborn
		CR++;
		Reaction[CR] = Alpha[i] * IJ[i] * (TotRat[i]);//Competition kills an infected juvenile
		CR++;
		Reaction[CR] = Alpha[i] * RJ[i] * (TotRat[i]);//Competition kills a recovered juvenile
		CR++;
		Reaction[CR] = Alpha[i] * IA[i] * (TotRat[i]);//Competition kills an infected adult
		CR++;
		Reaction[CR] = Alpha[i] * RA[i] * (TotRat[i]);//Competition kills a recovered adult
		CR++;

		Reaction[CR] = MNJ * IN[i];//A baby infected rat matures
		CR++;
		Reaction[CR] = MNJ * RN[i];//A baby recovered rat matures
		CR++;
		Reaction[CR] = MJA * IJ[i];//A juvenile infected rat matures
		CR++;
		Reaction[CR] = MJA * RJ[i];//A juvenile recovered rat matures
		CR++;

		//New migration based on probability of departing
		Reaction[CR] = mj * SJ[i];//A juvenile susceptible rat departs for greener pastures
		CR++;
		Reaction[CR] = ma * SA[i];//An adult susceptible rat departs for greener pastures
		CR++;
		Reaction[CR] = mj * IJ[i];//A juevnile infected rat departs for greener pastures
		CR++;
		Reaction[CR] = ma * IA[i];//An adult infected rat departs for greener pastures
		CR++;
		Reaction[CR] = mj * RJ[i];//A Juvenile resistant rat departs for greener pastures
		CR++;
		Reaction[CR] = ma * RA[i];//An adult resistant rat departs for greener pastures
		CR++;
	}
	NumReacts = CR;
	/*
	cout << "Here is a list of reaction probabilities:\n";
	cout << Reaction[0] << "," << Reaction[1] << ", " << Reaction[1] << "\n";
	cin >> Test;
	*/

	TotalReactions = 0;
	for (i = 0; i < NumReacts; i++)
	{
		TotalReactions = TotalReactions + Reaction[i];
	}
	for (i = 0; i <= Sites; i++)//for the number of habitats/sites run these reactions
	{
		if (TotalReactions == 0&&TotRat[i]>0)//if there are no more reactions...
		{
			cout << "WARNING: No REACTION PROBABILITY\n";
		}
	}

	if (TotalReactions < 0)//if reactions are negative throw an error and await acknowledgement
	{
		cout << "WARNING: Negative reactions!!!";
		cout << b[0] << "," << b[1] << "\n";
		cin >> Test;
	}


	if (TotalReactions > 0)//if total reactions is positive, use it to calculate probability of each reaction
	{
		for (i = 0; i < NumReacts; i++)
		{
			ReactionProbability[i] = Reaction[i] / TotalReactions;
			//cout << ReactionProbability[i] << "\n";
		}
	}
	for (i = 0; i <Sites; i++)//for the number of habitats/sites run these reactions
	{
		if (TotRat[i] <= 0)//end simulation if rats or infected are extinct
		{
			cout << "EXTINCTION\n";
			Extinct = 1;
			CurrentTime = TimeToRun + 1;
		}
	}
}

void PickReaction()
{
	double Low, Hi, GlobalRats;
	int CheckReact,k;

	
	//Draw a random number to identfy which reaction occurs
	Low = 0;
	Hi = 0;
	Picker = Big_Rand(0, 1);
	CheckReact = -1;
	for (i = 0; i < NumReacts; i++)
	{
		Hi = Low + ReactionProbability[i];

		if (Picker >= Low && Picker < Hi)
		{
			WhichReaction = i;
			CheckReact = 1;
		}

		Low = Hi;
	}
	GlobalRats = 0;
	for (k = 0; k < Sites; k++)
	{
		GlobalRats = GlobalRats + SN[k] + SJ[k] + SA[k] + IN[k] + IJ[k] + IA[k] + RN[k] + RJ[k] + RA[k];
	}
	//cout << "GlobalRats:" << GlobalRats << "\n";
	if (CheckReact == -1 && GlobalRats > 0)//if no reaction occurs, and the rats are not extinct, throw an error message and wait for acknowledgement
	{
		cout << "No reaction occurred\n";
		cout << Picker << "\n";
		cout << Picker - Hi << "\n";
		cin >> Test;
	}
	
	
	//	cout << WhichReaction << "\n";
	//	cin >> Test;

}

void ImplementReaction()
{
	int CR;

	CR = 0;//CR is just a reaction counter
	for (i = 0; i <= Sites; i++)//for the number of habitats run these reactions
	{
		//Standard MOUSE reactions
		if (WhichReaction == CR)//Birth a rat
		{
			SN[i]++;
		}
		CR++;
		if (WhichReaction == CR)//Kill a baby rat (DI)
		{
			SN[i]--;
		}
		CR++;
		if (WhichReaction == CR)//Competition kills a baby rat
		{
			SN[i]--;
		}
		CR++;
		if (WhichReaction == CR)//Mature a baby rat
		{
			SN[i]--;
			SJ[i]++;
		}
		CR++;
		if (WhichReaction == CR)//Juvenile dies
		{
			SJ[i]--;
		}
		CR++;
		if (WhichReaction == CR)//Juvenile dies (DD)
		{
			SJ[i]--;
		}
		CR++;
		if (WhichReaction == CR)//Juvenile matures
		{
			SJ[i]--;
			SA[i]++;
		}
		CR++;
		if (WhichReaction == CR)//Adult dies (DI)
		{
			SA[i]--;
		}
		CR++;
		if (WhichReaction == CR)//Adult dies (DD)
		{
			SA[i]--;
		}
		CR++;

		//Below are LaMOUSE specific reactions
		if (WhichReaction == CR)//An adult is infected by an infectious adult
		{
			SA[i]--;
			IA[i]++;
		}
		CR++;
		if (WhichReaction == CR)//An adult is infected by an infectious juvenile
		{
			SA[i]--;
			IA[i]++;
		}
		CR++;
		if (WhichReaction == CR)//A juvenile is infected by an infectious adult
		{
			SJ[i]--;
			IJ[i]++;
		}
		CR++;
		if (WhichReaction == CR)//A juvenile is infected by an infectious juvenile
		{
			SJ[i]--;
			IJ[i]++;
		}
		CR++;
		if (WhichReaction == CR)//A newborn is infected by an infectious adult (presumably its mother)
		{
			SN[i]--;
			IN[i]++;
		}
		CR++;
		if (WhichReaction == CR)//A newborn is infected by an infectious newborn (presumably its nest mate)
		{
			SN[i]--;
			IN[i]++;
		}
		CR++;
		if (WhichReaction == CR)//An infected newborn recovers
		{
			IN[i]--;
			RN[i]++;
		}
		CR++;
		if (WhichReaction == CR)//An infected juvenile recovers
		{
			IJ[i]--;
			RJ[i]++;
		}
		CR++;
		if (WhichReaction == CR)//An infected adult recovers
		{
			IA[i]--;
			RA[i]++;
		}
		CR++;
		if (WhichReaction == CR)//An infected newborn dies
		{
			IN[i]--;
		}
		CR++;
		if (WhichReaction == CR)//An infected juvenile dies
		{
			IJ[i]--;
		}
		CR++;
		if (WhichReaction == CR)//An infected adult dies
		{
			IA[i]--;
		}
		CR++;
		if (WhichReaction == CR)//A recovered newborn dies
		{
			RN[i]--;
		}
		CR++;
		if (WhichReaction == CR)//A recovered juvenile dies
		{
			RJ[i]--;
		}
		CR++;
		if (WhichReaction == CR)//A recovered adult dies
		{
			RA[i]--;
		}
		CR++;
		if (WhichReaction == CR)//An infected adult gives birth to an infected baby (vertical trans)
		{
			IN[i]++;
		}
		CR++;
		if (WhichReaction == CR)//An infected adult gives birth to an UNinfected baby (no vertical trans)
		{
			SN[i]++;
		}
		CR++;
		if (WhichReaction == CR)//A recovered adult gives birth to a recovered baby (vertical trans of antibody)
		{
			RN[i]++;
		}
		CR++;
		if (WhichReaction == CR)//A recovered adult gives birth to a susceptible baby (no vertical trans of antibody)
		{
			SN[i]++;
		}
		CR++;
		if (WhichReaction == CR)//Competition kills an infected baby
		{
			IN[i]--;
		}
		CR++;
		if (WhichReaction == CR)//Competition kills a recovered baby
		{
			RN[i]--;
		}
		CR++;
		if (WhichReaction == CR)//Competition kills an infected juvenile
		{
			IJ[i]--;
		}
		CR++;
		if (WhichReaction == CR)//Competition kills a recovered juvenile
		{
			RJ[i]--;
		}
		CR++;
		if (WhichReaction == CR)//Competition kills an infected adult
		{
			IA[i]--;
		}
		CR++;
		if (WhichReaction == CR)//Competition kills a recovered adult
		{
			RA[i]--;
		}
		CR++;
		if (WhichReaction == CR)//A baby infected rat matures
		{
			IN[i]--;
			IJ[i]++;
		}
		CR++;
		if (WhichReaction == CR)//A baby recovered rat matures
		{
			RN[i]--;
			RJ[i]++;
		}
		CR++;
		if (WhichReaction == CR)//A juvenile infected rat matures
		{
			IJ[i]--;
			IA[i]++;
		}
		CR++;
		if (WhichReaction == CR)//A juvenile recovered rat matures
		{
			RJ[i]--;
			RA[i]++;
		}
		CR++;
		//IMMIGRATION
		if (WhichReaction == CR)//SJ leaves
		{
			SJ[i]--;
			if (i == 0)
			{
				SJ[1]++;
			}
			if (i == 1)
			{
				SJ[0]++;
			}
		}
		CR++;
		if (WhichReaction == CR)//SA leaves
		{
			SA[i]--;
			if (i == 0)
			{
				SA[1]++;
			}
			if (i == 1)
			{
				SA[0]++;
			}
		}
		CR++;
		if (WhichReaction == CR)//IJ leaves
		{
			IJ[i]--;
			if (i == 0)
			{
				IJ[1]++;
			}
			if (i == 1)
			{
				IJ[0]++;
			}
		}
		CR++;
		if (WhichReaction == CR)//IA leaves
		{
			IA[i]--;
			if (i == 0)
			{
				IA[1]++;
			}
			if (i == 1)
			{
				IA[0]++;
			}
		}
		CR++;
		if (WhichReaction == CR)//RJ leaves
		{
			RJ[i]--;
			if (i == 0)
			{
				RJ[1]++;
			}
			if (i == 1)
			{
				RJ[0]++;
			}
		}
		CR++;
		if (WhichReaction == CR)//RA leaves
		{
			RA[i]--;
			if (i == 0)
			{
				RA[1]++;
			}
			if (i == 1)
			{
				RA[0]++;
			}
		}
		CR++;
	}



}

void RecordData()//This function records simulated data every day and backfills any days that were missed due to Gillespie jumping over days. Guarantees we have daily data to compare with real data
{
	double CaptureRate;
	double JuvFrac, AduFrac;

	//This was changed 5/23/19 to only record daily simulated data and backfill between missing days

	for (i = 0; i < Sites; i++)
	{
		//Record the current values of the variables and instantaneous (simplified) R0
		//newborns
		TotalSimSN[i][CurrentDay] = SN[i];
		TotalSimIN[i][CurrentDay] = IN[i];
		TotalSimRN[i][CurrentDay] = RN[i];
		//juveniles
		TotalSimSJ[i][CurrentDay] = SJ[i];
		TotalSimIJ[i][CurrentDay] = IJ[i];
		TotalSimRJ[i][CurrentDay] = RJ[i];
		//adults
		TotalSimSA[i][CurrentDay] = SA[i];
		TotalSimIA[i][CurrentDay] = IA[i];
		TotalSimRA[i][CurrentDay] = RA[i];
		//groups
		TotalSimJuv[i][CurrentDay] = SJ[i] + IJ[i] + RJ[i];
		TotalSimAdu[i][CurrentDay] = SA[i] + IA[i] + RA[i];
		//THE BELOW INCLUDES NEWBORNS SOLELY FOR DEBUGGING PURPOSES!!! TO BE ERASED LATER...
		TotalSimRats[i][CurrentDay] = SN[i] + IN[i] + RN[i] + SJ[i] + IJ[i] + RJ[i] + SA[i] + IA[i] + RA[i];
		
		//Instantenous R0 calculated using traditional formula
		R0OS[i][CurrentDay] = BA[i] * (SJ[i] + IJ[i] + RJ[i] + SA[i] + IA[i] + RA[i]) / (((dJ+dA)/2.0) + g+Alpha[i]*(IA[i] + IJ[i] + RA[i] + RJ[i] + SA[i] + SJ[i] + SN[i] + IN[i] + RN[i]));
		//Instantaneous R0 calculated using the next generation method (see M'ca notebooks)
		R0NG[i][CurrentDay] = (BA[i] * (dN + g + MNJ + Alpha[i] * (IA[i] + IJ[i] + RA[i] + RJ[i] + SA[i] + SJ[i] + SN[i] + IN[i] + RN[i]))*
			(dJ*(IA[i] + RA[i] + SA[i]) + dA * (IJ[i] + RJ[i] + SJ[i]) +
			(IA[i] + IJ[i] + RA[i] + RJ[i] + SA[i] + SJ[i])*(g + MJA + Alpha[i] * (IA[i] + IJ[i] + RA[i] + RJ[i] + SA[i] + SJ[i] + SN[i] + IN[i] + RN[i]))) +
			b[i] * MJA*MNJ*V + sqrt(pow(BA[i], 2)*pow(dN + g + MNJ +
				Alpha[i] * (IA[i] + IJ[i] + RA[i] + RJ[i] + SA[i] + SJ[i] + SN[i] + IN[i] + RN[i]), 2)*
				pow(dJ*(IA[i] + RA[i] + SA[i]) + dA * (IJ[i] + RJ[i] + SJ[i]) +
				(IA[i] + IJ[i] + RA[i] + RJ[i] + SA[i] + SJ[i])*(g + MJA + Alpha[i] * (IA[i] + IJ[i] + RA[i] + RJ[i] + SA[i] + SJ[i] + SN[i] + IN[i] + RN[i])), 2)
				+ 2 * b[i] * BA[i] * MNJ*(dN + g + MNJ + Alpha[i] * (IA[i] + IJ[i] + RA[i] + RJ[i] + SA[i] + SJ[i] + SN[i] + IN[i] + RN[i]))*
				(pow(MJA, 2)*(IA[i] + IJ[i] + RA[i] + RJ[i] + SA[i] + SJ[i]) +
					2 * (IA[i] + RA[i] + SA[i])*(dA + g + Alpha[i] * (IA[i] + IJ[i] + RA[i] + RJ[i] + SA[i] + SJ[i] + SN[i] + IN[i] + RN[i]))*
					(dJ + g + Alpha[i] * (IA[i] + IJ[i] + RA[i] + RJ[i] + SA[i] + SJ[i] + SN[i] + IN[i] + RN[i])) +
					MJA * (dJ*(IA[i] + RA[i] + SA[i]) + dA * (IJ[i] + RJ[i] + 2 * (IA[i] + RA[i] + SA[i]) + SJ[i]) +
					(IJ[i] + RJ[i] + 3 * (IA[i] + RA[i] + SA[i]) + SJ[i])*(g + Alpha[i] * (IA[i] + IJ[i] + RA[i] + RJ[i] + SA[i] + SJ[i] + SN[i] + IN[i] + RN[i]))))
				*V + pow(b[i], 2)*pow(MJA, 2)*pow(MNJ, 2)*pow(V, 2))) /
				(2.*(dA + g + Alpha[i] * (IA[i] + IJ[i] + RA[i] + RJ[i] + SA[i] + SJ[i] + SN[i] + IN[i] + RN[i]))*
			(dJ + g + MJA + Alpha[i] * (IA[i] + IJ[i] + RA[i] + RJ[i] + SA[i] + SJ[i] + SN[i] + IN[i] + RN[i]))*
					(dN + g + MNJ + Alpha[i] * (IA[i] + IJ[i] + RA[i] + RJ[i] + SA[i] + SJ[i] + SN[i] + IN[i] + RN[i])));
		//***********END of next gen method calcs****************************
		

		while (DaysElapsed > 0)//Fill in the skipped days with their values prior to this update. This is required because Gillespie can jump over days and then potentially miss days for which data is available preventing later comparison
		{
			//newborns
			TotalSimSN[i][T] = TotalSimSN[i][LastRecordDay];
			TotalSimIN[i][T] = TotalSimIN[i][LastRecordDay];
			TotalSimRN[i][T] = TotalSimRN[i][LastRecordDay];
			//juveniles
			TotalSimSJ[i][T] = TotalSimSJ[i][LastRecordDay];
			TotalSimIJ[i][T] = TotalSimIJ[i][LastRecordDay];
			TotalSimRJ[i][T] = TotalSimRJ[i][LastRecordDay];
			//adults
			TotalSimSA[i][T] = TotalSimSA[i][LastRecordDay];
			TotalSimIA[i][T] = TotalSimIA[i][LastRecordDay];
			TotalSimRA[i][T] = TotalSimRA[i][LastRecordDay];
			//groups
			TotalSimJuv[i][T] = TotalSimSJ[i][T] + TotalSimIJ[i][T] + TotalSimRJ[i][T];
			TotalSimAdu[i][T] = TotalSimSA[i][T] + TotalSimIA[i][T] + TotalSimRA[i][T];
			//THE BELOW INCLUDES NEWBORNS SOLELY FOR DEBUGGING PURPOSES!!! TO BE ERASED LATER...
			TotalSimRats[i][T] = TotalSimSN[i][T] + TotalSimIN[i][T] + TotalSimRN[i][T] + TotalSimSJ[i][T] + TotalSimIJ[i][T] + TotalSimRJ[i][T] + TotalSimSA[i][T] + TotalSimIA[i][T] + TotalSimRA[i][T];
			//R0
			R0NG[i][T] = R0NG[i][LastRecordDay];
			R0OS[i][T] = R0OS[i][LastRecordDay];

			T++;
			DaysElapsed--;
		}
	}
}


void SummaryStats()
{
	int k, FreqDataDays[10], AbundDataDays[10];
	double DT[10], DJ[10], DA[10], DSJ[10], DIJ[10], DRJ[10], DSA[10], DIA[10], DRA[10], TotalAbundDist, TotalFreqDist, TotDJDist, TotDADist, TotDSJDist, TotDIJDist, TotDRJDist, TotDSADist, TotDIADist, TotDRADist, ThreshEco, ThreshEpi;
	int DataPoints;

	//NOTE: I THINK THIS SHOULD BE MODIFIED TO ONLY OUTPUT DATA FROM DAY 365 (start day) ONWARDS!

	for (i = 0; i < Sites; i++)
	{
		if (TimeToRun > 10000)//Throw an error if data collection period exceeds local memory allocation of 10000 days
		{
			cout << "WARNING: MEMORY ALLOCATION EXCEEDED!!!\n";
			cin >> Test;
		}

		//GO TRAPPING DURING ONLY THOSE DAYS FOR WHICH DATA IS AVAILABLE
		//go trapping in the current time. Note that this could be a sample from a time warp into future. Thus, we backfill gaps linearly below in the while loop
		//We know how many nights Fichet-Calvet sampled
		//We assume all classes are equally likely to be trapped (infection classes and juvenile and adult. Newborns cannot be caught)
		//We assume each class is captured in proportion to its frequency within the pop
		//First calculate frequencies of types
		DataPoints = 0;
		for (k = StartTime; k < TimeToRun; k++)//calculate summaries throughout the time period in the data set
		{
			if (TrapNights[i][k] > 0)//if data is available for this day
			{
				//USE ANDREWS SAMPLING MODEL ON BOARD. INCLUDE DENSITY AND TRAPPING RATE. FREQUNCY IS NOT NECESSARY USING ANDREW's MODEL
				//Sampling model: dP_i/dt=-r*N_i where P is the probability of a trap being empty, r is the trapping rate and N_i is the density of class i
				//Solution to sampling model: P[Trap fails to catch a rat in a single day]=Exp[-r*N_i]
				//Thus, the probability of a trap catching a rat in a single day is equal to: 1-Exp[-r*N_i]

				//cout << TrapNights[k] << "," << 1.0 - exp(-TrapRate * TotalSimSN[k]) << "," << 1.0 - exp(-TrapRate * TotalSimSJ[k]) << "," << 1.0 - exp(-TrapRate * TotalSimSA[k]) << "\n";
				//cin >> Test;

				SNCaptured[i][k] = Bin_Rand(TrapNights[i][k], 1.0 - exp(-TrapRateN * TotalSimSN[i][k]));
				if (SNCaptured[i][k] > TotalSimSN[i][k])
				{
					SNCaptured[i][k] = TotalSimSN[i][k];
				}
				SJCaptured[i][k] = Bin_Rand(TrapNights[i][k], 1.0 - exp(-TrapRateJ * TotalSimSJ[i][k]));
				if (SJCaptured[i][k] > TotalSimSJ[i][k])
				{
					SJCaptured[i][k] = TotalSimSJ[i][k];
				}
				SACaptured[i][k] = Bin_Rand(TrapNights[i][k], 1.0 - exp(-TrapRateA * TotalSimSA[i][k]));
				if (SACaptured[i][k] > TotalSimSA[i][k])
				{
					SACaptured[i][k] = TotalSimSA[i][k];
				}

				INCaptured[i][k] = Bin_Rand(TrapNights[i][k], 1.0 - exp(-TrapRateN * TotalSimIN[i][k]));
				if (INCaptured[i][k] > TotalSimIN[i][k])
				{
					INCaptured[i][k] = TotalSimIN[i][k];
				}
				IJCaptured[i][k] = Bin_Rand(TrapNights[i][k], 1.0 - exp(-TrapRateJ * TotalSimIJ[i][k]));
				if (IJCaptured[i][k] > TotalSimIJ[i][k])
				{
					IJCaptured[i][k] = TotalSimIJ[i][k];
				}
				IACaptured[i][k] = Bin_Rand(TrapNights[i][k], 1.0 - exp(-TrapRateA * TotalSimIA[i][k]));
				if (IACaptured[i][k] > TotalSimIA[i][k])
				{
					IACaptured[i][k] = TotalSimIA[i][k];
				}

				RNCaptured[i][k] = Bin_Rand(TrapNights[i][k], 1.0 - exp(-TrapRateN * TotalSimRN[i][k]));
				if (RNCaptured[i][k] > TotalSimRN[i][k])
				{
					RNCaptured[i][k] = TotalSimRN[i][k];
				}
				RJCaptured[i][k] = Bin_Rand(TrapNights[i][k], 1.0 - exp(-TrapRateJ * TotalSimRJ[i][k]));
				if (RJCaptured[i][k] > TotalSimRJ[i][k])
				{
					RJCaptured[i][k] = TotalSimRJ[i][k];
				}
				RACaptured[i][k] = Bin_Rand(TrapNights[i][k], 1.0 - exp(-TrapRateA * TotalSimRA[i][k]));
				if (RACaptured[i][k] > TotalSimRA[i][k])
				{
					RACaptured[i][k] = TotalSimRA[i][k];
				}

				TotalJuvCap[i][k] = SJCaptured[i][k] + IJCaptured[i][k] + RJCaptured[i][k];
				TotalAduCap[i][k] = SACaptured[i][k] + IACaptured[i][k] + RACaptured[i][k];

				DataPoints++;
			}
		}
		//DONE TRAPPING



		DT[i] = 0;
		DJ[i] = 0;
		DA[i] = 0;
		DSJ[i] = 0;
		DIJ[i] = 0;
		DRJ[i] = 0;
		DSA[i] = 0;
		DIA[i] = 0;
		DRA[i] = 0;
		TotalAbundDist = 0;
		TotalFreqDist = 0;

		FreqDataDays[i] = 0;
		for (k = StartTime; k < TimeToRun; k++)//calculate summaries throughout the time period of data collection
		{
			//if (TotalSJ[i][k] > 0 && TotalIJ[i][k] > 0 && TotalRJ[i][k] > 0 && TotalSA[i][k] > 0 && TotalIA[i][k] > 0 && TotalRA[i][k] > 0)//If frequency data for infecion classes is available for this day, include it in the summary statistics
			//{
			if (TotalJuvCap[i][k] > 0 && TotalAduCap[i][k] > 0)
			{
				//Corrected 7/31/19 
				DSJ[i] = DSJ[i] + fabs(TotalSJ[i][k] - (SJCaptured[i][k] / (1.0*TotalJuvCap[i][k]))); //REMEMBER: TotalSJ is a frequency and is imported as such from the real and simulated data!
				DIJ[i] = DIJ[i] + fabs(TotalIJ[i][k] - (IJCaptured[i][k] / (1.0*TotalJuvCap[i][k])));
				DRJ[i] = DRJ[i] + fabs(TotalRJ[i][k] - (RJCaptured[i][k] / (1.0*TotalJuvCap[i][k])));
				DSA[i] = DSA[i] + fabs(TotalSA[i][k] - (SACaptured[i][k] / (1.0*TotalAduCap[i][k])));
				DIA[i] = DIA[i] + fabs(TotalIA[i][k] - (IACaptured[i][k] / (1.0*TotalAduCap[i][k])));
				DRA[i] = DRA[i] + fabs(TotalRA[i][k] - (RACaptured[i][k] / (1.0*TotalAduCap[i][k])));
				FreqDataDays[i]++;
			}
			//}
		}
		AbundDataDays[i] = 0;
		for (k = StartTime; k < TimeToRun; k++)//calculate summaries throughout the time period of data collection
		{
			if (TotalJuv[i][k] >= 0 && TotalAdu[i][k] >= 0)//If abundance data is available for this day, include it in the summary statistics
			{
				//cout << TotalJuv[k] << "," << TotalJuvCap[k] <<","<< TotalAdu[k] <Biol< "," << TotalAduCap[k] << "\n";
				//cin >> Test;
				//The code in the following block just calculates average distance between data and simulation across all time points for which data is available
				//EDITED JUN 28 TO MAKE PROPORTIONAL TO OBSERVED DENSITIES
				DJ[i] = DJ[i] + fabs(TotalJuv[i][k] - TotalJuvCap[i][k]) / (1.0*TotalJuv[i][k]);
				DA[i] = DA[i] + fabs(TotalAdu[i][k] - TotalAduCap[i][k]) / (1.0*TotalAdu[i][k]);
				DT[i] = DT[i] + fabs(TotalJuv[i][k] + TotalAdu[i][k] - TotalJuvCap[i][k] - TotalAduCap[i][k]) / (TotalJuv[i][k] + TotalAdu[i][k]);//ignores age strctr
				AbundDataDays[i]++;
			}

		}

		//Distance Summaries
		DT[i] = DT[i] / (1.0*AbundDataDays[i]);
		DJ[i] = DJ[i] / (1.0*AbundDataDays[i]);
		DA[i] = DA[i] / (1.0*AbundDataDays[i]);
		DSJ[i] = DSJ[i] / (1.0*FreqDataDays[i]);
		DIJ[i] = DIJ[i] / (1.0*FreqDataDays[i]);
		DRJ[i] = DRJ[i] / (1.0*FreqDataDays[i]);
		DSA[i] = DSA[i] / (1.0*FreqDataDays[i]);
		DIA[i] = DIA[i] / (1.0*FreqDataDays[i]);
		DRA[i] = DRA[i] / (1.0*FreqDataDays[i]);

		//Summarize R0
		MeanR0NG[i] = 0;
		for (j = StartTime; j < TimeToRun; j++)
		{
			MeanR0NG[i] = MeanR0NG[i] + R0NG[i][j] / (1.0*(TimeToRun - StartTime));
		}


	}//end of sites loop

	/*
	for (i = 0; i < Sites; i++)
	{
		cout << DJ[i] << "," << DA[i] << "," << DSJ[i] << "," << DIJ[i] << "," << DRJ[i] << "," << DSA[i] << "," << DIA[i] << "," << DRA[i] << "\n";
	}
	cin >> Test;
	*/


	//Calculate total distance between simulation and data weighting total abundance differently than differences in frequencies
	TotalAbundDist = 0;
	TotalFreqDist = 0;
	TotDJDist = 0;
	TotDADist = 0;
	TotDSJDist = 0;
	TotDIJDist = 0;
	TotDRJDist = 0;
	TotDSADist = 0;
	TotDIADist = 0;
	TotDRADist = 0;
	for (i = 0; i < Sites; i++)
	{
		TotalAbundDist = TotalAbundDist + DT[i];
		TotalFreqDist = TotalFreqDist + (DRJ[i] + DIJ[i] + DRA[i] + DIA[i]);//Changed August 27 to only count R and I since S is redundant in frequency space...
		TotDJDist = TotDJDist + DJ[i];
		TotDADist = TotDADist + DA[i];
		TotDSJDist = TotDSJDist + DSJ[i];
		TotDIJDist = TotDIJDist + DIJ[i];
		TotDRJDist = TotDRJDist + DRJ[i];
		TotDSADist = TotDSADist + DSA[i];
		TotDIADist = TotDIADist + DIA[i];
		TotDRADist = TotDRADist + DRA[i];

	}
	TotalAbundDist = TotalAbundDist / (1.0*Sites);
	TotalFreqDist = TotalFreqDist / (1.0*Sites);
	TotDJDist = TotDJDist / (1.0*Sites);
	TotDADist = TotDADist / (1.0*Sites);
	TotDSJDist = TotDSJDist / (1.0*Sites);
	TotDIJDist = TotDIJDist / (1.0*Sites);
	TotDRJDist = TotDRJDist / (1.0*Sites);
	TotDSADist = TotDSADist / (1.0*Sites);
	TotDIADist = TotDIADist / (1.0*Sites);
	TotDRADist = TotDRADist / (1.0*Sites);
	//cout << DJ << "," << DA << "," << DSJ << "," << DIJ << "," << DRJ << "," << DSA << "," << DIA << "," << DRA << "\n";
	//cin >> Test;


	//out_ErrorLog[0] << TotalAbundDist << "," << TotalFreqDist << "\n";
	//out_ErrorLog[1] << TotDJDist << "," << TotDADist << "," << TotDIJDist << "," << TotDRJDist << "," << TotDIADist << "," << TotDRADist << "\n";

	if (TotalAbundDist < 0.25&&TotalFreqDist < 0.5)//was 0.25 and 0.5 on first cluster run (0.225 and 0.45 works on cluster)
	{
		//cout << "trials:" << Trial << " Hits:" << Hits << "\n";//output status of run
		//output hits to posterior file
		out_Posterior << bWindow << "," << bTau << "," << bThresh << "," << bMin[0] << "," << bMin[1] << "," << bMax << "," << MNJ << "," << MJA << "," << Alpha[0] << "," << Alpha[1] << "," << dN << "," << dJ << "," << dA << "," << BA[0] << "," << BA[1] << "," << BM << "," << BN << "," << g << "," << V << "," << M << "," << mj << "," << ma << "," << TrapRateJ << "," << TrapRateA << "," << Lambda << "," << MeanR0NG[0] << "," << MeanR0NG[1] << "\n";
		out_Posterior.flush();
		//output data and simulated data that resulted in this hit (this is a debugging tool only)
		if (Hits < 10)
		{
			for (i = 0; i < Sites; i++)
			{
				for (k = 0; k < TimeToRun; k++)
				{
					//out_DSCompare << Day[k] << "," << Rain[k] << "," << LST[k] << "," << NDVI[k] << "," << TotalJuv[k] << "," << TotalSimJuv[k] << "," << TotalAdu[k] << "," << TotalSimAdu[k] << "\n";
					if (TrapNights[i][k] > 0)//output real data and simulated capture data if real data exists for a time point
					{
						//out_DSCompare << Day[k] << "," << Rain[k] << "," << LST[k] << "," << NDVI[k] << "," << TotalJuv[k] << "," << TotalSimJuv[k] << "," << TotalAdu[k] << "," << TotalSimAdu[k] << "," << TotalIJ[k] / (1.0*TotalJuv[k] + .0000001) << "," << TotalSimIJ[k] / (1.0*TotalSimJuv[k] + .0000001) << "," << TotalRJ[k] / (1.0*TotalJuv[k] + .0000001) << "," << TotalSimRJ[k] / (1.0*TotalSimJuv[k] + .0000001) << "," << TotalIA[k] / (1.0*TotalAdu[k] + .0000001) << "," << TotalSimIA[k] / (1.0*TotalSimAdu[k] + .0000001) << "," << TotalRA[k] / (1.0*TotalAdu[k] + .0000001) << "," << TotalSimRA[k] / (1.0*TotalSimAdu[k] + .0000001) << "\n";
						//out_DSCompare << Day[k] << "," << Rain[k] << "," << LST[k] << "," << NDVI[k] << "," <<TotalJuv[k]<<","<<TotalSimJuv[k]<<"," << TotalAdu[k] << "," << TotalSimAdu[k] << "," << TotalSJ[k] << "," << TotalSimSJ[k] / (TotalSimSJ[k] + TotalSimIJ[k] + TotalSimRJ[k]) << "," << TotalIJ[k]  << "," << TotalSimIJ[k] / (TotalSimSJ[k] + TotalSimIJ[k] + TotalSimRJ[k]) << "," << TotalRJ[k]  << "," << TotalSimRJ[k] / (TotalSimSJ[k] + TotalSimIJ[k] + TotalSimRJ[k]) << "," << TotalSA[k] << "," << TotalSimSA[k] / (TotalSimSA[k] + TotalSimIA[k] + TotalSimRA[k]) << "," << TotalIA[k] << "," << TotalSimIA[k] / (TotalSimSA[k] + TotalSimIA[k] + TotalSimRA[k]) << "," << TotalRA[k]<< "," << TotalSimRA[k] / (TotalSimSA[k] + TotalSimIA[k] + TotalSimRA[k]) << "\n";
						//out_DSCompare[i] << Day[i][k] << "," << Rain[i][k] << "," << LST[i][k] << "," << NDVI[i][k] << "," << TotalJuv[i][k] << "," << TotalSimJuv[i][k] << "," << TotalJuvCap[i][k] << "," << TotalAdu[i][k] << "," << TotalSimAdu[i][k] << "," << TotalAduCap[i][k] << "," << TotalSJ[i][k] << "," << TotalSimSJ[i][k] / (TotalSimSJ[i][k] + TotalSimIJ[i][k] + TotalSimRJ[i][k]) << "," << (SJCaptured[i][k] / (1.0*TotalJuvCap[i][k])) << "," << TotalIJ[i][k] << "," << TotalSimIJ[i][k] / (TotalSimSJ[i][k] + TotalSimIJ[i][k] + TotalSimRJ[i][k]) << "," << (IJCaptured[i][k] / (1.0* TotalJuvCap[i][k])) << "," << TotalRJ[i][k] << "," << TotalSimRJ[i][k] / (TotalSimSJ[i][k] + TotalSimIJ[i][k] + TotalSimRJ[i][k]) << "," << (RJCaptured[i][k] / (1.0*TotalJuvCap[i][k])) << "," << TotalSA[i][k] << "," << TotalSimSA[i][k] / (TotalSimSA[i][k] + TotalSimIA[i][k] + TotalSimRA[i][k]) << "," << (SACaptured[i][k] / (1.0*TotalAduCap[i][k])) << "," << TotalIA[i][k] << "," << TotalSimIA[i][k] / (TotalSimSA[i][k] + TotalSimIA[i][k] + TotalSimRA[i][k]) << "," << (IACaptured[i][k] / (1.0* TotalAduCap[i][k])) << "," << TotalRA[i][k] << "," << TotalSimRA[i][k] / (TotalSimSA[i][k] + TotalSimIA[i][k] + TotalSimRA[i][k]) << "," << (RACaptured[i][k] / (1.0* TotalAduCap[i][k])) << "\n";
					}
					/*
					else//otherwise output only simulated data
					{
						//out_DSCompare << Day[k] << "," << Rain[k] << "," << LST[k] << "," << NDVI[k] << ",," << TotalSimJuv[k] << ",," << TotalSimAdu[k] << ",," << TotalSimIJ[k] / (1.0*TotalSimJuv[k] + .0000001) << ",," << TotalSimRJ[k] / (1.0*TotalSimJuv[k] + .0000001) << ",," << TotalSimIA[k] / (1.0*TotalSimAdu[k] + .0000001) << ",," << TotalSimRA[k] / (1.0*TotalSimAdu[k] + .0000001) << "\n";
						out_DSCompare[i] << Day[i][k] << "," << Rain[i][k] << "," << LST[i][k] << "," << NDVI[i][k] << ",," << TotalSimJuv[i][k] << ",,," << TotalSimAdu[i][k] << ",,," << TotalSimSJ[i][k] / (TotalSimSJ[i][k] + TotalSimIJ[i][k] + TotalSimRJ[i][k]) << ",,," << TotalSimIJ[i][k] / (TotalSimSJ[i][k] + TotalSimIJ[i][k] + TotalSimRJ[i][k]) << ",,," << TotalSimRJ[i][k] / (TotalSimSJ[i][k] + TotalSimIJ[i][k] + TotalSimRJ[i][k]) << ",,," << TotalSimSA[i][k] / (TotalSimSA[i][k] + TotalSimIA[i][k] + TotalSimRA[i][k]) << ",,," << TotalSimIA[i][k] / (TotalSimSA[i][k] + TotalSimIA[i][k] + TotalSimRA[i][k]) << ",,," << TotalSimRA[i][k] / (TotalSimSA[i][k] + TotalSimIA[i][k] + TotalSimRA[i][k]) << "\n";
					}*/
				}

			}
		}
		Hits++;
	}

}
void Prior()
{
	double k, Mode, SetbMinHomo, bMinMaster, BetaMaster;

	//SetbMinHomo = Big_Rand(0.0, 0.2);//use this if you want to force the minimal birth rate to be static across space
	//First draw spatially variable parameters
	for (i = 0; i < Sites; i++)
	{
		//draw strength of density dependence from a uniform distribution
		Alpha[i] = Big_Rand(0.000008, 0.000035);// Big_Rand(0.000008, 0.00002); //Gam_Rand(k, Mode / (k - 1.0));
	}


	bTau = 0;// Rand_Int(1, 60);// Rand_Int(1, 60);// 0;//Currently no lag for linear rianpolation Rand_Int(1, 60);
	bWindow = Rand_Int(30, 365);//Define the window (in days) over which rainfall is averaged to determine the birth rate of mice
	bThresh = 0;//Not in use

	//Draw maximum feasible birth rate from a gamma distribution (note the name implies this is a minimal birth rate, but this is just a hangover from long ago)
	k = 50.0;
	Mode = 0.4296/2.0;//the mode is the mean rate of offspring production from Kyle. 10.74 pups every 25 days. 
	bMinMaster = Gam_Rand(k, Mode / (k - 1.0));
	bMin[0] = bMinMaster;// Gam_Rand(k, Mode / (k - 1.0)) / 2.0;//It is divided by two to account for only 1/2 the population being female
	bMin[1] = bMinMaster;// Gam_Rand(k, Mode / (k - 1.0)) / 2.0;//It is divided by two to account for only 1/2 the population being female


	//Draw sensitvity of birth rate to rainfall from a uniform distribution
	bMax = Big_Rand(0, 5.0); ;// Big_Rand(0, 1.0);// This parameter defines the sensitvity of birth rate to rainfall. Low values indicate a steady birth rate of bMin year round; large values indicate rainfall is required to approach bMin

	//Draw death rates from Gamma distributions with Shape k and mode M
	k = 20.0;
	Mode = 0.001;
	dN = Gam_Rand(k, Mode / (k - 1.0)); //DI Death rate for newborns
	k = 5.0;
	Mode = 0.003;
	dJ = Gam_Rand(k, Mode / (k - 1.0)); //DI Death rate for juveniles
	k = 5.0;
	Mode = 0.005;
	dA = Gam_Rand(k, Mode / (k - 1.0)); //DI Death rate for adults

	//Draw newborn to juvenile maturation rate from a gamma distribution
	k = 30.0;
	Mode = 0.05;
	MNJ = Gam_Rand(k, Mode / (k - 1.0));

	//Draw juvenile to adult maturation rate from a gamma distribution
	k = 20.0;
	Mode = 0.009;
	MJA = Gam_Rand(k, Mode / (k - 1.0));

	//Draw juvenile movement rate among populations from an exponential
	k = 1.0;
	Mode = 0.001;
	mj = Gam_Rand(k, Mode);//for exponential Mode is not actually the mode (obviously)

	//Draw adult immigration rate among populations from an exponential
	k = 1.0;
	Mode = 0.001;
	ma = Gam_Rand(k, Mode);//for exponential Mode is not actually the mode (obviously)

	//Draw trapping rates from gamma distributions
	k = 5.0;
	Mode = 0.000035;
	TrapRateN = 0;								//Newborns are not trapped
	TrapRateJ = Gam_Rand(k, Mode / (k - 1.0)); 	//Juveniles
	TrapRateA = Gam_Rand(k, Mode / (k - 1.0));	//Adults

	//Draw spatially uniform transmission rate from a uniform
	BetaMaster = Big_Rand(0.00007, 0.0003);//was 0.00007 to 0.00020 in original CRC run
	BA[0] = BetaMaster; //Free-living to free-living transmission
	BA[1] = BetaMaster; //Free-living to free-living transmission

	//Draw recovery rate from a gamma distribution
	k = 20.0;
	Mode = 0.0476;
	g = Gam_Rand(k, Mode / (k - 1.0));

	BM = 0;// Big_Rand(0, 0.02); //Adult to baby transmission
	BN = 0;// Big_Rand(0, 0.02); //Baby to baby transmsision

	//Draw vertical transmission parameters from a clipped exponential

	V = Gam_Rand(1.0, 0.05); //Vertical transmission of antibody
	if (V > 1.0)
	{
		V = 1.0;
	}


	//Draw vertical transmission of antbodies from a clipped exponential
	M = Gam_Rand(1.0, 0.05); //Vertical transmission of antibody
	if (M > 1.0)
	{
		M = 1.0;
	}

	//Draw density vs. frequency dependence from a uniform
	Lambda = 0;// Big_Rand(0, 1.0); //zero is density dependence, one is frequency dependence

	//Output the draw from the prior
	for (i = 0; i < Sites; i++)
	{
		if (Trial < 10000)//limit the output of priors to 10,000 to prevent mass storage of useless crap.
		{
			//out_Pars[i] << bWindow << "," << bTau << "," << bThresh << "," << bMinMaster << "," << bMax << "," << MNJ << "," << MJA << "," << Alpha[0] << "," << Alpha[1] << "," << dN << "," << dJ << "," << dA << "," << BA[0] << "," << BA[1] << "," << BM << "," << BN << "," << g << "," << V << "," << M << "," << mj << "," << ma << "," << TrapRateJ << "," << TrapRateA << "," << Lambda << "\n";
		}
	}
}



void IngestData()
{
	int Row, Column, k;


	for (k = 0; k < Sites; k++)
	{

		i = 0;
		Row = 0;
		Column = 0;

		while (in_Data[k] >> Data[k][Row][Column])
		{
			i++;
			Column = i % 16;// i % 9;//This modulus operator indicates when to switch to a new row and thus represents the number of columns (e.g., 2 indicates 2 columns)
			//cout << Data[Row][Column] << "\n";
			//cin >> Test;
			if (i > 0 && Column == 0)
			{
				Row++;
			}

		}


		TimeToRun = Row - 1;
		//cout << "TIMETORUN:" << TimeToRun << "\n";

		//Transfer the data into individual arrays for each variable
		//These are used only in SummaryStats() for comparison of simulated and observed data
		for (i = 0; i <= TimeToRun; i++)
		{
			Day[k][i] = i + 1;
			Rain[k][i] = Data[k][i][0];
			LST[k][i] = Data[k][i][1];// Data[TimeChunk][0];//record the rainfall at this particular point in time
			NDVI[k][i] = Data[k][i][2];// Data[TimeChunk][0];//record the rainfall at this particular point in time
			TotalSJ[k][i] = Data[k][i][3];//;//record the total number of trapped Mas at this point in time
			TotalIJ[k][i] = Data[k][i][4];//;//record the total number of trapped Mas at this point in time
			TotalRJ[k][i] = Data[k][i][5];//;//record the total number of trapped Mas at this point in time
			TotalSA[k][i] = Data[k][i][6];//;//record the total number of trapped Mas at this point in time
			TotalIA[k][i] = Data[k][i][7];//;//record the total number of trapped Mas at this point in time
			TotalRA[k][i] = Data[k][i][8];//;//record the total number of trapped Mas at this point in time
			TotalJuv[k][i] = Data[k][i][9]; //TotalSJ[i] + TotalIJ[i] + TotalRJ[i];
			TotalAdu[k][i] = Data[k][i][10]; //TotalSA[i] + TotalIA[i] + TotalRA[i];
			TotalRats[k][i] = TotalJuv[k][i] + TotalAdu[k][i]; //TotalSJ[i] + TotalIJ[i] + TotalRJ[ i] + TotalSA[i] + TotalIA[i] + TotalRA[i];
			JPCR[k][i] = Data[k][i][11];
			APCR[k][i] = Data[k][i][12];
			JSero[k][i] = Data[k][i][13];
			ASero[k][i] = Data[k][i][14];
			TrapNights[k][i] = Data[k][i][15];
			//cout << i<<","<<TotalRats[i] << "\n";
		}
		/*
		for (i = 0; i < TimeToRun; i++)
		{
			cout << Data[k][i][0] << "," << Data[k][i][1] << "," << Data[k][i][2] << "," << Data[k][i][3] << "," << Data[k][i][4] << "," << Data[k][i][5] << "," << Data[k][i][6] << "," << Data[k][i][7] << "," << Data[k][i][8] << "," << Data[k][i][9] << "," << Data[k][i][10] << "," << Data[k][i][11] << "," << Data[k][i][12] << "," << Data[k][i][13] << "," << Data[k][i][14] << "," << Data[k][i][15] << "\n";
		}
		cout << "LINEAR DATA FOLLOWS\n";
		*/
		for (i = 0; i < TimeToRun; i++)
		{
			//cout << Day[k][i] << "," << TotalJuv[k][i] << "," << TotalAdu[k][i] << "," << TrapNights[k][i] << "\n";
		}
		//cout << "Data ingestion look ok?\n";
		//cin >> Test;
	}
}

void Initialize()
{
	
	CurrentTime = StartTime;//Start the simulation 1 year in to allow back averaging of rainfall for a period of up to one year
	if (bWindow + bTau > 365)//If the priors are set to look at too much history, throw an error and wait for acknowledgement
	{
		cout << "WARNING: You have set your priors with too much back averaging and are overunning the rainfall data\n";
		cin >> Test;
	}
	if (bWindow + bTau >365)//If the priors are set to look at too much history, throw an error and wait for acknowledgement
	{
		cout << "WARNING: You have set your priors with too much back averaging and are overunning the rainfall data\n";
		cin >> Test;
	}

	for (i = 0; i < Sites; i++)
	{
		SN[i] = Rand_Int(75, 125); //Set initial numbers of newborns
		SJ[i] = Rand_Int(75, 125); //Set initial numbers of juveniles or sub adults
		SA[i] = Rand_Int(75, 125); //Set initial number of adults
		IN[i] = Rand_Int(75, 125);
		IJ[i] = Rand_Int(75, 125);
		IA[i] = Rand_Int(75, 125);
		RN[i] = Rand_Int(75, 125);
		RJ[i] = Rand_Int(75, 125);
		RA[i] = Rand_Int(75, 125);
	}
}



void IntegrateBinaryEnvironment()//No lag, just continuous back windowing, linear density dependence, and linear maturation rate
{
	int k;

	for (k = 0; k < Sites; k++)
	{
		//Calculate the cumulative rainfall over the previous "Window" days and use it to calcuate current birth rates
		AverageRain = 0;
		for (i = 0; i < bWindow; i++)//calculate cumulative rainfall in the "Window" days previous as a rolling window
		{
			AverageRain = AverageRain + Data[k][CurrentDay - i - bTau][0] / (1.0*bWindow);
		}
		//NOW USING Mckalis Menton Rainfall
		b[k] = (bMin[k]*AverageRain)/(bMax+AverageRain);
	

		if (b[k] < 0)//if a negative birth rate is generated throw an error and await acknowledgement
		{
			cout << "WARNING: Negative birth rate\n";
			cout << b[k] << "\n";
			cout << AverageRain << "," << bMin << "," << bMax << "\n";
			cin >> Test;
		}
	}
	
	/*
	//Calculate the cumulative rainfall over the previous "Window" days and use it to calculate the current intensity of competition
	AverageRain = 0;
	for (i = 0; i < AlphaWindow; i++)//calculate cumulative rainfall in the "Window" days previous as a rolling window
	{
		AverageRain = AverageRain + Data[TimeChunk - i][0] / (1.0*AlphaWindow);
	}
	//How strong is competition given the current food availability? 
	if (AverageRain > AlphaThresh)// If previous rainfall above a threshold value, competition takes "maximal" value
	{
		Alpha = AlphaMax;
	}
	if (AverageRain <= AlphaThresh)// If previous rainfall below a threshold value, competition takes "minimal" value
	{
		Alpha = AlphaMin;
	}
	*/

}