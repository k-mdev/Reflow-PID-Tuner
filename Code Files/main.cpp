#include <iostream>
#include <fstream>
#include <vector>
#include <algorithm>
#include <string>
#include "GNUPLOT_H.h";
#include <Windows.h>
#include <time.h>
#include <random>

using namespace std;

void setPID(double cTime, double PERIOD, double* tempSetPoint)
{
	if (cTime < (90000 / PERIOD))
	{
		if (cTime <= 0)
		{
			*tempSetPoint = 26;
		}
		*tempSetPoint += (0.001 * PERIOD);
	}
	else if (cTime < (180000 / PERIOD))
	{
		*tempSetPoint += (0.0005 * PERIOD);
	}
	else if (cTime < (225000 / PERIOD))
	{
		*tempSetPoint += 0.0014 * PERIOD;
	}
	else if (cTime < (255000 / PERIOD))
	{
		*tempSetPoint = 225;
	}
	else
	{
		*tempSetPoint -= 0.004 * PERIOD;
		if (*tempSetPoint < 26)
		{
			*tempSetPoint = 26;
		}
	}
}

double fRand(double fMin, double fMax)
{
	return fMin + (rand() / (RAND_MAX / (fMax - fMin)));
}

double getTemp(double cOTime, double PERIOD, double* ovenTemp, double* PID_total)
{
	double tempNoise = fRand(0, (0.0005 * PERIOD));
	if (cOTime < (255000 / PERIOD))
	{
		*ovenTemp += (*PID_total) - tempNoise;
		if (*ovenTemp > 300) {
			*ovenTemp = 300;
		}
		else if (*ovenTemp < 25) {
			*ovenTemp = 25;
		}
	}
	if (cOTime > (255000 / PERIOD))
	{
		*ovenTemp += (*PID_total) + tempNoise;
		if (*ovenTemp < 25)
		{
			*ovenTemp = 25;
		}
	}

	return tempNoise;
}

void PID(double* tempError, double* tempSetPoint, double* ovenTemp, double* prevTempError,
	double* PID_p, double* PID_i, double* PID_d, double* PID_total, double kp, double ki, double kd, double PERIOD)
{
	*tempError = *tempSetPoint - *ovenTemp;
	double temp_difference = *tempError - *prevTempError;

	*PID_p = kp * *tempError;
	*PID_d = kd * ((temp_difference) / PERIOD);

	if (-3 < *tempError && *tempError < 3)
	{
		*PID_i = *PID_i + (ki * *tempError);
	}
	else
	{
		*PID_i = 0;
	}
	*PID_total = *PID_p + *PID_i + *PID_d;
	*prevTempError = *tempError;
}

vector<vector<double>> runPIDErrorSim(double PERIOD, double pGain, double iGain, double dGain)
{
	//SIMULATION VARIABLES
	double tempSetPoint = 26;		//Reflow curve temperature (time-based)
	double ovenTemp = 25;
	double tempError;		//Deviation from reflow curve
	double prevTempError = 0;		//Last second's deviation
	double PID_p = 0, PID_i = 0, PID_d = 0, PID_total = 0;
	double noise = 0;

	double cookTime = 0;
	vector<double> error;
	vector<double> ovTemp;
	vector<double> idealTemp;
	vector<vector<double>> simData;

	while (cookTime < 360)
	{
		noise = getTemp(cookTime, PERIOD, &ovenTemp, &PID_total);
		setPID(cookTime, PERIOD, &tempSetPoint);
		PID(&tempError, &tempSetPoint, &ovenTemp, &prevTempError, &PID_p, &PID_i, &PID_d, &PID_total,
			pGain, iGain, dGain, PERIOD);

		error.push_back(tempError);
		ovTemp.push_back(ovenTemp);
		idealTemp.push_back(tempSetPoint);

		cookTime += (PERIOD / 1000); //SET TICK RATE HERE
	}

	simData.push_back(error);	//0 - Error Data
	simData.push_back(ovTemp);	///1 - Oven data
	simData.push_back(idealTemp); 	//2 - idealTemp chart

	return simData;
}

vector<vector<double>> generatePopulation(int numSims, double MUTATION, double PERIOD, double pGain, double iGain, double dGain) {
	
	//Vectors to store population data and PID values at injection
	vector<vector<double>> allSimPIDs;
	vector<vector<vector<double>>> allSimData;
	allSimPIDs.resize(numSims);
	allSimData.resize(numSims);

	//Initialize PID Variables
	double kp = 0;
	double ki = 0;
	double kd = 0;

	//Generate population of reflow sims with varying PID values
	for (int i = 0; i < numSims; i++) {
		//Mutate PID values for each reflow SIM
		if (i == 0) {
			kp = pGain;
			ki = iGain;
			kd = dGain;
		}
		else {
			kp = pGain * fRand(1 - MUTATION, 1 + MUTATION);
			ki = iGain * fRand(1 - MUTATION, 1 + MUTATION);
			kd = dGain * fRand(1 - MUTATION, 1 + MUTATION);
		}


		//Store PID values injected into simulation
		allSimPIDs.at(i).push_back(kp);
		allSimPIDs.at(i).push_back(ki);
		allSimPIDs.at(i).push_back(kd);

		//Store returned sim data into respective slot
		allSimData.at(i) = (runPIDErrorSim(PERIOD, kp, ki, kd));
	}

	//Store error reduction rating for each sim
	vector<double> errReductScores;

	//Calculate error reduction rating for sims
	for (int i = 0; i < numSims; i++) {

		vector<double> errSquared;
		double sumErrSquared = 0;

		//Square all errors
		for (int j = 0; j < allSimData.at(0).size(); j++) {
			errSquared.push_back(allSimData.at(i).at(0).at(j) * allSimData.at(i).at(0).at(j));
		}
		//Sum all squared errors
		for (double sqErr : errSquared) {
			sumErrSquared += sqErr;
		}
		//Store sum as score for sim #i
		errReductScores.push_back(sumErrSquared);
	}

	//Bulk data vector to store data on most error-reducing sim
	vector<vector<double>> minErrorSimData;

	//Get index of simulation with minimum error score
	int minErrorSimIndex = distance(errReductScores.begin(), min_element(errReductScores.begin(), errReductScores.end()));
	
	//Get error rating to return
	double errRating = errReductScores.at(minErrorSimIndex);
	errReductScores.clear();
	errReductScores.resize(0);
	errReductScores.push_back(errRating);
	
	//Get PID values of sim and its temp curves for graphing
	minErrorSimData.push_back(allSimPIDs.at(minErrorSimIndex)); //PID Gain Values
	minErrorSimData.push_back(errReductScores); //Error reduction rating
	minErrorSimData.push_back(allSimData.at(minErrorSimIndex).at(1)); //Oven temp curve
	minErrorSimData.push_back(allSimData.at(minErrorSimIndex).at(2)); //Ideal reflow temp curve

	//Return data for graphing and next population iteration
	return minErrorSimData;
	}

void graphCurves(vector<double> ovTemp, vector<double> idealTemp) {
	
	ofstream myfile;
	myfile.open("temperatures.txt");
	if (myfile.is_open()) {
	for (int i = 0; i < ovTemp.size(); i++) {
			myfile << i << "  " << ovTemp.at(i) << "\n";
		}
		myfile << "\n\n";

		for (int i = 0; i < idealTemp.size(); i++) {
			myfile << i << "  " << idealTemp.at(i) << "\n";
		}
	}
	myfile.close();

	Gnuplot p;

	//cin.ignore();
	Sleep(500);
}

int main()
{
	srand(time(NULL));

	int numSims = 30;
	double MUTATION = 0.2;
	double kp = 20; //0.58
	double ki = 20; //1.24
	double kd = 20; //7.22
	double PERIOD = 1000;

	vector<vector<double>> mostFitSim;

	cout << "Press enter to begin" << endl;
	cin.ignore();

	mostFitSim = generatePopulation(numSims, MUTATION, PERIOD, kp, ki, kd);
	graphCurves(mostFitSim.at(2), mostFitSim.at(3));

	while (mostFitSim.at(1).at(0) > 5) {

		kp = mostFitSim.at(0).at(0);
		ki = mostFitSim.at(0).at(1);
		kd = mostFitSim.at(0).at(2);

		cout << "Best Kp: " << kp << endl;
		cout << "Best Ki: " << ki << endl;
		cout << "Best Kd: " << kd << endl;

		mostFitSim = generatePopulation(numSims, MUTATION, PERIOD, kp, ki, kd);
		graphCurves(mostFitSim.at(2), mostFitSim.at(3));

	}

	cout << "FINAL PID:" << endl;
	cout << mostFitSim.at(0).at(0) << endl;
	cout << mostFitSim.at(0).at(1) << endl;
	cout << mostFitSim.at(0).at(2) << endl;

	return 0;
}

