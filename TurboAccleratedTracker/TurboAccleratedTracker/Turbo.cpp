#include <cstdlib>
#include <cstdio>
#include <iostream>
#include <cmath>
#include <vector>
#include <algorithm>
#include <random>
#include <numeric>
#include <chrono>
#include <time.h>

static double a[4] = { 2.50662823884, -18.61500062529, 41.39119773534, -25.44106049637 };
static double b[4] = { -8.47351093090, 23.08336743743, -21.06224101826, 3.13082909833 };
static double c[9] = { 0.3374754822726147, 0.9761690190917186, 0.1607979714918209,
0.0276438810333863, 0.0038405729373609, 0.0003951896511919,
0.0000321767881768, 0.0000002888167364, 0.0000003960315187 };

using namespace std;

void cholesky_deco(vector<vector<double> > &cor, int n)
{
	double cho_sum = 0;
	int i, j, k = 0;

	for (i = 0; i < n; i++)
	{
		for (j = i; j < n; j++)
		{
			for (cho_sum = cor[i][j], k = i - 1; k >= 0; k--)
			{
				cho_sum -= cor[i][k] * cor[j][k];
			}
			if (i == j) cor[i][i] = sqrt(cho_sum);
			else
			{
				cor[j][i] = cho_sum / cor[i][i];
				cor[i][j] = 0;
			}
		}
	}
}

int GetRequiredDimension(vector<vector<double> > &Table, int no_watchday, int no_asset)
{
	int i, j = 0;
	for (i = no_watchday; i >= 0; i--) {
		if (Table[0][i] == 0) j++;
	}
	return j*no_asset + 1;
}

double cndev(double u)
{
	double x, r;

	x = u - 0.5;
	if (fabs(x) < 0.42) {
		r = x*x;
		r = x*(((a[3] * r + a[2])*r + a[1])*r + a[0]) / ((((b[3] * r + b[2])*r + b[1])*r + b[0])*r + 1.0);
		return r;
	}

	r = u;
	if (x > 0.0) r = 1.0 - u;
	r = log(-log(r));
	r = c[0] + r*(c[1] + r*(c[2] + r*(c[3] + r*(c[4] + r*(c[5] + r*(c[6] + r*(c[7] + r*c[8])))))));

	if (x < 0.0) r = -r;

	return r;
}

void MoroGauss(int n_gauss, vector<double> &rv)
{
	for (int i = 0; i < n_gauss; i++) rv[i] = cndev(rv[i]);
}

void CholeskyTrans(vector<vector<double> > &cho_matrix, vector<vector<double> > &corrv, vector<double> &indrv, int dimension, int time)
{
	int m1, m2, m11 = 0;

	for (m11 = 0; m11 < time; m11++) {
		for (m1 = 0; m1 < dimension; m1++) {
			corrv[m1][m11] = 0;
			for (m2 = 0; m2 <= m1; m2++) {
				corrv[m1][m11] = corrv[m1][m11] + indrv[m11*dimension + m2] * cho_matrix[m1][m2];
			}
		}
	}
}

long _stdcall isKnockOut(long type, long no_observation, long double y, double r, double ParticipationRate,
	double put_k, double H, double bonuscouponRate,
	double* spot, double* spot2, double* vol, double* stk, double t, double FX0)
{
	long no_asset = 1;
	int no_period = no_observation + 1;  //總期數=觀察期總數+1
	int n_past = 0;
	long ko_flag = 0;

	vector<double> R_final(no_period + 1);
	vector<vector<double> > DailyArray(no_asset, vector<double>(no_period));  //FixingDate Price(no_observation*no_asset)
	vector<vector<double> > StockTable(no_asset, vector<double>(no_period + 1));  //for simulation use
	vector<vector<double> > ReturnTable(no_asset, vector<double>(no_period));     //for simulation use

																				  //DailyArray設定:由外部傳入資料
	for (int i = 0; i < no_asset; i++) {
		for (int j = 0; j < no_period; j++) {
			//DailyArray[j][i] = stk[i * (no_period)+j];
			DailyArray[i][j] = stk[i * (no_observation)+j];
		}
	}

	//StockTable設定:內部計算
	for (int i = 0; i < no_asset; i++) { StockTable[i][0] = spot[i]; }  //  //-1.the first period
	for (int j = 0; j < no_period; j++) {
		for (int i = 0; i < no_asset; i++) {
			StockTable[i][j + 1] = DailyArray[i][j];
		}
	}
	//ReturnTable設定:內部計算
	for (int i = 0; i < no_asset; i++) {
		for (int j = 0; j < no_period; j++) {
			ReturnTable[i][j] = DailyArray[i][j] / spot[i];
		}
	}

	//-2.the following period
	for (int j = 1; j <= no_observation; j++) {
		if (StockTable[0][j] > 0.0) {
			R_final[j] = StockTable[0][j] / StockTable[0][0];
			for (int m1 = 0; m1 < no_asset; m1++) {
				if (StockTable[m1][j] / StockTable[m1][0] < R_final[j]) {
					R_final[j] = StockTable[m1][j] / spot[m1];
				}
			}
			if (R_final[j] >= H) {  //R_final[j] - H <= err  => R_final[j] >= H
									//value1 = put_k*re_coupon[j] + put_k*(1 - exp(-y*(no_period - j)*dt));

				switch (type)
				{
				case 0:  //Normal
						 //value1 = put_k*(re_coupon[j]) + put_k*(1 - exp(-y*(no_period - j)*dt));
						 //value1 = put_k * re_coupon[j];  //提前出場配息
					ko_flag = 1;
					break;
					/*
					case 1:  //Floating
					value1 = (put_k*(re_coupon[j]) + put_k*(1 - exp(-y*(no_period - j)*dt)))*FXt;  //final_payoff幣別=USD
					break;

					case 2:  //Compo
					value1 = (put_k*(re_coupon[j]) + put_k*(1 - exp(-y*(no_period - j)*dt)))*FX0;  //final_payoff幣別=USD
					break;

					case 3:  //Quanto
					value1 = (put_k*(re_coupon[j]) + put_k*(1 - exp(-y*(no_period - j)*dt)))*FXT;  //final_payoff幣別=USD
					break;
					*/
				}

				return ko_flag;
			}
			/* 因為最後一期不配息, 所以此區塊不判斷
			else if (j == no_observation && R_final[j] < put_k) {  //R_final[j] - put_k > err => R_final[j] < put_k
			double value1 = 0;
			//value1 = R_final[j] - put_k;  //value1<0

			switch (type)
			{
			case 0:  //Normal
			value1 = -1 * max(put_k - R_final[j], 0.0);  //final_payoff幣別=HKD
			break;

			case 1:  //Floating
			value1 = -1 * max(put_k - R_final[j], 0.0)*FXt;  //final_payoff幣別=USD
			break;

			case 2:  //Compo
			value1 = -1 * max(FX0*put_k - FXt*R_final[j], 0.0);  //final_payoff幣別=USD
			break;

			case 3:  //Quanto
			value1 = -1 * max(put_k - R_final[j], 0.0)*FXT;  //final_payoff幣別=USD
			break;
			}

			return value1;
			}
			*/
			n_past++;
		}
		else
		{
			ko_flag = 0;
			return ko_flag;
		}
		//for (i = 0; i < no_asset; i++) { st_m[j][i] = StockTable[j][i]; }
	}
}

double _stdcall myTurboNote_Price(long type, long no_observation, long double y, double r, double ParticipationRate,
	double put_k, double H, double bonuscouponRate,
	double* spot, double* spot2, double* vol, double* stk, double t, double FX0)
{
	long seed = -1;  //亂數種子
	std::mt19937 generator(seed);
	std::uniform_real_distribution<double> unif(0.0, 1.0);
	//std::normal_distribution<double> norm(0.0, 1.0);
	//int i, j, k;
	double temp = -1;
	int temp1 = -2;
	long n_trial = 10000;
	int n_past = 0;
	int n_left = 0;
	int ko_flag = 0;
	int ko_period = 0;
	double err = 0.000001;
	double dt1 = 0.0;
	int currentPeriod = 0;
	int no_period = no_observation + 1;  //總期數=觀察期總數+1
	double FXt = 1;
	double FXT = FXt;
	double dt = 1.0 / no_period;
	long period_now = 0;

	long no_asset = 1;
	double coupon_sum = 0.0;
	double discount_time = 0.0;
	double discount_payoff = 0.0;
	double	payoff_sum = 0.0;
	double	payoff_average = 0.0;
	int required_dimension = 0;

	vector<int> count(3);  //count[0]=提前出場次數;count[1]=到期<K次數;count[2]=到期>=K次數
						   //int count[no_observation] = {0};

						   //Aarray Parameter
						   //vector<double> dt(no_observation);
	vector<double> R_final(no_period + 1);
	vector<double> div(no_asset);
	vector<double> final_payoff(n_trial);
	//vector<int> obDay(no_observation);   //觀察期間Vector
	vector<double> re_coupon(no_period);  //觀察日KO的利息補償
	vector<double> DriftRate(no_asset);

	for (int i = 0; i < no_observation; i++) { re_coupon[i] = bonuscouponRate * i; }  //第i期KO, coupon = couponRate * i * Notional
	re_coupon[no_observation] = 0.0;    //最後一期KO, coupon=0

										//Double Aarray Parameter
	vector<vector<double> > DailyArray(no_asset, vector<double>(no_period));  //FixingDate Price(no_observation*no_asset)
	vector<vector<double> > StockTable(no_asset, vector<double>(no_period + 1));  //for simulation use
	vector<vector<double> > ReturnTable(no_asset, vector<double>(no_period));     //for simulation use
	vector<vector<double> > corre(no_asset, vector<double>(no_asset));                 //for simulation use
	vector<vector<double> > corrv(no_asset, vector<double>(no_period + 1));       //for simulation use

																				  //int no_day = std::accumulate(obDay.begin(), obDay.end(), 0);  //契約總天數
																				  //計算dt[i]
																				  //for (i = 0; i < no_observation; i++) { dt[i] = static_cast<double>(observationDay[i]) / no_day; }
																				  //dividend設定
	for (int i = 0; i < no_asset; i++)
	{
		div[i] = 0.0;
		switch (type)
		{
		case 0:  //Normal
			DriftRate[i] = r - div[i];   //r_HKD-q
			break;
			/*
			case 1:  //Floating
			DriftRate[i] = r - div[i];   //r_HKD-q
			break;

			case 2:  //Compo
			DriftRate[i] = rf - div[i];   //r_USD-q
			vol[i] = sqrt(vol[i] * vol[i] + volFX * volFX + 2 * corrFX * vol[i] * volFX);
			break;

			case 3:  //Quanto
			DriftRate[i] = r - div[i] - corrFX * vol[i] * volFX;   //r_TWD-q-corr_FX*vol*volFX
			break;
			*/
		}
	}

	//DailyArray設定:由外部傳入資料
	for (int i = 0; i < no_asset; i++) {
		for (int j = 0; j < no_period; j++) {
			//DailyArray[j][i] = stk[i * (no_period)+j];
			DailyArray[i][j] = stk[i * (no_observation)+j];
		}
	}

	//StockTable設定:內部計算
	for (int i = 0; i < no_asset; i++) { StockTable[i][0] = spot[i]; }  //  //-1.the first period
	for (int j = 0; j < no_period; j++) {
		for (int i = 0; i < no_asset; i++) {
			StockTable[i][j + 1] = DailyArray[i][j];
		}
	}

	//ReturnTable設定:內部計算
	for (int i = 0; i < no_asset; i++) {
		for (int j = 0; j < no_period; j++) {
			ReturnTable[i][j] = DailyArray[i][j] / spot[i];
		}
	}
	//corre設定:內部計算
	for (int i = 0; i < no_asset; i++) {
		for (int j = 0; j < no_asset; j++) {
			//corre[i][j] = corr[i*no_asset + j];
			corre[i][j] = 1.0;
		}
	}

	//-2.the following period
	for (int j = 1; j <= no_period; j++) {
		if (StockTable[0][j] > 0.0) {
			R_final[j] = StockTable[0][j] / StockTable[0][0];
			for (int m1 = 0; m1 < no_asset; m1++) {
				R_final[j] = StockTable[m1][j] / spot[0] < R_final[j] ? StockTable[m1][j] / spot[m1] : R_final[j];
			}
			if (R_final[j] >= H && j != no_period) { //最後一期不能KO
				double value1 = 0;
				//value1 = put_k*re_coupon[j] + put_k*(1 - exp(-y*(no_period - j)*dt));

				switch (type)
				{
				case 0:  //Normal
						 //value1 = put_k*(re_coupon[j]) + put_k*(1 - exp(-y*(no_period - j)*dt));
					value1 = put_k * re_coupon[j];  //提前出場配息
					break;
					/*
					case 1:  //Floating
					value1 = (put_k*(re_coupon[j]) + put_k*(1 - exp(-y*(no_period - j)*dt)))*FXt;  //final_payoff幣別=USD
					break;

					case 2:  //Compo
					value1 = (put_k*(re_coupon[j]) + put_k*(1 - exp(-y*(no_period - j)*dt)))*FX0;  //final_payoff幣別=USD
					break;

					case 3:  //Quanto
					value1 = (put_k*(re_coupon[j]) + put_k*(1 - exp(-y*(no_period - j)*dt)))*FXT;  //final_payoff幣別=USD
					break;
					*/
				}

				return value1;
			}
			else if (j == no_period) {
				double value1 = 0;
				if (R_final[j] < put_k) {   //ST< K 就會賠錢
					switch (type)
					{
					case 0:  //Normal
						value1 = -exp(-r*(j - no_period)*dt)*(max(put_k - R_final[j], 0.0));  //final_payoff幣別=HKD
						break;
						/*
						case 1:  //Floating
						value1 = - exp(-rf*(currentPeriod - n_past)*dt)*(max(put_k - R_final[currentPeriod], 0.0))*FXt;  //final_payoff幣別=USD
						break;

						case 2:  //Compo
						value1 = - exp(-rf*(currentPeriod - n_past)*dt)*(max(FX0*put_k - FXt*R_final[currentPeriod], 0.0));  //final_payoff幣別=USD
						break;

						case 3:  //Quanto
						value1 = - exp(-rf*(currentPeriod - n_past)*dt)*(max(put_k - R_final[currentPeriod], 0.0))*FXT;  //final_payoff幣別=USD
						break;
						*/
					}

					//fprintf(fout, "%d : %.8lf\n", n_final, R_final[m11]);
				}
				else {
					//ST>=K =>拿100%+參與率*上漲幅度
					switch (type)
					{
					case 0:  //Normal
						value1 = exp(-r*(j - no_period)*dt)*ParticipationRate*max(R_final[j] - 1, 0.0);  //final_payoff幣別=HKD
						break;
						/*
						case 1:  //Floating
						value1= exp(-rf*(currentPeriod - n_past)*dt)*ParticipationRate*(R_final[currentPeriod] - spot[0])* FXt;  //final_payoff幣別=USD
						break;

						case 2:  //Compo
						value1 = exp(-rf*(currentPeriod - n_past)*dt)*ParticipationRate*((FXt*R_final[currentPeriod] - FX0*spot[0]));  //final_payoff幣別=USD
						break;

						case 3:  //Quanto
						value1 = exp(-rf*(currentPeriod - n_past)*dt)*ParticipationRate*(R_final[currentPeriod] - spot[0]) * FXT;  //final_payoff幣別=USD
						break;
						*/
					}
				}
				return value1;
			}
			n_past++;
		}
		else break;

		//for (i = 0; i < no_asset; i++) { st_m[j][i] = StockTable[j][i]; }
	}

	//將最近一日的股價當作最近一期的價格=>以便後續進行模擬試算
	for (int i = 0; i < no_asset; i++) {
		StockTable[i][n_past] = spot2[i];
	}

	/*
	//判斷目前是第幾期
	for (int i = 0; i < no_asset; i++) {
	for (int j = 0; j < no_period; j++) {
	if (StockTable[i][j + 1] == 0) {
	period_now = j+1;  //period_now=n_past+1
	break;
	}
	}
	}
	*/
	double min_dt1 = 1.0 / 250.0;
	n_left = no_period - n_past;
	dt1 = max(t - dt*(n_left - 1), dt);

	/*2-處理correlation table*/
	cholesky_deco(corre, no_asset);

	/*3-定義所需的Random Vector大小*/
	required_dimension = GetRequiredDimension(StockTable, no_observation, no_asset);
	vector<double> randomvector(required_dimension);  //randomvector

	if (ko_flag == 0)    //若尚未KO,則開始模擬程序
	{
		//--------------開始進行模擬-----------------------------------------//
		for (int k = 1; k <= n_trial; k++)
		{
			ko_flag = 0;
			//int m1 = 0;
			// m111 = 0;
			for (int i = n_past + 1; i <= no_period; i++) { R_final[i] = 100.0; }  //重設R_final

			for (int j = 0; j < required_dimension; j++) {
				//randomvector[j] = ran1(&seed);      //取均勻分配亂數
				randomvector[j] = unif(generator);    //梅森旋轉演算法:產生均勻分配亂數
													  //randomvector[j] = norm(generator);   //C++11內建產生常態分配亂數
			}
			MoroGauss(required_dimension, randomvector);  //Moro's取常態分配亂數
			CholeskyTrans(corre, corrv, randomvector, no_asset, n_left);    //generate correlated normal random vector

			if (ko_flag == 0)  //開始模擬後,若有KO,則停止模擬
			{
				//for (int m1 = 0; m1 < no_asset; m1++) {StockTable[m1][n_past] = StockTable[m1][0];}
				for (int m111 = n_past + 1; m111 <= no_period; m111++)
				{
					for (int m1 = 0; m1 < no_asset; m1++)  //m1=標的個數
					{
						StockTable[m1][m111] = StockTable[m1][m111 - 1] * exp((DriftRate[m1] - vol[m1] * vol[m1] / 2.0)*(dt)+vol[m1] * sqrt(dt)*corrv[m1][m111 - (n_past + 1)]);
						//判斷最小報酬標的及最小報酬
						R_final[m111] = StockTable[m1][m111] / spot[m1] < R_final[m111] ? StockTable[m1][m111] / spot[m1] : R_final[m111];
					}  //end of m1 loop
					   //KO判斷
					if (ko_flag == 0 && m111 <= no_observation) {  //尚未KO且期數<=觀察期總數,才進來判斷
						for (int m1 = 0; m1 < no_asset; m1++) {
							if (StockTable[m1][m111] / spot[m1] >= H)
							{
								ko_flag = 1;
								ko_period = m111 - n_past;
								break;
							}
						}
					}

					currentPeriod = m111;
				}  //end of m111 loop

				   //到期判斷(decide the final payoff)
				if (ko_flag == 0 && currentPeriod == no_period) {  //沒KO,到期  => 都拿100%
					if (R_final[currentPeriod] < put_k) {   //ST< K 就會賠錢
						switch (type)
						{
						case 0:  //Normal
							final_payoff[k - 1] = final_payoff[k - 1] - exp(-r*(currentPeriod - n_past)*dt)*(max(put_k - R_final[currentPeriod], 0.0));  //final_payoff幣別=HKD
							break;
							/*
							case 1:  //Floating
							final_payoff[k - 1] = final_payoff[k - 1] - exp(-rf*(currentPeriod - n_past)*dt)*(max(put_k - R_final[currentPeriod], 0.0))*FXt;  //final_payoff幣別=USD
							break;

							case 2:  //Compo
							final_payoff[k - 1] = final_payoff[k - 1] - exp(-rf*(currentPeriod - n_past)*dt)*(max(FX0*put_k - FXt*R_final[currentPeriod], 0.0));  //final_payoff幣別=USD
							break;

							case 3:  //Quanto
							final_payoff[k - 1] = final_payoff[k - 1] - exp(-rf*(currentPeriod - n_past)*dt)*(max(put_k - R_final[currentPeriod], 0.0))*FXT;  //final_payoff幣別=USD
							break;
							*/
						}

						discount_payoff = discount_payoff + final_payoff[k - 1];
						count[1] += 1;
						//fprintf(fout, "%d : %.8lf\n", n_final, R_final[m11]);
					}
					else {
						//ST>=K =>拿100%+參與率*上漲幅度
						switch (type)
						{
						case 0:  //Normal
							final_payoff[k - 1] = final_payoff[k - 1] + exp(-r*(currentPeriod - n_past)*dt)*ParticipationRate*max(R_final[currentPeriod] - 1, 0.0);  //final_payoff幣別=HKD
							break;
							/*
							case 1:  //Floating
							final_payoff[k - 1] = final_payoff[k - 1] + exp(-rf*(currentPeriod - n_past)*dt)*ParticipationRate*(R_final[currentPeriod] - spot[0])* FXt;  //final_payoff幣別=USD
							break;

							case 2:  //Compo
							final_payoff[k - 1] = final_payoff[k - 1] + exp(-rf*(currentPeriod - n_past)*dt)*ParticipationRate*((FXt*R_final[currentPeriod] - FX0*spot[0]));  //final_payoff幣別=USD
							break;

							case 3:  //Quanto
							final_payoff[k - 1] = final_payoff[k - 1] + exp(-rf*(currentPeriod - n_past)*dt)*ParticipationRate*(R_final[currentPeriod] - spot[0]) * FXT;  //final_payoff幣別=USD
							break;
							*/
						}

						discount_payoff = discount_payoff + final_payoff[k - 1];
						count[2] += 1;
					}
				}
			}  //end of if (ko_flag == 0)  //開始模擬後,若有KO,則停止模擬

			   //觀察期結束前判斷是否已KO
			if (ko_flag == 1) {  //KO則拿利息
				coupon_sum += re_coupon[ko_period];
				count[0] += 1;
			}
			else
			{
				coupon_sum += 0;
			}
		}  //end of for (k = 1; k <= n_trial; k++)
	}  //end of if (ko_flag == 0)  //若尚未KO,則開始模擬程序

	switch (type)
	{
	case 0:  //Normal
		payoff_sum = (discount_payoff + coupon_sum * put_k);  //加總全部
		break;
		/*
		case 1:  //Floating
		payoff_sum = (discount_payoff + coupon_sum * put_k * FXt);  //加總全部;幣別=USD
		break;

		case 2:  //Compo
		payoff_sum = (discount_payoff + coupon_sum * put_k * FX0);  //加總全部;幣別=USD
		break;

		case 3:  //Quanto
		payoff_sum = (discount_payoff + coupon_sum * put_k * FXT);
		break;
		*/
	}

	//payoff_sum = (discount_payoff + coupon_sum * put_k);  //加總全部
	payoff_average = payoff_sum / n_trial;

	return payoff_average;   //payoff_average=選擇權價格$  => payoff_average/put_k=選擇權價格%
}

double _stdcall myTurboNote_Delta(long type, long no_observation, long double y, double r, double ParticipationRate,
	double put_k, double H, double bonuscouponRate,
	double* spot, double* spot2, double* vol, double* stk, double t, double FX0)
{
	long no_asset = 1;
	long ko_flag_greeks = 0;
	long no_period = no_observation + 1;
	long no = 0;

	int i = 0;
	double p1 = 0.0;
	double p2 = 0.0;
	double delta;
	//double ds = 0.0001*spot[no];
	double ds = 0.0001;

	double* spot_u;
	double* spot_d;
	spot_u = new double[no_asset];
	spot_d = new double[no_asset];
	for (i = 0; i < no_asset; i++) {
		spot_u[i] = spot2[i];
		spot_d[i] = spot2[i];
	}

	spot_u[no] = spot_u[no] + ds;
	spot_d[no] = spot_d[no] - ds;

	ko_flag_greeks = isKnockOut(type, no_observation, y, r, ParticipationRate, put_k, H, bonuscouponRate, spot, spot2, vol, stk, t, FX0);

	if (ko_flag_greeks == 0) {
		p1 = myTurboNote_Price(type, no_observation, y, r, ParticipationRate, put_k, H, bonuscouponRate, spot, spot_u, vol, stk, t, FX0);
		p2 = myTurboNote_Price(type, no_observation, y, r, ParticipationRate, put_k, H, bonuscouponRate, spot, spot_d, vol, stk, t, FX0);

		delta = (p1 - p2) / (2 * ds);
	}
	else {
		delta = 0.0;
	}
	return delta;
}

double _stdcall myTurboNote_Gamma(long type, long no_observation, long double y, double r, double ParticipationRate,
	double put_k, double H, double bonuscouponRate,
	double* spot, double* spot2, double* vol, double* stk, double t, double FX0)
{
	long no_asset = 1;
	long ko_flag_greeks = 0;
	long no_period = no_observation + 1;
	long no = 0;

	int i = 0;
	double p1 = 0.0;
	double p2 = 0.0;
	double gamma;
	//double ds = 0.0001*spot[no];
	double ds = 0.0001;

	double* spot_u;
	double* spot_d;
	spot_u = new double[no_asset];
	spot_d = new double[no_asset];
	for (i = 0; i < no_asset; i++) {
		spot_u[i] = spot2[i];
		spot_d[i] = spot2[i];
	}

	spot_u[no] = spot_u[no] + ds;
	spot_d[no] = spot_d[no] - ds;

	ko_flag_greeks = isKnockOut(type, no_observation, y, r, ParticipationRate, put_k, H, bonuscouponRate, spot, spot2, vol, stk, t, FX0);

	if (ko_flag_greeks == 0) {
		p1 = myTurboNote_Delta(type, no_observation, y, r, ParticipationRate, put_k, H, bonuscouponRate, spot, spot_u, vol, stk, t, FX0);
		p2 = myTurboNote_Delta(type, no_observation, y, r, ParticipationRate, put_k, H, bonuscouponRate, spot, spot_d, vol, stk, t, FX0);

		gamma = (p1 - p2) / (2 * ds);
	}
	else {
		gamma = 0.0;
	}
	return gamma;
}
double _stdcall myTurboNote_Vega(long type, long no_observation, long double y, double r, double ParticipationRate,
	double put_k, double H, double bonuscouponRate,
	double* spot, double* spot2, double* vol, double* stk, double t, double FX0)
{
	long no_asset = 1;
	long ko_flag_greeks = 0;
	long no_period = no_observation + 1;
	long no = 0;

	int i = 0;
	double p1 = 0.0;
	double p2 = 0.0;
	double vega;
	double dv = 0.0001;

	double* vol_u;
	double* vol_d;
	vol_u = new double[no_asset];
	vol_d = new double[no_asset];
	for (i = 0; i < no_asset; i++) {
		vol_u[i] = vol[i];
		vol_d[i] = vol[i];
	}

	vol_u[no] = vol_u[no] + dv;
	vol_d[no] = vol_d[no] - dv;

	ko_flag_greeks = isKnockOut(type, no_observation, y, r, ParticipationRate, put_k, H, bonuscouponRate, spot, spot2, vol, stk, t, FX0);

	if (ko_flag_greeks == 0) {
		p1 = myTurboNote_Price(type, no_observation, y, r, ParticipationRate, put_k, H, bonuscouponRate, spot, spot2, vol_u, stk, t, FX0);
		p2 = myTurboNote_Price(type, no_observation, y, r, ParticipationRate, put_k, H, bonuscouponRate, spot, spot2, vol_d, stk, t, FX0);

		vega = (p1 - p2) / (2 * dv);
	}
	else {
		vega = 0.0;
	}
	return vega;
}
long _stdcall myTurboNote_isKnockOut(long type, long no_observation, long double y, double r, double ParticipationRate,
	double put_k, double H, double bonuscouponRate,
	double* spot, double* spot2, double* vol, double* stk, double t, double FX0)
{
	long ans = isKnockOut(type, no_observation, y, r, ParticipationRate, put_k, H, bonuscouponRate, spot, spot2, vol, stk, t, FX0);

	return ans;
}