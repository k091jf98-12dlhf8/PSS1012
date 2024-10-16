#include <set>
#include <map>
#include <unordered_map>
#include <ctime>
#include <time.h>
#include <cmath>
#include <algorithm>
#include "class.h"
#include "class2.h"
#include "LF.h"
#include "PISketch.h"
#include "strawman.h"
#include "para.h"

using namespace std;



uint32_t insert_data[MAX_INSERT_PACKAGE]; 
uint32_t query_data[MAX_INSERT_PACKAGE];

int load_data(){
	FILE *pf = fopen(filename, "rb");
	ground_truth_f.clear();
	ground_truth_p.clear();
	ground_truth_w.clear();
	char ip[13];
    int ret = 0;
	int now_w = 0;
    while (fread(ip, 1, 13, pf)) {  //13 bytes for a packet
        uint32_t key = *(uint32_t *) ip;  //use the first 4 bytes(32bit) as the key (copy from ColdFilter)
        insert_data[ret] = key;   //a new insert
        ground_truth_f[key]++;
		if(ground_truth_w.find(key) == ground_truth_w.end()){  //no recent windows info, new flow
			ground_truth_p[key]++;
			ground_truth_w[key]=now_w;
		}
		else if(ground_truth_w[key]<now_w){ //find recent window info, but this is the first time in current window
			ground_truth_p[key]++;
			ground_truth_w[key]=now_w;
		}
		else{;}  //find recent window and the flow has come in current window, do nothing
		ret++;
		if(ret%WindowSize==0)now_w++;  //a window passed
        if (ret == MAX_INSERT_PACKAGE)
            break;
    }
    fclose(pf);

    int i = 0;
    for (auto itr: ground_truth_f) {
        query_data[i++] = itr.first;
    }

    printf("Total stream size = %d\n", ret);
    printf("Distinct item number = %d\n", ground_truth_f.size());

    return ret;
}

inline double ABS(double a, double b){
	if(a>b)return a-b;
	return b-a;
}






int main(int argc,char **argv){   

	int pkt_num = load_data();
  

	//Pres_Filter() PF;
	//auto FS=Freq_Sketch(L2_recourd_num);   //L2构造
	//auto LF=P_Filter(L1_bucket_num,L1_entry_num,L1_ID_len,L1_F_len,L1_P_len,L1_w_len,P_thr,Seed);  //L1构造
	
	//PSSketch<P_Filter,Freq_Sketch> PSSketch_Model(LF,FS);   //PSSketch构造，并指定L1和L2

	auto SLayer = S_Sketch();
	auto LLayer = L_Sketch();
	LPSSketch LPS(SLayer,LLayer); //LPSSketch构造
	
	BloomFilter BF(3,BF_len);    //Bloom Filter构造
	PISketch PISketch_Model(PI_X,PI_Y); //PISketch构造
	Strawman Strawman_Model; //Starwman构造
	


printf("\n\n-------------------------------------------------------------------------\n");

printf("\n\n****** PSSketch ******\n");

printf("\n\n-------------------------------------------------------------------------\n\n");

//para
	printf("Current Version:V1.0\n");
	printf("Data:%s  Pkt_num:%d  FLow_num:%d\n",filename,pkt_num,ground_truth_f.size());
	printf("P_thr:%d  D_thr:%d  W_thr:%d  L:%d  TopN:%d\n",P_thr,D_thr,W_thr,L,TopN);
	int PI_num=0;

	
	int* P_sta = (int*)malloc(9999999 * sizeof(int));
	int* f_sta = (int*)malloc(9999999 * sizeof(int));
	int* d_sta = (int*)malloc(9999999 * sizeof(int));
	
	if (P_sta == NULL || f_sta == NULL || d_sta == NULL) {
		fprintf(stderr, "Memory allocation failed\n");
		exit(1);
	}

	FILE *P_dist = fopen("p.log","a");
	FILE *F_dist = fopen("f.log","a");
	FILE *D_dist = fopen("den.log","a");
	for (auto itr: ground_truth_f) {
		uint32_t k = itr.first;
		double den = (double)(ground_truth_f[k])/ground_truth_p[k];
		P_sta[ground_truth_p[k]]++;
		f_sta[itr.second]++;
		d_sta[(int)((den-1)*1000)]++;
		if(ground_truth_p[k]>P_thr && den<D_thr){
			//printf("Flow %d (f=%d,p=%d): %f\n",k,ground_truth_f[k],ground_truth_p[k],den);
			PI_num++;
		}
	}
	printf("Totally %d PI-flows in %d (%f\%) \n",PI_num,ground_truth_f.size(),(double)PI_num/ground_truth_f.size()*100);
	
	printf("\n\n-------------------------------------------------------------------------\n\n");

	for(int i=0;i<99999;i++){
			if(f_sta[i])fprintf(F_dist,"%d\t%d\n",i,f_sta[i]);
			if(P_sta[i])fprintf(P_dist,"%d\t%d\n",i,P_sta[i]);
			if(d_sta[i])fprintf(D_dist,"%f\t%d\n",(double)i/1000+1,d_sta[i]);
		}

	fclose(P_dist);
	fclose(F_dist);
	fclose(D_dist);



/*****************************
          LPSSketch
******************************/ 

//Entrance & throuput
	clock_t Tstart_LPS = clock();
	for(int i=0;i<pkt_num;i++){
		LPS.insert(insert_data[i]);
	}
	clock_t Tend_LPS = clock();
	double time_taken_LPS = ((double)(Tend_LPS - Tstart_LPS)) / CLOCKS_PER_SEC;
	printf("LPS-Memory:%f (KB) --- %f+%f\n",((double)(S_bucket_num * S_entry_num * S_size + L_record_num * L_size)/1024),
		  									  (double)(S_bucket_num * S_entry_num * S_size)/1024,
											  (double)L_record_num * L_size/1024
											);
	printf("LPS-Throughput: %f (Mpps)\n",pkt_num/time_taken_LPS/1000000);
	
//F1
	int TP_LPS = 0;
	int FP_LPS = 0;
	int TF_LPS = 0;
	vector<report_item> ReportList_LPS =  LPS.getPI(1);
 	for(report_item Item : ReportList_LPS){
		if(Item.p>P_thr && Item.D<D_thr){
			int true_P = ground_truth_p[Item.ID];
			double true_D = (double)(ground_truth_f[Item.ID])/ground_truth_p[Item.ID];
			if(true_P>P_thr && true_D<D_thr){TP_LPS++;}
			else{FP_LPS++;}
		}
	}
	double Acc_LPS = (double)TP_LPS/(TP_LPS+FP_LPS);
	double Rec_LPS = (double)TP_LPS/PI_num;
	printf("LPS-Acc:%f\n",Acc_LPS);	
	printf("LPS-Rec:%f\n",Rec_LPS);  //REC
	printf("LPS-F1:%f\n",2*Acc_LPS*Rec_LPS/(Acc_LPS+Rec_LPS));
//ARE
	int fk_LPS=0;
	double errsum_LPS = 0.0;
	vector<report_item> LList_LPS =  LPS.getPI(0);
	for(report_item Item : LList_LPS){
		double report_D = Item.D;
		double true_D = (double)(ground_truth_f[Item.ID])/ground_truth_p[Item.ID];
		errsum_LPS+=ABS(report_D,true_D);
		fk_LPS++;
	}
	printf("LPS-Errsum:%f fk:%d\n",errsum_LPS,fk_LPS);
	printf("LPS-ARE:%f\n",errsum_LPS/fk_LPS);

	printf("\n\n-------------------------------------------------------------------------\n\n");
	

/*****************************
           PSSketch
******************************/ 

/*

//Entrance & throuput
	clock_t Tstart_PSS = clock();
	int i;
	for(i=0;i<pkt_num;i++){
		PSSketch_Model.insert(insert_data[i]);
	}
	clock_t Tend_PSS = clock();
	double time_taken_PS = ((double)(Tend_PSS - Tstart_PSS)) / CLOCKS_PER_SEC;
	//printf("PS-Memory:%f (KB)\n",(L1_bucket_num*L1_entry_num*L1_size+L2_recourd_num*L2_size)/1024);
	//printf("PS-Throughput: %f (Mpps)\n",pkt_num/time_taken_PS/1000000);
// F1
	int TP=0;
	int FP=0;
	int TF=0;
	vector<report_item> ReportList =  PSSketch_Model.getL2PI_report();
	for(report_item Item : ReportList){
		int true_p = ground_truth_p[Item.ID];
		double true_D = (double)(ground_truth_f[Item.ID])/ground_truth_p[Item.ID];
		if(true_p>P_thr && true_D<D_thr){TP++;}
		else{
			//printf("判断失误！流%d的实际答案为(%d,%f),而上报的结果位(%d,%f)\n",Item.ID,true_p,true_D,Item.p,Item.D);
			FP++;
		}
	}
	double Acc = (double)TP/(TP+FP);
	double Rec = (double)TP/PI_num;
	//printf("PS-Acc:%f\n",Acc);	
	//printf("PS-Rec:%f\n",Rec);  //REC
	//printf("PS-F1:%f\n",2*Acc*Rec/(Acc+Rec));
// ARE
	int fk;
	double errsum = 0.0;
	int L2_item_num = PSSketch_Model.query(3);
	for(fk=0;fk<L2_item_num;fk++){
		uint32_t findkey = PSSketch_Model.getL2fp(fk);
		double ans = (double)(ground_truth_f[findkey])/ground_truth_p[findkey];
		errsum+=ABS(PSSketch_Model.getL2den(fk),ans);

		
		printf("Flow %d -- f=%d p=%d d=%f\n",findkey,
								ground_truth_f[findkey],
								ground_truth_p[findkey],
								(double)(ground_truth_f[findkey])/ground_truth_p[findkey]);
		
		for(int x=0;x<L1_bucket_num;x++){
			for(int y=0;y<L1_entry_num;y++){
					if(LF.Filter[x].ID[y]==findkey){
						printf("INFO in L1: f=%d p=%d d=%f\n",LF.Filter[x].f[y],LF.Filter[x].p[y],LF.Filter[x].den[y]);
					}
				}
		}
		
	}
	//printf("PS-Errsum:%f fk:%d\n",errsum,fk);
	//printf("PS-ARE:%f\n",errsum/fk);
	
	//printf("\n\n-------------------------------------------------------------------------\n\n");

*/



/*****************************
           PISketch
******************************/ 

	int i;
//Entrance & throuput
	clock_t Tstart_PIS = clock();
	int w_cnt = 0;
	for(i=0;i<pkt_num;i++){
		int w;
		//printf("%d\n",insert_data[i]);
        if(BF.query(insert_data[i]))
            w = -1;
        else {
            BF.insert(insert_data[i]);
            w = L;  //init L
        }
		PISketch_Model.insert(insert_data[i],w);
		w_cnt++;
		if(w_cnt%WindowSize==0){
			BF.clear();
		}
	}
	clock_t Tend_PIS = clock();
	double time_taken_PI = ((double)(Tend_PIS - Tstart_PIS)) / CLOCKS_PER_SEC;
	printf("PI-Memory:%f (KB)\n",(BF_size*BF_len+PI_X*PI_Y*PI_size)/1024);
	printf("PI-Throughput: %f (Mpps)\n",pkt_num/time_taken_PI/1000000);

	int TP_PIS=0;
	int FP_PIS=0;
// F1
 	for(int pix=0;pix<PI_X;pix++){
 		for(int piy=0;piy<PI_Y;piy++){
			int tmp_p = PISketch_Model.bucket[pix].p[piy];    
			int tmp_f = PISketch_Model.bucket[pix].p[piy]*(L+1)-PISketch_Model.bucket[pix].w[piy]; 
			double tmp_d = (double)tmp_f/tmp_p;
			if(tmp_p>P_thr && tmp_d<D_thr){  //calc the p and d with PISketch to find PI flow but not W
			//if(PISketch_Model.bucket[pix].w[piy]>W_thr){
				uint32_t findkey = PISketch_Model.getfp(pix,piy);
				if(ground_truth_p[findkey]>P_thr && (double)(ground_truth_f[findkey])/ground_truth_p[findkey]<D_thr){TP_PIS++;}  //报告正确
				else{
					//printf("判断失误！流%d的实际答案为(%d,%f)\n",reportID,true_p,true_D);
					FP_PIS++;}  //报告错误
				}
 			
 		}
 	}
	double Acc_PI = (double)TP_PIS/(TP_PIS+FP_PIS);
	double Rec_PI = (double)TP_PIS/PI_num;
	printf("PI-Acc:%f\n",Acc_PI); 
	printf("PI-Rec:%f\n",Rec_PI);
	printf("PI-F1:%f\n",2*Acc_PI*Rec_PI/(Acc_PI+Rec_PI));
// ARE  
	int fk_PI=0;
	double errsum_PI = 0.0;
	for(int pix=0;pix<PI_X;pix++){
		for(int piy=0;piy<PI_Y;piy++){
			fk_PI++;
			// W = p*L-(f-p) = pL-f+p = (L+1)p-f
			// f = (L-1)p-w 
			int PIS_f = (L+1)*PISketch_Model.bucket[pix].p[piy]- PISketch_Model.bucket[pix].w[piy];
			double PIS_den = (double)PIS_f/PISketch_Model.bucket[pix].p[piy];
			uint32_t findkey = PISketch_Model.getfp(pix,piy);
			double ans_den = (double)(ground_truth_f[findkey])/ground_truth_p[findkey];
			errsum_PI+=ABS(PIS_den,ans_den);
		}
	}
	
	printf("PI-Errsum:%f fk:%d\n",errsum_PI,fk_PI);
	printf("PI-ARE:%f\n",errsum_PI/fk_PI);

	printf("\n\n-------------------------------------------------------------------------\n\n");

	



/*****************************
           Strawman
******************************/ 

//Entrance & throuput
	clock_t Tstart_SM = clock();
	for(i=0;i<pkt_num;i++){
		Strawman_Model.insert(insert_data[i]);
	}
	clock_t Tend_SM = clock();
	double time_taken_SM = ((double)(Tend_SM - Tstart_SM)) / CLOCKS_PER_SEC;
	
	//3 shows the num of hash functions,OO and CM have the same len, PI_list should be calculaated.
	printf("SM-Memory:%f (KB)\n",(3*OO_PE_Len*OO_size+3*OO_PE_Len*CM_size+Strawman_Model.PI_List.size()*12)/1024);  
	printf("SM-Throughput: %f (Mpps)\n",pkt_num/time_taken_SM/1000000);

// F1
	int TP_SM=0;
	int FP_SM=0;
	for(const auto& reportItem : Strawman_Model.PI_List){
		int true_p = ground_truth_p[reportItem.first];
		double true_D = (double)(ground_truth_f[reportItem.first])/ground_truth_p[reportItem.first];
		if(true_p>P_thr && true_D<D_thr){TP_SM++;}
		else{FP_SM++;}
	}
	double Acc_SM = (double)TP_SM/(TP_SM+FP_SM);
	double Rec_SM = (double)TP_SM/PI_num;
	printf("SM-Acc:%f\n",Acc_SM);	
	printf("SM-Rec:%f\n",Rec_SM);  //REC
	printf("SM-F1:%f\n",2*Acc_SM*Rec_SM/(Acc_SM+Rec_SM));	

// ARE
	int fk_SM=0;
	double errsum_SM = 0.0;
	for(const auto& reportItem : Strawman_Model.PI_List){
		fk_SM++;
		double SM_den = reportItem.second;
		double ans_den =  (double)(ground_truth_f[reportItem.first])/ground_truth_p[reportItem.first];
		errsum_SM += ABS(SM_den,ans_den);
	}

	printf("SM-Errsum:%f fk:%d\n",errsum_SM,fk_SM);
	printf("SM-ARE:%f\n",errsum_SM/fk_SM);

	printf("\n\n-------------------------------------------------------------------------\n\n");







/************************************
                Top N
*************************************/ 
	
/*
	printf("Save TopN to log file...\n");

	FILE *TopN_PIS = fopen("PIS.log","a");
	FILE *TopN_PSS = fopen("PSS.log","a");
	FILE *TopN_STM = fopen("STM.log","a");
//PISketch TopN
	vector<pair<uint32_t,int>> PISketch_Top;
	for(int pix=0;pix<10;pix++){
 		for(int piy=0;piy<1000;piy++){
 			uint32_t findkey = PISketch_Model.getfp(pix,piy);
			int w = PISketch_Model.getw(pix, piy);
			if(w>=W_thr){
				PISketch_Top.push_back({findkey,w});
			}
 		}
 	}
	std::sort(PISketch_Top.begin(), PISketch_Top.end(), 
              [](const std::pair<uint32_t, int>& a, const std::pair<uint32_t, int>& b) {
                  return a.second > b.second; 
              });
	fprintf(TopN_PIS,"PISketch Top%d\n",TopN);
	fprintf(TopN_PIS,"ID \t w \t f \t p \t d \n");
	for(int cnt=0;cnt<TopN;cnt++){
		fprintf(TopN_PIS,"%d \t %d \t %d \t %d \t %f \n",PISketch_Top[cnt].first,PISketch_Top[cnt].second,
											  ground_truth_f[PISketch_Top[cnt].first],
											  ground_truth_p[PISketch_Top[cnt].first],
											  (double)ground_truth_f[PISketch_Top[cnt].first]/ground_truth_p[PISketch_Top[cnt].first]
											  );
	}

//PSSketch TopN
	vector<pair<uint32_t,double>> PSSketch_Top;
	for(int x=0;x<L2_recourd_num;x++){
		if(FS.PS_record[x].fp!=0 && FS.PS_record[x].den<D_thr){
			PSSketch_Top.push_back({FS.PS_record[x].fp,FS.PS_record[x].den});
			//printf("%d -- %f\n",FS.PS_record[x].fp,FS.PS_record[x].den);
		}
	}
	std::sort(PSSketch_Top.begin(), PSSketch_Top.end(), 
              [](const std::pair<uint32_t, double>& a, const std::pair<uint32_t, double>& b) {
                  return a.second < b.second; 
              });
	fprintf(TopN_PSS,"PSSketch Top%d\n",TopN);
	fprintf(TopN_PSS,"ID \t f \t p \t d \n");
	for(int cnt=0;cnt<TopN;cnt++){
		fprintf(TopN_PSS,"%d \t %d \t %d \t %f \n",PSSketch_Top[cnt].first,
										ground_truth_f[PSSketch_Top[cnt].first],
										ground_truth_p[PSSketch_Top[cnt].first],
										(double)ground_truth_f[PSSketch_Top[cnt].first]/ground_truth_p[PSSketch_Top[cnt].first]
										);
		//printf("%d \t %d \t %d \t %f \n",PSSketch_Top[cnt].first,
		//								PSSketch_Model.getL1f(PSSketch_Top[cnt].first),
		//								PSSketch_Model.getL1p(PSSketch_Top[cnt].first),
		//								PSSketch_Top[cnt].second
		//								);
	}
	//printf("Cur_ptr: %d\n",PSSketch_Model.L2.cur_ptr);

//Strawman TopN
	vector<pair<uint32_t,double>> Strawman_Top = Strawman_Model.PI_List;
	std::sort(Strawman_Top.begin(), Strawman_Top.end(), 
              [](const std::pair<uint32_t, double>& a, const std::pair<uint32_t, double>& b) {
                  return a.second < b.second; 
              });
	fprintf(TopN_STM,"Strawman Top%d\n",TopN);
	fprintf(TopN_STM,"ID \t f \t p \t d \n");
	for(int cnt=0;cnt<TopN;cnt++){
		fprintf(TopN_STM,"%d \t %d \t %d \t %f \n",Strawman_Top[cnt].first,
										ground_truth_f[Strawman_Top[cnt].first],
										ground_truth_p[Strawman_Top[cnt].first],
										(double)ground_truth_f[Strawman_Top[cnt].first]/ground_truth_p[Strawman_Top[cnt].first]
										);
	}
	

	fclose(TopN_PIS);
	fclose(TopN_PSS);
	fclose(TopN_STM);

	PIN('o');

*/

	printf("Done.\n");
	
}






