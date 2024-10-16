
#ifndef CLASS2_H
#define CLASS2_H

#include "hash.h"
#include "BOBHash32.h"
#include <cstring>
#include <algorithm>
#include <vector>
#include <stdint.h>
#include "para.h"
#include <time.h>
#include <stdlib.h>
#include <immintrin.h>

using namespace std;

class Bucket_S
{
public:
    S_ID_type *FP;  //fp of a flow
    S_F_type *f;  //freq 
	S_P_type *p;  //pers(6 bit) and flag(2bit)

    Bucket_S() {
        FP = new S_ID_type[S_entry_num];
        f = new S_F_type[S_entry_num];
		p = new S_P_type[S_entry_num];
        memset(FP, 0, sizeof(S_ID_type) * S_entry_num);
		memset(f, 0, sizeof(S_F_type) * S_entry_num);
		memset(p, 0, sizeof(S_P_type) * S_entry_num);
    }
};



class Bucket_L
{
public:
    L_ID_type *ID;  //full ID of a flow
    L_of_type *of;  //freq overflow & pers overflow 
    Bucket_L() {
        ID = new L_ID_type[L_record_num];
        of = new L_of_type[L_record_num];
        memset(ID, 0, sizeof(L_ID_type) * L_record_num);
		memset(of, 0, sizeof(L_of_type) * L_record_num);
    }
};




class L_virtual
{
public:
	virtual int report(uint32_t ID, uint8_t f,uint8_t p,char mode) = 0; 
};




class S_Sketch{
public:
    Bucket_S *Filter;  
    BOBHash32 *P_Hash;
	int w_cnt;

	L_virtual* L_vir;

    S_Sketch() {
        Filter = new Bucket_S[S_bucket_num];
        for (int i = 0; i < S_bucket_num; i++){
            Filter[i] = Bucket_S();
        }
        
        P_Hash = new BOBHash32(97);
		w_cnt = 0;
    }

/*
P: 000000 0 0
0xFD:111111 0 1
0xFE:111111 1 0
0x03:000000 1 1
*/

bool insert(uint32_t key, int count = 1){
    	w_cnt++;
		if (w_cnt == WindowSize) {
			w_cnt = 0;
			__m256i mask = _mm256_set1_epi8(0xFD);	
			for (int i = 0; i < S_bucket_num; i++) {
				__m256i data = _mm256_loadu_si256((__m256i*)Filter[i].p);  
				data = _mm256_and_si256(data, mask);  
				_mm256_storeu_si256((__m256i*)Filter[i].p, data); 
			}
		}
	
		
		S_ID_type keyfp = key%pow(2,S_ID_len);
		int hash_index = P_Hash->run((const char *)&key, 4) % S_bucket_num;
		Bucket_S Bid = Filter[hash_index];

		int match_index = -1; //自己的位置
		int first_empty=-1;  //第一个空位
		int p_min_pos=-1;    //最小的p所在位置		
		int p_min_value=99;    //最小的p所在位置		
		
	 
		__m256i target_FP = _mm256_set1_epi32(key); 
		__m256i target_empty = _mm256_set1_epi32(0);
   		for (int i = 0; i < 32; i += 8) {
			__m256i fp_data = _mm256_loadu_si256((__m256i*)&Bid.FP[i]);
			//find self
			__m256i fp_cmp = _mm256_cmpeq_epi32(fp_data, target_FP);
        	int fp_mask = _mm256_movemask_epi8(fp_cmp);
			if (fp_mask != 0) {
            	match_index = (_tzcnt_u32(fp_mask) / 4)+i;  //find!
            	goto search_over;
        	}
			//find empty
			__m256i empty_cmp = _mm256_cmpeq_epi32(fp_data, target_empty);
       		int empty_mask = _mm256_movemask_epi8(empty_cmp);
        	if (empty_mask != 0) {
	 	        first_empty = (_tzcnt_u32(empty_mask) / 4)+i;
				goto search_over; 
		    }
    	}

/*
			//find min_p
			__m256i min_p_tmp = _mm256_set1_epi8(0xFF); //maximum
			__m256i p_data = _mm256_loadu_si256((__m256i*)&Bid.p[0]);
			__m256i p_high6 = _mm256_srli_epi16(p_data, 2);  
		    p_high6 = _mm256_and_si256(p_high6, _mm256_set1_epi8(0x3F)); 
			__m256i p_min_cmp = _mm256_min_epu8(p_high6, min_p_tmp);
			int min_mask = _mm256_movemask_epi8(_mm256_cmpeq_epi8(p_min_cmp, p_high6));
		    if (min_mask != 0) {
            	min_p_tmp = p_min_cmp;
		        p_min_pos = (_tzcnt_u32(min_mask));
		    }
*/			
		//走到这里说明没有自己也没有空位，需要寻找最小的p
		//预备方案
		for(int index=0;index<S_entry_num;index++){
			if((Bid.p[index]>>2)<p_min_value && (Bid.p[index]&0x1)==0){
				p_min_value = (Bid.p[index]>>2);
				p_min_pos = index;
			}
		}
			

	search_over:
	
		if(match_index>=0){  //已存在
			int i = match_index;
			Bid.f[i]++;
			if(Bid.f[i]==pow(2,S_F_len)-1){  //blow the f!
				Bid.f[i]=0;  //reset f for per 'overflow'
				if((Bid.p[i]&0x1)==0){ //not a pers flow, a trash flow, kick!	
					Bid.FP[i]=0;
					Bid.f[i]=0;
					Bid.p[i]=0;
					return -1; //kicked
				}
				int res = L_vir->report(key,0,0,'f'); //only need f_of +1
				if(res<0){
					Bid.FP[i]=0;
					Bid.f[i]=0;
					Bid.p[i]=0;
					return -1;  //kicked by L_layer
				}				
			}
			if(((Bid.p[i]>>1)&0x1)==0){  //窗口 
				uint8_t p_tmp = (Bid.p[i]>>2)+1;
				Bid.p[i] = (Bid.p[i] & 0x03)|((p_tmp & 0x3F)<<2); //p+1
				Bid.p[i] |= 0x2; 
				if((Bid.p[i]>>2)>=P_thr){ //blow the p, report it!
					if((Bid.p[i]&0x1)==0){ //first time report
						Bid.p[i] |= 0x1; //r=1
						int res = L_vir->report(key,Bid.f[i],Bid.p[i],'o');   
						if(res<0){
							Bid.FP[i]=0;
							Bid.f[i]=0;
							Bid.p[i]=0;
							return -1;  //kicked by L_layer
						}
					}
					else{
						int res = L_vir->report(key,0,0,'p');  
						if(res<0){
							Bid.FP[i]=0;
							Bid.f[i]=0;
							Bid.p[i]=0;
							return -1;  //kicked by L_layer
						}
					}
					Bid.p[i] &= 0x3;  //reset p for per 'overflow'
				}
			}
			return 1;  //找到了自己，完成记录后结束
		}		
		
	    // No record found
		if(first_empty>=0){		//有空位
			Bid.FP[first_empty]=keyfp;
			Bid.f[first_empty]=1;
			Bid.p[first_empty]=0;
			Bid.p[first_empty]|=0x06;  //p=1,w=1,r=0   000001 1 0
	        return 1;  //success
		}
		else{  //没有空位
			bool kick_op = p_min_value<P_thr;  //如果>=，则一定不能踢
			static uint32_t rng_state = static_cast<uint32_t>(time(0));
			rng_state ^= rng_state << 13;
        	rng_state ^= rng_state >> 17;
        	rng_state ^= rng_state << 5;
       		uint32_t random_number = rng_state;
			double random_value = random_number / static_cast<double>(UINT32_MAX);
		
			if(LPS_POS){  //概率替换
				if (random_value > 1.0/p_min_value) {kick_op=0;}  
			}
			if(kick_op){   //踢！
				Bid.FP[p_min_pos] = keyfp; //replace it!
				Bid.f[p_min_pos]=1;
				Bid.p[p_min_pos]=0;
				Bid.p[p_min_pos]|=0x06;  //p=1,w=1,r=0   000001 1 0
				return 1;
			}
			return 0; //kick fail due to random  
		}

		//no one can be kick, just ingnore the new one
		return 0; //fail			
   }

	void nextlevel(L_virtual * nl){
        L_vir = nl;
    }
	
};





class L_Sketch: public L_virtual{ 
public:
    Bucket_L *PIRecord;
	int cur_ptr;
	int cur_max;
	uint32_t tmp_w;

    L_Sketch(){
        PIRecord = new Bucket_L();		
		cur_ptr=0;
		tmp_w = 0;
    }

#define Burst_Multi 1.1  //how to define burst

int hash_col = 0;
int no_col =0;

    int report(uint32_t ID, uint8_t f,uint8_t p,char mode) {

		int op = -1;
		no_col ++;
	
		if(mode=='o'){ //first time report
			if(LPS_debug){
				FILE *log = fopen("LPS.log","a");	
				fprintf(log,"流%u第一次上报到大层，被放在%d位置\n",ID,cur_ptr);
				fclose(log);
			}
			PIRecord->ID[cur_ptr]=ID;
			PIRecord->of[cur_ptr]=0x1;   //0..0 0..1
			cur_max++;
			while(PIRecord->ID[cur_ptr]!=0)cur_ptr++;
			//printf("    L_layer:put at %d\n",cur_ptr-1);
			//if(cur_ptr%10000==0){printf("%d\n",cur_ptr);}
			assert(cur_ptr<L_record_num-1);  //L layer full, kick?    //TODO
			return 1;  //success
		}
		else if(mode=='f'){ //f overflow
			for(int i=0;i<cur_max;i++){
				if(PIRecord->ID[i]==ID){    //Avoid Overflow!!! Unreasonable!!!
					if((PIRecord->of[i]>>L_f_of_len)>=pow(2,L_f_of_len)){
						hash_col++;
						if(LPS_debug){
							FILE *log = fopen("LPS.log","a");	
							fprintf(log,"【错误%d/%d】 L_layer:ID:%u (fp:%u) 在更新f时溢出或p已满\n",hash_col,no_col,ID,ID%pow(2,S_ID_len));
							fclose(log);
						}
						goto post;
					}
					uint8_t f_tmp = (PIRecord->of[i] >> L_p_of_len)+1;
					PIRecord->of[i] = (PIRecord->of[i] & 0xFF | ((f_tmp&0xFF)<<L_p_of_len)); //f+1
					op = i;
					goto post;
				}
				//assert((PIRecord->of[i]>>4)<15);  //check overflow!
			}
			hash_col++;
			if(LPS_debug){
				FILE *log = fopen("LPS.log","a");	
				fprintf(log,"【错误%d/%d】 L_layer:ID:%u (fp:%u) 在更新f时没有找到项\n",hash_col,no_col,ID,ID%pow(2,S_ID_len));
				fclose(log);
			}
			//assert(0); //impossible now, unless kicked     //TODO
		}
		else if(mode=='p'){  //p overflow
			for(int i=0;i<cur_max;i++){
				if(PIRecord->ID[i]==ID){
					if((PIRecord->of[i]&0xFF)>=pow(2,L_p_of_len)){
						hash_col++;
						if(LPS_debug){
							FILE *log = fopen("LPS.log","a");	
							fprintf(log,"【错误%d/%d】 L_layer:ID:%u (fp:%u) 在更新p时溢出\n",hash_col,no_col,ID,ID%pow(2,S_ID_len));
							fclose(log);
						}
						goto post;
					}
					PIRecord->of[i] = (PIRecord->of[i] & 0xFF00) | ((PIRecord->of[i] + 1) & 0x00FF); //p+1
					op = i;
					goto post;
				}
				//assert((PIRecord->of[i] & 0xF)<15);  //check overflow!
			}
			hash_col++;
			if(LPS_debug){
				FILE *log = fopen("LPS.log","a");
				fprintf(log,"【错误%d/%d】 L_layer:ID:%u (fp:%u) 在更新p时没有找到项\n",hash_col,no_col,ID,ID%pow(2,S_ID_len));
				fclose(log);
			}
			//assert(0); //impossible now, unless kicked     //TODO
		}
		else{
			assert(0); //impossible now, for future expand;
		}

	post:
		//update complete, do some post
		//assert(op>=0);
		if(op<0){return 0;}  //fail (hash collision, item not found)
		int cur_f = PIRecord->of[op]>>L_f_of_len;
		int cur_p = PIRecord->of[op]&0xFF;
		if(cur_f>=cur_p && cur_p<pow(2,L_p_of_len)){
			PIRecord->ID[op]=0;
			PIRecord->of[op]=0;
			if(cur_ptr>op) cur_ptr = op;   //cur_ptr shows the new smallest empty pos, and will find the next
			if(LPS_debug){
				FILE *log = fopen("LPS.log","a");
				fprintf(log,"流%u被大层踢掉了\n",ID);
				fclose(log);
			}
			return -1;  //kicked, inform S_layer       
		}		

		return 1;  //success
		
    }

};




/***********************************
	Definition for LPSSketch (Main2)
************************************/
class LPSSketch{
public:
	S_Sketch S_layer;
	L_Sketch L_layer;
	BOBHash32 *P_Hash;

	LPSSketch(const S_Sketch& SL, const L_Sketch& LL) : S_layer(SL), L_layer(LL){
		S_layer.nextlevel(&L_layer);
		P_Hash = new BOBHash32(97);
	}

	bool insert(uint32_t key){	
		S_layer.insert(key);
		return 0;
	}		

	vector<report_item> getPI(int mode){  //  mode=0 getAll   mode=1 getPI  
		vector<report_item> res;
		report_item new_PI;
		for(int i=0;i<L_record_num;i++){
			L_ID_type ID = L_layer.PIRecord->ID[i];
			if(ID){
				int f = ((L_layer.PIRecord->of[i])>>L_p_of_len)*(pow(2,S_F_len)-1);
				int p = ((L_layer.PIRecord->of[i])&0xFF)*P_thr;
				if(LPS_debug)printf("%d\t%u\t",i,ID);
				S_ID_type FP = ID%pow(2,L1_ID_len);
				int hash_index = P_Hash->run((const char *)&ID, 4) % S_bucket_num;
					
				Bucket_S Bid = S_layer.Filter[hash_index];
				//printf("在桶%d   ",hash_index);
				for(int j=0;j<S_entry_num;j++){
					if(Bid.FP[j]==ID){
						//printf("的位置%d找到了该流   ",j);
						f+= Bid.f[j];
						p+= (Bid.p[j]>>2); 
					}
				}
				if(LPS_debug)printf("%d\t%d\t",f,p);
				double den = (double)f/p;
				if(LPS_debug)printf("%d\t%d\t%f\t",ground_truth_f[ID],ground_truth_p[ID],(double)(ground_truth_f[ID])/ground_truth_p[ID]);
				if(mode){
					if(den<D_thr){  //Find a PI flow!
						if(LPS_debug)printf("1\n");
						new_PI = {ID,p,den};
						res.push_back(new_PI);
					}
					else{
						if(LPS_debug)printf("0\n");
					}
				}
				else{
					new_PI = {ID,p,den};
					res.push_back(new_PI);
				}
			}
		}
		return res;
	}
	
	bool clear();
	bool init();

};


#endif
