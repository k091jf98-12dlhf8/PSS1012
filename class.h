  #ifndef CLASS_H
#define CLASS_H

#include "hash.h"
#include "BOBHash32.h"
#include <cstring>
#include <algorithm>
#include <vector>
#include <stdint.h>
#include "para.h"
#include <time.h>
#include <stdlib.h>

using namespace std;

int kick_p = 0;
int kick_t = 0;
int kick_max=-1;

/***********************************
	Definition about hash
************************************/

// getFP: cut the hash (int32) got from BOB to suitable length
uint32_t getFP(uint32_t key, int key_len)
{
    static BOBHash32 fpHash(100);
    return fpHash.run((const char *)&key, 4) % 0xFFFF + 1;
}


/***********************************
	Definition for L1 Filter
************************************/

uint32_t cur_w=0;  // global var. show the cur window num, which is added by time


struct PSS_RES{
	int flag;   // 0: new insert   1: old update   2: kicked
	L1_F_type f;
	L1_P_type p;
	int index;
	int pos;
};

// EX bucket definition for PSS filter
class Bucket_PSS
{
public:
    L1_ID_type *ID;  //ID
    L1_F_type *f;  //freq
	L1_P_type *p;  //persi
	//uint16_t *w;  //latest window
	bool *wflag;
	double *den;  //den
	uint32_t *kick_time;  //times for kick
	bool *if_L2;   // 0:the first time to L2   1: not

    Bucket_PSS() {}

    Bucket_PSS(int entry_num)     //we have entry_num items in one hash index (bucket)
    {
        ID = new L1_ID_type[entry_num];
        f = new L1_F_type[entry_num];
		p = new L1_P_type[entry_num];
		//w = new uint16_t[entry_num];
		wflag = new bool[entry_num];
		den = new double[entry_num];
		kick_time = new uint32_t[entry_num];
		
		if_L2 = new bool[entry_num]();
        memset(ID, 0, sizeof(L1_ID_type) * entry_num);
		memset(f, 0, sizeof(L1_F_type) * entry_num);
		memset(p, 0, sizeof(L1_P_type) * entry_num);
		//memset(w, 0, sizeof(uint16_t) * entry_num);
		memset(wflag, 0, sizeof(bool) * entry_num);
		memset(den, 0, sizeof(double) * entry_num);
		memset(kick_time, 0, sizeof(uint32_t) * entry_num);
    }
/*
    void permutation(int p) // permute the p-th item to the first
    {
        for (int i = p; i > 0; --i)
        {
            swap(fp[i], fp[i - 1]);
            swap(counter[i], counter[i - 1]);
        }
    }
*/
};


int lock_b = 0;
int lock[L1_bucket_num]={0};

//class L2_virtual_PSS: A virtual for calling the insert() of L2 Sketch
class L2_virtual_PSS
{
public:
	virtual bool insert(uint32_t ID, int f,int p, bool if_L2) = 0; 
};

// Ex version of ladder filter for our task
class P_Filter{
public:
    Bucket_PSS *Filter;  //we only need one filter for now
    int *kick_time;  
	
    BOBHash32 *P_Hash;
    int bucket_num;
    int entry_num;
	int ID_len,f_len, p_len,w_len;  //digits for the four things
    int P_threshold;   //is_pers
	int w_cnt;

	L2_virtual_PSS* L2_vir;

    P_Filter() {}
    P_Filter(int _bucket_num, int _entry_num,
                 int _ID_len, int _f_len, int _p_len, int _w_len,
                 int _P_threshold,
                 int rand_seed)
        : bucket_num(_bucket_num), entry_num(_entry_num),
          ID_len(_ID_len), f_len(_f_len),p_len(_p_len),w_len(_w_len),
          P_threshold(_P_threshold)
    {
        Filter = new Bucket_PSS[bucket_num];
		kick_time = new int[bucket_num];
        for (int i = 0; i < bucket_num; ++i){
            Filter[i] = Bucket_PSS(entry_num);
			kick_time[i]=0;
        }
        
        P_Hash = new BOBHash32(rand_seed);
		w_cnt = 0;
    }

	inline bool report(int ID,int f,int p, bool if_L2){ 
		if(if_L2){
			//printf("Report %d as an old PS flow: f=%d p=%d\n",ID,f,p);
			return L2_vir->insert(ID,f,p,if_L2);  //old item back
		}  
		//printf("Report %d as a new PS flow: f=%d p=%d\n",ID,f,p);
		return L2_vir->insert(ID,f,p,if_L2);		
	}

	void check_full(int h,int cnt){
		Bucket_PSS Bid = Filter[h];
		bool flag=1;
		for(int i=0;i<entry_num;i++){
			if(Bid.p[i]<P_thr){flag=0;}  //还有可以踢的项目
		}
		if(flag && lock[h]==0){
			lock_b++;
			//printf("桶%d被锁定，当前报文数%d,当前被锁定数:%d\n",h,cnt,lock_b);
			lock[h]=1;
		}
	}


	PSS_RES insert(uint32_t key, int count = 1) 
    {
    	w_cnt++;
		if(w_cnt%WindowSize==0){
			cur_w++;
			for(int i=0;i<bucket_num;i++){
				for(int j=0;j<entry_num;j++){
					Filter[i].wflag[j]=0;
				}	
			}
		}
		//uint32_t keyfp=key;
		L1_ID_type keyfp = key%pow(2,L1_ID_len);
		
        //auto keyfp = getFP(key, ID_len);   //change the uint32 to a counter_len-int
		PSS_RES re = {-1,0,0,0,0};

        // if the item is in B1
        int hash_index[5] = {
        	P_Hash->run((const char *)&key, 4) % bucket_num,
			(P_Hash->run((const char *)&key, 4)+1) % bucket_num,
			(P_Hash->run((const char *)&key, 4)+2) % bucket_num,
			(P_Hash->run((const char *)&key, 4)+3) % bucket_num,
			(P_Hash->run((const char *)&key, 4)+4) % bucket_num
		};
		Bucket_PSS Bid[5] = {
			Filter[hash_index[0]],
			Filter[hash_index[1]],
			Filter[hash_index[2]],
			Filter[hash_index[3]],
			Filter[hash_index[4]]
		};

		//多hash策略：访问5个桶，找到自己就更新，有空位就直接放，如果没有，踢掉5个桶里p最小的

		for(int h=0;h<5;h++){	

		//首先找5个桶中有没有自己
	        for (int i = 0; i < entry_num; ++i)  //in this bucket, we have many pos for item
	            if (Bid[h].ID[i] == keyfp)  //find it
	            {
	            	Bid[h].f[i]++;
					//if(Bid.f[i]>f_max){f_max=Bid.f[i];printf("fmax=%d\n",f_max);}
					assert(Bid[h].f[i]<pow(2,L1_F_len)-1);  //F overflow
					//if(Bid.w[i]!=cur_w){
					if(Bid[h].wflag[i]==0){
						Bid[h].p[i]+=1;
						//if(Bid.p[i]>p_max){p_max=Bid.p[i];printf("pmax=%d\n",p_max);}
						assert(Bid[h].p[i]<pow(2,L1_P_len)-1);  //F overflow
						Bid[h].wflag[i]=1;
						//Bid.p[i]+=(Bid.w[i]^cur_w)&0x1;
						//Bid.w[i]=cur_w;
						
						if(Bid[h].p[i]>P_threshold){  //a pers fow!
							bool report_res = report(key,Bid[h].f[i],Bid[h].p[i],Bid[h].if_L2[i]);  //report to L2
							Bid[h].if_L2[i] = report_res?1:Bid[h].if_L2[i];  // res=1: success report.   res=0: fail, do not change
						}
					}
	                //B1.permutation(i);  //make it to the first	

					re = {1,Bid[h].f[i],Bid[h].p[i],keyfp,i};
					Bid[h].den[i]=(double)(Bid[h].f[i])/Bid[h].p[i];
					//printf("item (%d,%d) update its info f:%d , d:%d\n",hash_index,i,Bid.f[i],Bid.p[i]);

					check_full(hash_index[h],w_cnt);
					return re;  //在5个桶中某一个找到了自己，插入结束
	            }
		}

		//5个桶中都没有，找空位
		for(int h=0;h<5;h++){		
	        // Not found, not full 
	        for (int i = 0; i < entry_num; ++i)
	            if (Bid[h].f[i] == 0)   //an empty pos
	            {
	            	Bid[h].ID[i]=keyfp;
					Bid[h].f[i]=1;
					Bid[h].p[i]=1;
					//Bid.w[i]=cur_w;
					Bid[h].wflag[i]=1;
					
					re = {0,1,1,keyfp,i};
					Bid[h].den[i]=(double)(Bid[h].f[i])/Bid[h].p[i];
					//printf("an item (%d,%d) is added\n",hash_index,i);

	                return re;  //没有找到自己，在5个桶中某一个找到了空位，插入结束
	            }
		}


		//5个桶都没有空位，找最小的p
		uint32_t p_min = (1u<<31)-1;
		int kick[2] = {-1,-1};
		for(int h=0;h<5;h++){		
	        // Not found and full
			for(int i=0;i<entry_num;++i){  //find the min_p (need low than P_thres)
				if(Bid[h].p[i]<p_min && Bid[h].p[i]<P_threshold){
					p_min = Bid[h].p[i];
					kick[0] = h;
					kick[1] = i;
				}
			}
		}
		
		srand(time(0));  //&& rand()%p_min==0
		if(kick[0]>=0 && kick[1]>=0){ //在五个桶中找到了一个可踢的对象

			re = {2,Bid[kick[0]].f[kick[1]],Bid[kick[0]].p[kick[1]],keyfp,kick[1]};

			kick_p += Bid[kick[0]].p[kick[1]];
			kick_t ++;
			if(Bid[kick[0]].p[kick[1]]>kick_max){kick_max=Bid[kick[0]].p[kick[1]];}
			//printf("kick_time:%d  kick_p:%d  kick_max:%d\n",kick_t,kick_p,kick_max);

			Bid[kick[0]].ID[kick[1]] = keyfp;  //replace it
			Bid[kick[0]].f[kick[1]]=1;
			Bid[kick[0]].p[kick[1]]=1;  //继承f和p,不继承就取消掉注释
			Bid[kick[0]].wflag[kick[1]]=1;
			Bid[kick[0]].kick_time[kick[1]]++;  //one kick time

			Bid[kick[0]].den[kick[1]]=(double)(Bid[kick[0]].f[kick[1]])/Bid[kick[0]].p[kick[1]];
			kick_time[kick[0]]++;
			//printf("an item (%d,%d) with f:%d , p:%d is kicked!\n",hash_index,kick,Bid.f[kick],Bid.p[kick]);
	        return re;
		}
		
		//no one can be kick, just ingnore the new one
		//printf("new item with fp:%d insert fail due to fullness!\n",hash_index);
		return re;
		
    }


/*
	PSS_RES insert(uint32_t key, int count = 1) 
		{
			w_cnt++;
			if(w_cnt%WindowSize==0){
				cur_w++;
				for(int i=0;i<bucket_num;i++){
					for(int j=0;j<entry_num;j++){
						Filter[i].wflag[j]=0;
					}	
				}
			}
			//uint32_t keyfp=key;
			L1_ID_type keyfp = key%pow(2,L1_ID_len);
			
			//auto keyfp = getFP(key, ID_len);	 //change the uint32 to a counter_len-int
			PSS_RES re = {-1,0,0,0,0};
	
			// if the item is in B1
			int hash_index = P_Hash->run((const char *)&key, 4) % bucket_num;
			auto &Bid = Filter[hash_index]; //get a bucket in filter
			for (int i = 0; i < entry_num; ++i)  //in this bucket, we have many pos for item
				if (Bid.ID[i] == keyfp)  //find it
				{
					Bid.f[i]++;
					//if(Bid.f[i]>f_max){f_max=Bid.f[i];printf("fmax=%d\n",f_max);}
					assert(Bid.f[i]<pow(2,L1_F_len)-1);  //F overflow
					//if(Bid.w[i]!=cur_w){
					if(Bid.wflag[i]==0){
						Bid.p[i]+=1;
						//if(Bid.p[i]>p_max){p_max=Bid.p[i];printf("pmax=%d\n",p_max);}
						assert(Bid.p[i]<pow(2,L1_P_len)-1);  //F overflow
						Bid.wflag[i]=1;
						//Bid.p[i]+=(Bid.w[i]^cur_w)&0x1;
						//Bid.w[i]=cur_w;
						
						if(Bid.p[i]>P_threshold){  //a pers fow!
							bool report_res = report(key,Bid.f[i],Bid.p[i],Bid.if_L2[i]);  //report to L2
							Bid.if_L2[i] = report_res?1:Bid.if_L2[i];  // res=1: success report.   res=0: fail, do not change
						}
					}
					//B1.permutation(i);  //make it to the first	
	
					re = {1,Bid.f[i],Bid.p[i],keyfp,i};
					Bid.den[i]=(double)(Bid.f[i])/Bid.p[i];
					check_full(hash_index,w_cnt);
					//printf("item (%d,%d) update its info f:%d , d:%d\n",hash_index,i,Bid.f[i],Bid.p[i]);
					return re;
				}
			
			// Not found, not full 
			for (int i = 0; i < entry_num; ++i)
				if (Bid.f[i] == 0)	 //an empty pos
				{
					Bid.ID[i]=keyfp;
					Bid.f[i]=1;
					Bid.p[i]=1;
					//Bid.w[i]=cur_w;
					Bid.wflag[i]=1;
					
					re = {0,1,1,keyfp,i};
					Bid.den[i]=(double)(Bid.f[i])/Bid.p[i];
					//printf("an item (%d,%d) is added\n",hash_index,i);
					return re;
				}
	
			// Not found and full
			uint32_t p_min = (1u<<31)-1;
			int kick = -1;
			for(int i=0;i<entry_num;++i){  //find the min_p (need low than P_thres)
				if(Bid.p[i]<p_min && Bid.p[i]<P_threshold){
					p_min = Bid.p[i];
					kick = i;
				}
			}
			srand(time(0));  //&& rand()%p_min==0
			if(kick>=0 ){ // this item will be kicked out in possibility of 1/p 
	
				re = {2,Bid.f[kick],Bid.p[kick],keyfp,kick};
	
				kick_p += Bid.p[kick];
				kick_t ++;
				if(Bid.p[kick]>kick_max){kick_max=Bid.p[kick];}
				printf("kick_time:%d  kick_p:%d  kick_max:%d\n",kick_t,kick_p,kick_max);
				
	
				Bid.ID[kick] = keyfp;  //replace it
				//Bid.f[kick]=1;
				//Bid.p[kick]=p_min;  //继承f和p,不继承就取消掉注释
				//Bid.w[kick]=cur_w;
				Bid.wflag[kick]=1;
				Bid.kick_time[kick]++;	//one kick time
	
				Bid.den[kick]=(double)(Bid.f[kick])/Bid.p[kick];
				kick_time[hash_index]++;
				//printf("an item (%d,%d) with f:%d , p:%d is kicked!\n",hash_index,kick,Bid.f[kick],Bid.p[kick]);
				return re;
			}
			
			//no one can be kick, just ingnore the new one
			//printf("new item with fp:%d insert fail due to fullness!\n",hash_index);
			return re;
			
		}
*/




	void nextlevel(L2_virtual_PSS * nl)
    {
        L2_vir = nl;
    }


	
};



/***********************************
	Definition for Freq_Sketch (L2)
************************************/

//Bucket_L2: A bucket for FS, it contains only one entry (one fp and one den)
struct Bucket_L2{
    uint32_t fp=0;
    L2_D_type den=0.0;  //density
    int kick_time=0;
};


//class Freq_Sketch: An implement of L2 Sketch, it simplely stores num bucket(entry) for top PS flows
class Freq_Sketch: public L2_virtual_PSS{ 
public:
    Bucket_L2 *PS_record;
	int record_num;
	
	int cur_ptr;
	L2_D_type shared_den;

	uint32_t tmp_w;

    Freq_Sketch() {}
    Freq_Sketch(int _record_num)
        : record_num(_record_num)
    {
        PS_record = new Bucket_L2[record_num];		
		cur_ptr=0;
		shared_den=0.0;
		tmp_w = 0;
    }

	vector<uint32_t> getPI_report(){    //output all PI in L2
		vector<uint32_t> PIList;
		int sum=0;
		for(int i=0;i<cur_ptr;i++){
			if(PS_record[i].den<D_thr){
				PIList.push_back(PS_record[i].fp);
				sum++;
			}
		}
		printf("PI Report:%d  Totally L2:%d\n",sum,cur_ptr);
		return PIList;
	}

	int Rec(){   //output correct num for REC
		int sum=0;
		for(int i=0;i<cur_ptr;i++){
			if(PS_record[i].den<D_thr){sum++;}
		}
		return sum;
	}

	int ARE(){  //output cur_ptr for ARE
		return cur_ptr;
	}

	uint32_t getfp(int pos){
		return PS_record[pos].fp;
	}

	double getden(int pos){
		return PS_record[pos].den;
	}

	int getPI(){
		uint32_t res;
		for(int i=0;i<cur_ptr;i++){
			res = PS_record[i].fp;
		}
		//printf("getPI:%d\n",res);
		return cur_ptr;
	}


#define Burst_Multi 1.1  //how to define burst


    bool insert(uint32_t ID, int f,int p,bool if_kicked) {

		double new_den = (if_kicked&&shared_den>0)?shared_den:((double)f/p);
			
		for(int pos=0;pos<cur_ptr;pos++){  //serach for history
			if(PS_record[pos].fp==ID){
				/*************** Busrt Report/Ignore ***************/
				if((double)f/p>Burst_Multi*PS_record[pos].den){
					//printf("Flow %d Busrt !!!\n",ID);    //burst, report it and do not update
					return 1;
				}
				/***************************************************/
				PS_record[pos].den=(double)f/p;  //find it, just update! nothing about new_den
				if(PS_record[pos].den==0){printf("????\n");}
				return 1;
			}
		}
	
		int new_pos = cur_ptr;   //not found, get a new pos
		if(cur_ptr==record_num){  //the pos is invalid due to fullness
			double den_max = 0.0;
			for(int i=0;i<cur_ptr;i++){
				//printf("%f>%f?\n",PS_record[i].den,den_max);
				if(PS_record[i].den>den_max){
					den_max=PS_record[i].den;  //find the max
					//printf("den_max<-%f\n",PS_record[i].den);
					new_pos = i;
				}
			}
			
			/*
			Insert fail! There're two possible situation:
			1. new_den = den_max -- Two old flow (kicked) using the same share_den and battle for one position. Meaningless op.
			2. new_den > denmax -- New item has a higher den than any of current item, kicked.
			That will make the share_den less and less
			*/
			//printf("newden = %f \t cur_den_max = %f\n",new_den,den_max);
			if(new_den >= den_max){return 0;} 

			//printf("newden = %f \t cur_den_max = %f\n",new_den,den_max);
			shared_den += den_max;  //kicked, update the shared_den
			shared_den /= 2;
			cur_ptr--;
			PS_record[new_pos].kick_time++;
			//printf("Flow %d is kicked and shared_den = (%f+%f)/2 = %f\n",ID,tmp,den_max,shared_den);
			//printf("shared_den changes -> %f\n",shared_den);
		}
		PS_record[new_pos].fp=ID;
		PS_record[new_pos].den=new_den;

		cur_ptr++;
		if(PS_record[new_pos].den==0){printf("!!!!!\n");}
		//printf("Flow %d is stored in %d with den=%f\n",ID,new_pos,new_den);

		return 1;
		
    }

};


/***********************************
	Definition for PSSketch (Main)
************************************/



/*
PSSketch: 
Template: you need to point out what type to use for L1 and L2
Class: Implement of main data structure with 2 level
For convinient, we give L2 to L1 for calling L2 insert()
*/
template<typename L1_Class, typename L2_Class>
class PSSketch{
public:
	P_Filter L1;
	Freq_Sketch L2;

	PSSketch() {}
	PSSketch(const P_Filter& L1_arg, const Freq_Sketch& L2_arg) : L1(L1_arg), L2(L2_arg){
		L1.nextlevel(&L2);
	}

	bool insert(uint32_t key){
		L1.insert(key);
		return 0;
	}
	vector<report_item> getL2PI_report(){   //get all PI in L2 with a vector
		vector<report_item> res;
		vector<uint32_t> reportID = L2.getPI_report();
		for(uint32_t item : reportID){
			report_item rep;
			rep.ID = item;
			//find p in L1
			for(int i=0;i<L1_bucket_num;i++){
				for(int j=0;j<L1_entry_num;j++){
					if(item%pow(2,L1_ID_len)==L1.Filter[i].ID[j]){rep.p=L1.Filter[i].p[j];goto findD;}
				}
			}
			//find d in L2
		findD:
			for(int k=0;k<L2_recourd_num;k++){
				if(item == L2.PS_record[k].fp){rep.D=L2.PS_record[k].den;break;}
			}
			res.push_back(rep);
		}
		return res;
	}
	
	int query(int op){
		if(op==1){  // (no use)	
			printf("No use\n");
		}
		else if(op==2){ //(no use)
			return L2.Rec();
		}
		else if(op==3){
			return L2.ARE();
		}
		else{
			printf("No such op!\n");
		}
		return 0;
	}
	int getL1f(uint32_t key){
		for(int i=0;i<L1.bucket_num;i++){
			for(int j=0;j<L1.entry_num;j++){
				if(key%pow(2,L1_ID_len)==L1.Filter[i].ID[j]){return L1.Filter[i].f[j];}
			}
		}
		return -1;
	}
	int getL1p(uint32_t key){
		for(int i=0;i<L1.bucket_num;i++){
			for(int j=0;j<L1.entry_num;j++){
				if(key%pow(2,L1_ID_len)==L1.Filter[i].ID[j]){return L1.Filter[i].p[j];}
			}
		}
		return -1;
	}
	uint32_t getL2fp(int pos){
		return L2.getfp(pos);
	}
	double getL2den(int pos){
		return L2.getden(pos);
	}
	int getPI(){
		return L2.getPI();
	}
	
	bool clear();
	bool init();

};






#endif
