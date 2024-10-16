
#include"BOBHash32.h"
#include<queue>
#include<unordered_map>
#include<iostream>
#include<map>
#include<ctime>
#include<cstring>
#include<stdlib.h>
#include<vector>
using namespace std;

class BloomFilter{
    int k,m;
    BOBHash32 *bh1,*bh2,*bh3;
    bool *S;
public:
    BloomFilter(int _k,int _m):k(_k),m(_m){
        bh1 = new BOBHash32(11); //seed=11
        bh2 = new BOBHash32(17);
		bh3 = new BOBHash32(23);
        S = new bool[m];
        memset(S,0,sizeof(bool)*m);
    }
    ~BloomFilter(){
        //delete []bh;
        delete []S;
    }
    void insert(uint32_t x){
    	int pos1 = bh1->run((const char*)&x,4) % m;
		int pos2 = bh2->run((const char*)&x,4) % m;
		int pos3 = bh3->run((const char*)&x,4) % m;
        S[pos1]=1;
		S[pos2]=1;
		S[pos3]=1;
    }
    bool query(uint32_t x){
    //P_Hash->run((const char *)&key, 4) % bucket_num;
    	int pos1 = bh1->run((const char*)&x,4) % m;
		int pos2 = bh2->run((const char*)&x,4) % m;
		int pos3 = bh3->run((const char*)&x,4) % m;
		return S[pos1]&&S[pos2]&&S[pos3];
    }
    void clear(){
        memset(S,0,sizeof(bool)*m);
    }
    int getMem(){
        return 1*m;
    }
};

struct Bucket{
    int *w;
	uint16_t *p;
    uint32_t *fp;
    bool *valid;
    int maxnum;
    Bucket(){}
    void init(int maxnum_){
        maxnum = maxnum_;
        valid = new bool[maxnum];
        memset(valid,0,sizeof(bool));
        w = new int[maxnum];
        p = new uint16_t[maxnum];
        fp = new uint32_t[maxnum];
    }
    void insert(uint32_t d,int wi){   //d是当前流，wi是当前窗口的状态，进来时要么-1要么L，且不会改变
        int mi = -1;
        for(int i=0;i<maxnum;++i){
            if(!valid[i]) continue;
            if(fp[i] == d){   //搜索所有的bucket位
                w[i] += wi;    //找到了捏，改变w
                if(wi > 0) p[i] ++;  //-1的时候不加，L的时候加，假阳错误会导致p偏小
                return;
            }
            if(mi < 0 || w[i] < w[mi]) mi = i;   //谁是权重最小的仔
        }
        for(int i=0;i<maxnum;++i){  //走到这里说明没找到f，再遍历一次
            if(!valid[i]){  //找到空位，添加一个新项
                fp[i] = d;
                w[i] = wi;  //假阳错误可能导致这里是-1，代码里并没有处理
                p[i] = 1;  //持续性位1
                valid[i] = 1;
                return;
            }
        }
        w[mi] --;   //走到这里说明既没有找到f，也没有空位，直接惩罚最小的那个仔
       
        if(w[mi] < 0){   //最小的仔寄了，踢掉
            fp[mi] = d;  //上位！
            w[mi] = wi + 1;
            p[mi] = 1;
        }
    }
	
    //pair<int,int> query(uint32_t d){
    //    for(int i=0;i<maxnum;++i){
    //        if(valid[i] && fp[i] == d)
    //            return make_pair(w[i],p[i]);
    //    }
    //    return make_pair(-1,-1);
    //}
	

	
};
class PISketch{
public:
    Bucket *bucket;
    int bucket_num;
    int num_per_bucket;
    BOBHash32 *bh;
    PISketch(int bn,int num_per_bucket_):bucket_num(bn){
        num_per_bucket = num_per_bucket_;
        bucket = new Bucket[bucket_num];
        for(int i=0;i<bucket_num;++i) bucket[i].init(num_per_bucket);
        bh = new BOBHash32(37); //seed=37
    }
    void insert(uint32_t d,int wi){
        int pos = bh->run((const char*)&d,4) %bucket_num;
        bucket[pos].insert(d,wi);
    }
    //pair<int,int> query(uint32_t d){
    //    int pos = bh.run((const char*)d,4) %bucket_num;
    //    return bucket[pos].query(d);
    //}

	uint32_t getfp(int x,int y){
		return bucket[x].fp[y];		
	}

	int getw(int x,int y){
		return bucket[x].w[y];		
	}
	
    //int getMem(){
    //    return bucket_num * num_per_bucket * ( sizeof(int) + sizeof(int) + sizeof(Data) + sizeof(bool) );
    //}
};