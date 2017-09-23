#include <bits/stdc++.h>
#include "Signed_reversal_distance/invdist.c"

using namespace std;
const int MAXN = 242;
const int PROB = 136;                                                                /// 1/136

int n;
vector<int> perm; //the permutation which reversal distance (RD) we are computing
double av_rd=0.0; //the average RD after the runs of the program

void read()                                                                          ///read
{
    cin >> n; //the length of the permuatation
}

void init()                                                                          ///init
{
    perm.clear();
    perm.resize(n);
    for (int i=0; i<n; ++i) perm[i]=i+1; //assigning the permutation 1, 2, 3, ...

    random_shuffle(perm.begin(), perm.end()); //getting a random permutation
}

vector<int> gen_random_signed_permutation(vector<int> perm)                          ///get random signs for the elements of the permutation
{
    vector<int> res;
    for(int i = 0; i < (int)perm.size(); i++)
    {
        int prob = (rand() % 2); //1/2 probability for the sign to be either - or +
        if(prob) res.push_back(-perm[i]);
        else res.push_back(perm[i]);
    }

    return res;
}

vector<int> mutate_permutation(vector<int> p, int prob)                              ///mutate
{
    if(rand() % prob) return p;

    int idx = rand() % p.size(); //small probability for mutation
    p[idx] *= -1; //changing the sign of the element, which we are mutating

    return p;
}

vector<int> combine(vector<int> l, vector<int> r)                                   ///combine
{
    vector<int> res; //the permutation-child of the parents l and r
    int v1 = rand() % r.size(); //a random point of crossing is choosen

    for(int i = 0; i < v1; i++) res.push_back(l[i]); //the first part of the resulting permutation we take from the first parrent
    for(int i = v1; i < (int)r.size(); i++) res.push_back(r[i]); //the second part of the resulting permutation we take from the second parrent

    return res;
}

int reversal_distance(vector<int> perm)                                             ///reversal distance
{
    struct genome_struct g1;
    struct genome_struct g2;

    g1.genes = ( int * ) malloc ( ( n + 2 ) * sizeof ( int ) );
    g2.genes = ( int * ) malloc ( ( n + 2 ) * sizeof ( int ) );

    for(int i = 0; i < n; i++) //the permuation
        g1.genes[i] = perm[i];

    for(int i = 0; i < n; i++) //the identical permuation
        g2.genes[i] = i + 1;

    int result = invdist_noncircular_nomem(&g1, &g2, n, n); //calculating the RD betweeen the permutation and its identical permutation

    free ( g1.genes );
    free ( g2.genes );

    return result;
}

vector<vector<int> > listt; //with the permutations
int dist[2 * MAXN * MAXN]; //just the distances

double get_rand01()
{
    double low = 0, high = 1;
    for(int ops_bin_search = 0; ops_bin_search < 20; ops_bin_search++)
        if(rand() % 2) low = (low + high) / 2.0;
        else high = (low + high) / 2.0;

    return high;
}

int roulette(vector<double> pref)                                                      ///roulette wheel method for the selection
{
    double rand_el=get_rand01(); //selecting random integer between 0 and 1

    int ind=pref.size()-1;
    for(int i=0; i<(int)pref.size()-1; i++)
    {
        if(pref[i]<=rand_el&&rand_el<pref[i+1]) ind=i; //finding the index corresponding to it
    }

    return ind;
}

void solve()                                                                        ///MAIN PART:
{
    listt.clear();
	for(int i = 0; i < MAXN*MAXN; i++) dist[i]=0;

    for(int i = 0; i < n * n; i++) //the array with the random-signed permutations
    {
        vector<int> tmp = gen_random_signed_permutation(perm);
        listt.push_back(tmp);
    }

    int cnt_same = 1; //within how many generations the reversal distance has not been changed
    double result=(double)MAXN; //the smallest reversal distance of those of the permuations from tmp

    while(true)
    {
        double nw_mn = result; //for comparing the smallest previously found reversal distance with the current one
        for(int i = 0; i < n * n; i++)
        {
            dist[i]= reversal_distance(listt[i]); //finding the signed reversal distance of the current permutation
            nw_mn=min((double)dist[i], nw_mn); //updating the current optimal RD if a smaller value was found
        }

        if(nw_mn < result)
        {
            result = nw_mn;
            cnt_same = 1; //a new best RD was found!
        }
        else cnt_same++; //another generation without more optimal reversal distance

        if(cnt_same == 3) //the smallest RD has not been changed within three generations - end of the algorithm
            break;

        vector<pair<int, vector<int> > > curr_li; //pair<RD, permuatation with RD, equal to that>
        vector<vector <int> > new_li; //the new generation

        double summ=0.0; //for normalization
        for(int i = 0; i < n * n; i++)
        {
            curr_li.push_back(make_pair(dist[i], listt[i]));
            summ+=(double)1/dist[i];
        }

        vector<double> fitness; //the fitness function for every individual
        for(int i = 0; i < n * n; i++)
        {
            fitness.push_back((double)((double)1/dist[i])/summ);
        }

        vector<double> pref; //preparing the array for the selection
        double prev=0.0;
        pref.push_back(prev);
        for(int i = 0; i < n * n; i++)
        {
            prev+=fitness[i];
            pref.push_back(prev);
        }

        for(int i = 0; i < n * n; i++) //combining using the roulette wheel method
        {
            int l_in=roulette(pref), r_in=roulette(pref);
            new_li.push_back(combine(curr_li[l_in].second, curr_li[r_in].second));
        }

        for(int i = 0; i < (int)new_li.size(); i++) //loop for mutating the new generation
            new_li[i] = mutate_permutation(new_li[i], PROB);                          ///mutation 1/PROB

        listt = new_li;
    }

    av_rd += result;
}

int main()
{
    srand(time(0));

    for(int i = 0; i < 15; i++) //testing for 15 values (could be more or less) - just an easier way for the testing
    {
        av_rd = 0; //summing the reversal distances to found the mean of 10 testings (to avoid working only with "bad" or "good" permuations) - for finding the mean value
        read(); //the length of the permutation
        for(int j = 0; j < 10; j++)
        {
            init(); //getting the random permutation which RD we are expecting
            solve(); //finding the reversal distance
        }
        cout << "Approximate reversal distance: " << av_rd / 10 << endl; //the mean value of all RDs for the testings
    }

    return 0;
}
