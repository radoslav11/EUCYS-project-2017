#include <bits/stdc++.h>
#include "invdist.c"

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
    for(int i = 0; i < perm.size(); i++)
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
    for(int i = v1; i < r.size(); i++) res.push_back(r[i]); //the second part of the resulting permutation we take from the second parrent

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

int count_mns[MAXN];
vector<vector <int> > li[MAXN];
double avg[MAXN];
int dist[MAXN][MAXN * MAXN];

void solve()                                                                        ///MAIN PART:
{
	for(int i = 0; i < MAXN; i++) li[i].clear(); //clearing everything - in case of a previous test
	for(int i = 0; i < MAXN; i++) avg[i] = 0;
	for(int i = 0; i < MAXN; i++)
        for(int j = 0; j < MAXN*MAXN; j++) dist[i][j] = 0;

    for(int i = 0; i < n * n; i++) //the array with the random-signed permutations
    {
        vector<int> tmp = gen_random_signed_permutation(perm);

        int cnt_m = 0;
        for(int j = 0; j < n; j++)
            cnt_m += tmp[j] < 0;

        li[cnt_m].push_back(tmp); //a list with the permutations with every count of negative signs (cnt_m)
    }

    int cnt_same = 1; //within how many generations the reversal distance has not been changed
    double result=(double)MAXN; //the smallest reversal distance of those of the permuations from tmp

    while(true)
    {
        double nw_mn = result; //for comparing the smallest previously found reversal distance with the current one
        for(int i = 0; i <= n; i++)
            for(int j = 0; j < li[i].size(); j++)
            {
                dist[i][j] = reversal_distance(li[i][j]); //finding the reversal distance of a permutation with i negative signs
                nw_mn = min(nw_mn, (double)dist[i][j]); //updating the current optimal RD if a smaller value was found
            }

        if(nw_mn < result)
        {
            result = nw_mn;
            cnt_same = 1; //a new best RD was found!
        }
        else cnt_same++; //another generation without more optimal reversal distance

        if(cnt_same == 3) //the smallest RD has not been changed within three generations - end of the algorithm
            break;

        for(int c = 0; c <= n; c++) //computing the RDs of the permutations with the different counts of negative signs
        {
            int sum = 0;
            for(int i = 0; i < li[c].size(); i++)
                sum += dist[c][i];

            avg[c] = (double)sum / (double)li[c].size(); //the average value of all the RDs of permutations with this count of negative signs
        }

        vector<pair <double, vector <int> > > curr_li; //pair<average distance of all permuations with i negative signs, permuatation with i negative signs>
        vector<vector <int> > new_li;

        for(int i = 0; i <= n; i++)
        {
            for(int j = 0; j < li[i].size(); j++)
                curr_li.push_back({avg[i], li[i][j]});
            li[i].clear();
        }

        sort(curr_li.begin(), curr_li.end()); //for finding the number of negative signs for which the mean of the RDs of permutations of this count of negative signs is the smallest
        for(int i = 0; i < n; i++) //combining the best individuals to form the new generation
        {
            for(int j = 0; j <= i; j++)
            {
                new_li.push_back(combine(curr_li[j].second, curr_li[i].second));
                if(new_li.size() >= n * n) break; //if the size of the population is already n*n, the generaton is completely filled in
            }
            if(new_li.size() >= n * n) break;
        }

        for(int i = 0; i < new_li.size(); i++) //loop for mutating the new generation
            new_li[i] = mutate_permutation(new_li[i], PROB);                          ///mutation 1/PROB

        for(int i = 0; i < new_li.size(); i++)
        {
            int cnt_m = 0;
            for(int j = 0; j < n; j++)
                cnt_m += new_li[i][j] < 0;

            li[cnt_m].push_back(new_li[i]); //a list with the permuatatios with every count of negative signs (cnt_m)
        }
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
