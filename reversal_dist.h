#ifndef REVERSAL_DIST_H
#define REVERSAL_DIST_H

#include <bits/stdc++.h>
#include "Signed_reversal_distance/invdist.c"

using namespace std;

const int RD_MAXN = 242;
const int RD_PROB = 25;

/** all functions find the reversal distance from a permutation to it's identity permutation (sorted reordering of the elements) **/

int signed_reversal_distance(vector<int> perm);
vector<pair<int, int> > reversal_distance(vector<int> perm);
vector<int> mutate_permutation(vector<int> p, int RD_PROB);
vector<int> combine_permutations(vector<int> l, vector<int> r);
vector<int> generate_random_signed_permutation(vector<int> perm);
void find_signed_order(vector<int> perm, vector<pair<int, int> > &answer);
bool is_sign_sorted(vector<int> perm);




/* a function that generates a random sign for every element of a permutation */
vector<int> generate_random_signed_permutation(vector<int> perm)
{
    vector<int> res;
    for(int i = 0; i < (int)perm.size(); i++)
    {
        int RD_PROB = (rand() % 2); //1/2 RD_PROBability for the sign to be either - or +
        if(RD_PROB) res.push_back(-perm[i]);
        else res.push_back(perm[i]);
    }

    return res;
}

/* a function to find the signed reversal distance of a permutation in linear time */
int signed_reversal_distance(vector<int> perm)
{
    int n = perm.size();
    struct genome_struct g1;
    struct genome_struct g2;

    g1.genes = ( int * ) malloc ( ( n + 2 ) * sizeof ( int ) );
    g2.genes = ( int * ) malloc ( ( n + 2 ) * sizeof ( int ) );

    for(int i = 0; i < n; i++) //the permuation
        g1.genes[i] = perm[i];

    for(int i = 0; i < n; i++) //the identical permuation
        g2.genes[i] = i + 1;

    int result = invdist_noncircular_nomem(&g1, &g2, 0, n); //calculating the RD betweeen the permutation and its identical permutation

    free ( g1.genes );
    free ( g2.genes );

    return result;
}

/* function to check if a signed permutation is sorted */
bool is_sign_sorted(vector<int> perm)
{
    for(int i = 0; i < (int)perm.size(); i++)
        if((perm[i] > 0 ? perm[i] : -perm[i]) != i + 1)
            return false;

    return true;
}

/* function to generate the reversal operations */
void find_signed_order(vector<int> perm, vector<pair<int, int> > &answer)
{
    if(is_sign_sorted(perm)) return;

    int current_reversal_dist = signed_reversal_distance(perm);

    for(int l = 0; l < (int)perm.size(); l++)
        for(int r = l; r < (int)perm.size(); r++)
        {
            for(int i = l; i <= r; i++) perm[i] *= -1;
            reverse(perm.begin() + l, perm.begin() + r + 1);

            if(signed_reversal_distance(perm) + 1 == current_reversal_dist)
            {
                answer.push_back({l + 1, r + 1});
                find_signed_order(perm, answer);
                return;
            }

            for(int i = l; i <= r; i++) perm[i] *= -1;
            reverse(perm.begin() + l, perm.begin() + r + 1);
        }
}

/* a function that finds the list of reversals to sort a permutation */
vector<pair<int, int> > reversal_distance(vector<int> permutation)
{
    int n = permutation.size();
    vector<vector<int> > listt;
    vector<int> dist;

    if(n <= 1) { return vector<pair<int, int> >(0, {0, 0}); }

    dist.assign(RD_MAXN * RD_MAXN, 0);
	listt.clear();

    for(int i = 0; i < n * n; i++) //the array with the random-signed permutations
    {
        vector<int> tmp = generate_random_signed_permutation(permutation);
        listt.push_back(tmp);
    }

    int cnt_same = 1; //within how many generations the reversal distance has not been changed
    int result = RD_MAXN; //the smallest reversal distance of those of the permuations from tmp
    vector<int> best_permutation, current_best_permutation;

    while(true)
    {
        int nw_mn = result; //for comparing the smallest previously found reversal distance with the current one
        for(int i = 0; i < n * n; i++)
        {
            dist[i] = signed_reversal_distance(listt[i]); //finding the signed reversal distance of the current permutation
            if(dist[i] < nw_mn) current_best_permutation = listt[i];
            nw_mn = min(dist[i], nw_mn); //updating the current optimal RD if a smaller value was found
        }

        if(nw_mn < result)
        {
            result = nw_mn;
            best_permutation = current_best_permutation;
            cnt_same = 1; //a new best RD was found!
        }
        else cnt_same++; //another generation without more optimal reversal distance

        if(cnt_same == 3) //the smallest RD has not been changed within three generations - end of the algorithm
            break;

        vector<pair<int, vector<int> > > curr_li; //pair<RD, permuatation with RD, equal to that>
        vector<vector <int> > new_li; //the new generation
        for(int i = 0; i < n * n; i++)
            curr_li.push_back(make_pair(dist[i], listt[i]));

        sort(curr_li.begin(), curr_li.end()); //for finding the number of negative signs for which the mean of the RDs of permutations of this count of negative signs is the smallest
        for(int i = 0; i < n; i++) //combining the best individuals to form the new generation
        {
            for(int j = 0; j <= i; j++)
            {
                new_li.push_back(combine_permutations(curr_li[j].second, curr_li[i].second));
                new_li.push_back(combine_permutations(curr_li[i].second, curr_li[j].second));
                if((int)new_li.size() >= n * n) break; //if the size of the population is already n*n, the generaton is completely filled in
            }
            if((int)new_li.size() >= n * n) break;
        }

        for(int i = 0; i < (int)new_li.size(); i++) //loop for mutating the new generation
            new_li[i] = mutate_permutation(new_li[i], RD_PROB);                          ///mutation with proability equal to 1/RD_PROB

        listt = new_li;
    }

    vector<pair<int, int> > answer;
    find_signed_order(best_permutation, answer);
    return answer;
}

/* the mutation is just changing the sign of one random element of the permutation with probability 1/RD_PROB */
vector<int> mutate_permutation(vector<int> p, int RD_PROB_F)                              ///mutate
{
    if(rand() % RD_PROB_F) return p;
    int idx = rand() % p.size();
    p[idx] *= -1;
    return p;
}

/* in the combining of two individuals we get a random position from 1 to n (lets call it v1) and then we merge the first v1 elements of the first permutation and last n - v1 of the second one */
vector<int> combine_permutations(vector<int> l, vector<int> r)                                   ///combine
{
    vector<int> res;
    int v1 = rand() % r.size();
    for(int i = 0; i < v1; i++) res.push_back(l[i]);
    for(int i = v1; i < (int)r.size(); i++) res.push_back(r[i]);

    return res;
}

#endif
