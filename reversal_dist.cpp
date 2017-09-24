#include <bits/stdc++.h>
#include "reversal_dist.h"
#define endl '\n'

using namespace std;
const int MAXN = (1 << 20);

int n;
vector<int> perm;

void read()
{
    cin >> n;
    perm.assign(n, 0);
    for(int i = 0; i < n; i++)
        cin >> perm[i];
}

string number_suffix(int x)
{
    if(x == 1) return "st";
    if(x == 2) return "nd";
    if(x == 3) return "rd";
    return "th";
}

void solve()
{
    vector<pair<int, int> > li = reversal_distance(perm);

    cout << endl;
    cout << "The reversal distance is equal to "<< li.size() << endl;
    for(int i = 0; i < (int)li.size(); i++)
    {
        cout << "After " << i + 1 << number_suffix((i + 1) % 10) << " reversal: "<< endl;
        cout << "Reverse between " << li[i].first << " and " << li[i].second << endl;
        reverse(perm.begin() + li[i].first - 1, perm.begin() + li[i].second);

        cout << "The permutation currently is: ";
        for(int x: perm) cout << x << " ";
        cout << endl << endl;
    }
}

void comparing_brute_force_and_gentetic_algorithm(int n)
{
	vector<int> perm;
	perm.resize(n);	
	for(int i = 0; i < n; i++) perm[i] = i + 1;

	double avg_ga = 0, avg_bf = 0; 
	for(int runs = 0; runs < 20; runs++)
	{
		random_shuffle(perm.begin(), perm.end());
		avg_ga += reversal_distance(perm).size();
		avg_bf += reversal_distance_brute_force(perm).size();
	}

	cout << "Brute force average:  " << setprecision(6) << fixed << avg_bf / 20.0 << endl;
	cout << "Genetic algorithm average:  " << setprecision(6) << fixed << avg_ga / 20.0 << endl;
}

void comparing_brute_force_and_gentetic_algorithm_maximum_difference(int n)
{
	vector<int> perm;
	perm.resize(n);	
	for(int i = 0; i < n; i++) perm[i] = i + 1;

	int max_difference = 0;
	for(int runs = 0; runs < 20; runs++)
	{
		random_shuffle(perm.begin(), perm.end());
		max_difference = max(max_difference, (int)reversal_distance_brute_force(perm).size() - (int)reversal_distance(perm).size());
	}

	cout << "Maximum difference: " << max_difference << endl;
}


int main()
{
	ios_base::sync_with_stdio(false);
	cin.tie(NULL);

	//read();
	//solve();

	/*
	//Monte Carlo (average of all)
	
	for(int len = 1; len <= 15; len++)
	{
		cout << len << " : " << endl;
		comparing_brute_force_and_gentetic_algorithm(len);
		cout << endl;
	}
	*/

	/*
	//Monte Carlo (maximum difference)
	for(int len = 1; len <= 15; len++)
	{
		cout << len << " : " << endl;
		comparing_brute_force_and_gentetic_algorithm_maximum_difference(len);
		cout << endl;
	}
	*/
	
	return 0;
}
