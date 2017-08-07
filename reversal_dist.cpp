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
    for(int i = 0; i < li.size(); i++)
    {
        cout << "After " << i + 1 << number_suffix((i + 1) % 10) << " reversal: "<< endl;
        cout << "Reverse between " << li[i].first << " and " << li[i].second << endl;
        reverse(perm.begin() + li[i].first - 1, perm.begin() + li[i].second);

        cout << "The permutation currently is: ";
        for(int x: perm) cout << x << " ";
        cout << endl << endl;
    }
}

int main()
{
	ios_base::sync_with_stdio(false);
	cin.tie(NULL);

	read();
	solve();
	return 0;
}
