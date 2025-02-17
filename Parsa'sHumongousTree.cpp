#include <bits/stdc++.h>
#include <ext/pb_ds/assoc_container.hpp>
#include <ext/pb_ds/tree_policy.hpp>
using namespace std;
using namespace __gnu_pbds;

#define int long long
#define ll long long
#define ld long double
#define double long double
#include <cstring>
int min(int a,int b){ if(a>b) return b; return a;}
int max(int a,int b){ if(a>b) return a; return b;}
// ll max(ll a,ll b){ if(a>b) return a; return b;}int min(int a,int b){ if(a>b) return b; return a;}

constexpr ll INF = 4e18;
constexpr double EPS = 4e-218;
#define eps  1e-8
// constexpr int mod = 998244353;
#define mod 1000000007
#define inf 4e18;




bool cmpsecond(pair<ll, ll>& a,pair<ll, ll>& b) { return a.second < b.second; }
bool cmpfirst(pair<ll, ll>& a,pair<ll, ll>& b) { return a.first < b.first; }
bool cmp(pair<ll, ll>& a,pair<ll, ll>& b) { if(a.second!=b.second) return a.second < b.second; else return a.first<b.first; }



#define maxpq priority_queue<int> // Default max-heap
#define minpq priority_queue<int, vector<int>, greater<int>> // Min-heap
typedef vector<int> vi;
typedef vector<vector<int>> vvi;
typedef vector<ll> vll;
typedef vector<vector<ll>> vvll;
typedef vector<string> vs;
typedef vector<vector<string>> vvs;
#define all(x) x.begin(), x.end()
#define rall(x) x.rbegin(), x.rend()
#define from(vec,a,b) vec.begin() + a, vec.begin() + b+1
#define mp make_pair
#define pb push_back
#define nline '\n'
#define yes cout << "YES\n"
#define no cout << "NO\n"
#define sqrt (long long)sqrtl
#define f for(int i=0;i<n;i++)
#define fi(a,b) for (int i = a; i < b; i++)
#define fj(a,b) for (int j = a; j < b;j++)
#define fk(a,b) for (int k = a; k < b;k++)
#define vvvi(name,a,b,c) vector<vector<vector<int>>> name(a+1, vector<vector<int>>(b, vector<int>(c, 0)));
vi ask(int i, int j)
{
    cout << "? " << i << " " << j << endl;
    vi temp;
    fk(0,(j-i+1)) { int x; cin>>x; temp.push_back(x);}
    return temp;
}
class dsu
{
private:
    vector<int> par;
    vector<int> sz;
 
public:
    dsu(int n)
    {
        par = vector<int>(n);
        iota(par.begin(), par.end(), 0);
        sz = vector<int>(n, 1);
    }
 
    int parent(int u)
    {
        // this optimisation was good.
        if (par[u] != par[par[u]])
            par[u] = parent(par[par[u]]);
        return par[u];
    }
 
    bool connected(int u, int v)
    {
        u = parent(u);
        v = parent(v);
        if (u == v)
            return true;
        return false;
    }
 
    bool join(int u, int v)
    {
        u = parent(u);
        v = parent(v);
        if (u == v)
            return false;
        if (sz[u] <= sz[v])
        {
            sz[v] += sz[u];
            par[u] = v;
        }
        else
        {
            sz[u] += sz[v];
            par[v] = u;
        }
        return true;
    }
};
 

// ll inv(ll i) {if (i == 1) return 1; return (mod - ((mod / i) * inv(mod % i)) % mod) % mod;}
 
// ll modmul(ll a, ll b) {a = a % mod; b = b % mod; return (((a * b) % mod) + mod) % mod;}
 
// ll modadd(ll a, ll b) {a = a % mod; b = b % mod; return (((a + b) % mod) + mod) % mod;}
 
// ll modsub(ll a, ll b) {a = a % mod; b = b % mod; return (((a - b + mod) % mod) + mod) % mod;}
  
// ll ceildiv(ll a, ll b) {return a % b == 0 ? a / b : a / b + 1;}
 
// ll pow(ll a, ll b) {a %= mod; ll res = 1; while (b > 0) {if (b & 1) res = res * a % mod; a = a * a % mod; b >>= 1;} return res;}
 
// vector<ll> sieve(int n) {int*arr = new int[n + 1](); vector<ll> vect; for (int i = 2; i <= n; i++)if (arr[i] == 0) {vect.push_back(i); for (int j = 2 * i; j <= n; j += i)arr[j] = 1;} return vect;}
int inv(int i) {if (i == 1) return 1; return (mod - ((mod / i) * inv(mod % i)) % mod) % mod;}
 
int modmul(int a, int b) {a = a % mod; b = b % mod; return (((a * b) % mod) + mod) % mod;}
 
int modadd(int a, int b) {a = a % mod; b = b % mod; return (((a + b) % mod) + mod) % mod;}
 
int modsub(int a, int b) {a = a % mod; b = b % mod; return (((a - b + mod) % mod) + mod) % mod;}
  
int ceildiv(int a, int b) {return a % b == 0 ? a / b : a / b + 1;}
 
int pow(int a, int b) {a %= mod; ll res = 1; while (b > 0) {if (b & 1) res = res * a % mod; a = a * a % mod; b >>= 1;} return res;}
vi fact(1000005);
vector<int> sieve(int n) {vector<int> vect; for (int i = 2; i <= n; i++) if (fact[i] == 0) {vect.push_back(i); for (int j = i; j <= n; j += i) if(fact[j]==0) fact[j] = i; } return vect;}

// Mathematical functions
int gcd(int a, int b)
{
    while (b)
    {
        a %= b;
        swap(a, b);
    }
    return a;
}
 
int GCD_extended(int a, int b, int &x, int &y)
{
    x = 1, y = 0;
    int x1 = 0, y1 = 1, a1 = a, b1 = b;
    while (b1)
    {
        int q = a1 / b1;
        tie(x, x1) = make_tuple(x1, x - q * x1);
        tie(y, y1) = make_tuple(y1, y - q * y1);
        tie(a1, b1) = make_tuple(b1, a1 - q * b1);
    }
    return a1;
}
int lcm(int a, int b)
{
    return ((ll)a * b) / gcd(a, b);
}
 
ll modpow(ll x, ll n, int m = mod)
{
    if (x == 0 && n == 0)
        return 0; // undefined case
    ll res = 1;
    while (n > 0)
    {
        if (n % 2)
            res = (res * x) % m;
        x = (x * x) % m;
        n /= 2;
    }
    return res;
}
 
int modinv(int x, int m = mod)
{
    return modpow(x, m - 2, m);
}
mt19937 rng;
int randomnumber(int l, int r)
{
    uniform_int_distribution<int> dist(l, r);
    return dist(rng);
}

 
int reverse(int n)
{
    int r=0;
    while(n>0)
    {
        r=r*10+n%10;
        n/=10;
    }
    return r;
}

bool ispalindrome(int n)
{
    return (reverse(n)==n); 
}
#define py cout<<"Yes"<<endl
#define pn cout<<"No"<<endl
int binarytodecimal(string s) { return bitset<64>(s).to_ullong(); }
string decimaltobinary(int a) { return bitset<64>(a).to_string(); }
 
ll andOperator(ll a, ll b)
{
    ll shiftcount = 0;
 
    while (a != b and a > 0)
    {
        shiftcount++;
        a = a >> 1;
        b = b >> 1;
    }
    return int64_t(a << shiftcount);
}
int factorial(int n){
    if (n==0){
        return 1;
    }
    ll ans=1;
    for (ll i=1;i<=n;i++){
         ans=modmul(ans,i);
    }
    return ans;
}
 
 
 
 
int power(int base, int exp)
{
    if (exp == 0)
       return 1;
    else if (exp == 1)
       return base;
    else
    {
       int calc;
       if (exp % 2 == 0)
       {
         calc = power(base, exp/2);
         calc *= calc;
       }
       else
       {
         calc = base*power(base, exp-1);
       }
       return calc;
    }
}
 
// int sqrt(ll n){
//     ll ans=0;
//     ll low=1;
//     ll high=1e9;
//     while(low<=high){
//         ll md=(low+high)/2;
//         if (md*md<=n){
//             ans=md;
//             low=md+1;
//         }
//         else{
//             high=md-1;
//         }
//     }
//     return ans;
// }
bool parity(int a,int b){
    if(a%2==0 && b%2==0) return true;
    if(a%2!=0 && b%2!=0) return true;
    return false;
}


struct Mint
{
    int value;
    static const int MOD_value = mod;
 
    Mint(long long v = 0)
    {
        value = v % mod;
        if (value < 0)
            value += mod;
    }
    Mint(long long a, long long b) : value(0)
    {
        *this += a;
        *this /= b;
    }
 
    Mint &operator+=(Mint const &b)
    {
        value += b.value;
        if (value >= mod)
            value -= mod;
        return *this;
    }
    Mint &operator-=(Mint const &b)
    {
        value -= b.value;
        if (value < 0)
            value += mod;
        return *this;
    }
    Mint &operator*=(Mint const &b)
    {
        value = (long long)value * b.value % mod;
        return *this;
    }
 
    friend Mint mexp(Mint a, long long e)
    {
        Mint res = 1;
        while (e)
        {
            if (e & 1)
                res *= a;
            a *= a;
            e >>= 1;
        }
        return res;
    }
    friend Mint inverse(Mint a) { return mexp(a, mod - 2); }
 
    Mint &operator/=(Mint const &b) { return *this *= inverse(b); }
    friend Mint operator+(Mint a, Mint const b) { return a += b; }
    friend Mint operator-(Mint a, Mint const b) { return a -= b; }
    friend Mint operator-(Mint const a) { return 0 - a; }
    friend Mint operator*(Mint a, Mint const b) { return a *= b; }
    friend Mint operator/(Mint a, Mint const b) { return a /= b; }
    friend std::ostream &operator<<(std::ostream &os, Mint const &a) { return os << a.value; }
    friend bool operator==(Mint const &a, Mint const &b) { return a.value == b.value; }
    friend bool operator!=(Mint const &a, Mint const &b) { return a.value != b.value; }
};
Mint ncr(int n, int r)
{
    assert(n >= r);
    assert(r >= 0);
    Mint a = 1;
    for (int i = n; i >= n - r + 1; i--)
    {
        a *= i;
    }
    Mint b = 1;
    for (int i = 1; i <= r; i++)
    {
        b *= i;
    }
    return a / b;
}

int leftshift(int i){
    if(i==0) return 1;
    int ans=1;
    for(int j=1;j<=i;j++){
        ans=ans*2;
    }
    return ans;
}

int customAND(int a, int b) {
    int result = 0, power = 1;
    while (a > 0 || b > 0) {
        if ((a % 2 == 1) && (b % 2 == 1)) {
            result += power;
        }
        a /= 2;
        b /= 2;
        power *= 2;
    }
    return result;
}

// Custom right shift function
long long customRightShift(long long num, int shift) {
    while (shift--) {
        num /= 2;
    }
    return num;
}

// Function to count set bits using custom AND and Right Shift
int countsetbits(long long num) {
    int count = 0;
    while (num > 0) {
        if (customAND(num, 1) == 1) { // Check if LSB is 1
            count++;
        }
        num = customRightShift(num, 1); // Right shift by 1
    }
    return count;
}
bool sendd(int n){
    while(n>0){
        int x=n%10;
        if(x==7) return true;
        n/=10;
    }
    return false;
}
int bfs(int n) {
    std::queue<int> q;
    // minpq q;
// if(sendd(node)) return 0;
    // visited[i] = 1;
    q.push(n);
    vector<int>v(10);
    int x=9;
    for(int i=1;i<11;i++){
        v[i-1]=x;
        x*=10;
        x+=9;
    }
    // for(auto i:v) cout<<i<<" ";
    
    int count=0;
    map<int,int>depth;
    depth[n]=0;
    while (!q.empty()) {

        int node = q.front();
        if(sendd(node)) return depth[node];
        // cout<<node<<" ";
        
        q.pop();
        // if(node>1e18) continue;
        count++;
        for (int i:v) {
            // count++;
            // if (s.find(node+i)==s.end()) {
            //     q.push(node+i);
            //     depth[node+i]=1+depth[node];
            // }
        }
    }
    return -1;
}


struct CompareSum {
    bool operator()(const pair<int, int>& a, const pair<int, int>& b) {
        return (a.first) < (b.first );
    }
};

int digitsum(int n) {
    int sum = 0;
    while (n > 0) {
        sum += n % 10; // Extract the last digit and add to sum
        n /= 10;       // Remove the last digit
    }
    return sum;
}



int dfs(int n,vector<pair<int,int>>&cost,vvi &v,int root,int p,vvi &dp,int parent) {
    // std::cout << root << " ";
    if(dp[root][p]!=-1) return dp[root][p];
    // if(p==1)
    // visited[root] = 1;
    int child=0,temp=0;
    int ans=0;
    int maxi=0;
    if(p==0){
        maxi=cost[root].first;
    }
    if(p==1){
        maxi=cost[root].second;
    }
    for (int neighbor : v[root]) {
        
        if (neighbor!=parent) {
            child++;
            // vi temp=visited;
            // case if i treat root as p
            int ans1=dfs(n,cost,v,neighbor,0,dp,root);

            int ans2=dfs(n,cost,v,neighbor,1,dp,root);
            // if(root==6 && p==0){
            //     cout<<endl;
            //     cout<<ans1<<" "<<ans2<<" hii"<<endl;
            // }

            ans+=(max(abs(maxi-cost[neighbor].first)+ans1,abs(maxi-cost[neighbor].second)+ans2));
            if(root==6 && p==1){
                // cout<<endl;
                // cout<<ans1<<" "<<ans2<<" hii "<<ans<<endl;
            }
        }
    }
    // if(child==0){
    //     dp[root][p]=maxi;
    //     return maxi;
    // }
    // if(p==0) visited[root]=0;
    return dp[root][p]=ans;
}
void solve(){
 int n;
 cin>>n;
 vector<pair<int,int>>cost(n+1);
 fi(1,n+1){
    int x,y;
    cin>>x>>y;
    cost[i]=mp(x,y);
 }
 vvi v(n+1);
 fi(0,n-1){
    int x,y;
    cin>>x>>y;
    v[x].pb(y);
    v[y].pb(x);
 }
 vi visited(n+1);
 vvi dp(n+1,vi (2,-1));
 // memset(dp,-1,sizeof,dp(dp,-1));
 int ans1=dfs(n,cost,v,1,0,dp,0);
    fi(0,n+1){
        visited[i]=0;
    }
    // if(n==6)
    // cout<<dp[6][0]<<" this is me ";
    // fi(0,n+1){
    //     dp[i][0]=-1;
    //     dp[i][1]=-1;
    // }

    int ans2=dfs(n,cost,v,1,1,dp,0);
    // int ans2=0;
    // for(int )
    // cout<<ans1<<" "<<ans2<<endl;
    cout<<max(ans1,ans2)<<endl;


    



}



signed main() {
    ios_base::sync_with_stdio(false);
    cin.tie(nullptr);


    int x=0;
    int t=1;
    cin >> t;
    // actually my answer is already fixed so i can find it here and then pass the answer

    

    while (x<t) {
        solve();
        x++;
    }

    return 0;
}
