#include <bits/stdc++.h>
#include <ext/pb_ds/assoc_container.hpp>
#include <ext/pb_ds/tree_policy.hpp>
using namespace std;
using namespace __gnu_pbds;

#define int long long
#define ll long long
#define ld long double
#define double long double

int min(int a,int b){ if(a>b) return b; return a;}
int max(int a,int b){ if(a>b) return a; return b;}
// ll max(ll a,ll b){ if(a>b) return a; return b;}int min(int a,int b){ if(a>b) return b; return a;}

constexpr ll INF = 4e18;
constexpr double EPS = 4e-218;
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

void bfs(int i, vvi &graph) {
    int n=graph.size()-1;
    vi visited(n + 1, 0);  // Mark all nodes as unvisited
    // std::queue<int> q;
    minpq q;

    visited[i] = 1;
    q.push(i);

    while (!q.empty()) {

        int node = q.top();
        cout<<node<<" ";
        q.pop();

        for (int neighbor : graph[node]) {
            if (!visited[neighbor]) {
                visited[neighbor] = 1;
                q.push(neighbor);
            }
        }
    }
}


struct CompareSum {
    bool operator()(const pair<int, int>& a, const pair<int, int>& b) {
        return (a.first) < (b.first );
    }
};
int ans=0;
int dfs(int node, vvi &graph, vi& visited,int parent) {
    // std::cout << node << " ";
    visited[node] = 1;
    // now i am at a node and i will check if i can break the edge with the child
    int child=0,temp=0;
    // cout<<node<<" ";
    for (int neighbor : graph[node]) {
        if (!visited[neighbor] ) {
            child++;
            int lll=dfs(neighbor, graph, visited,node);
            // temp==0 // it means vo kisi ke sath set ho gya
            temp+=lll;
            if(lll==0){
                ans++;
                child--;   }

        }
    }
    if(child>1 && (-3*child)== (temp)){
        // kuch bnde mile h
        // ye bnde 
        if(child%2!=0) return 0;
        else return -3;
    }
    if(child==1 && temp==-3) return 0;
    if(child==0) {
        return -3;
    }
    return 1;


    // if(dp[node]==0) dp[node]=-1;
    // if(child%2 ==0) {ans++; return 1;}
    // else return 1;
    
}
bool send(int n,vi &v,int mid,int k,int i){
    // now i will try to make ith position as mid
    if(v[i]>=mid) return true;
    int count=0,temp=mid;
    for(int j=i;j<n;j++){
        if(v[j]<temp){
            count+=(temp-v[j]);
            temp--;
            if(count>k) return false;
        }
        else return true;
    }
    return false;
}
void fillbinary(vi &v,int n){
    for(int i=72;i>=0;i--){
        if((leftshift(i)&n)) v[i]=1;
    }
}

void solve(){

    int n,m;
    cin>>n>>m;
    vi a(n),b(m);
    fi(0,n) cin>>a[i];
    fi(0,m) cin>>b[i];
    // fi(0,n) cout<<a[i]<<" ";
    // cout<<endl;
    sort(b.begin(),b.end());
    // now at every pos
    // a.push_back(1e13);
    for(int i=n-1;i>=0;i--){
        // if(i+1<n)
        // cout<<a[i+1]<<" ";
        // i will try to maximise v[i] such that it is less than v[i+1]
        int x=1e13;
        if(i!=n-1)
        x=(a[i]+a[i+1]);
        // cout<<x<<" ";
        // now i need x 
        int y=upper_bound(b.begin(),b.end(),x)-b.begin();
        y--;
        if(y==-1) continue;
        if(y==m-1) {
            a[i]=max(a[i],b[y]-a[i]);
        }
        // if(b[y+1]==x) y++;
        // cout<<y<<endl;
        if(i==n-1 || a[i]<=a[i+1])
            a[i]=max(a[i],b[y]-a[i]);
        else a[i]=b[y]-a[i];

    }

// cout<<endl;
    int p=0;
    fi(0,n-1){
        // cout<<a[i]<<" ";
        if(a[i]>a[i+1]) p=1;
    }
    if(p==0) cout<<"Yes"<<endl;
    else cout<<"No"<<endl;


}

signed main() {
    ios_base::sync_with_stdio(false);
    cin.tie(nullptr);


    int x=0;
    int t=1;
    cin >> t;
    while (x<t) {
        solve();
        x++;
    }

    return 0;
}

// all means v.begin(),v.end()
// rall means v.rbegin(),v.rend()
// from (v,a,b) means 0 based from a to b for example 0 to v.size()-1 means complete vector
// vi means vector<int>
// vvi means vector<vector<int>>
// vs means vector<string>
// vvs means vector<vector<string>>
// f means for(int i=0;i<n;i++)
// fi(a,b) means for(int i=a;i<b;i++)
/// fact[i] means smallest factor of i
// dfs need root graph and visited array 
// bfs need root and graph




// int n;
//     cin>>n;
//     string s;
//     cin>>s;
//     int ans=INT_MAX;
//     // cout<<s<<" ";
//     vi v(26);
//     fi(0,n){
//         // cout<<s[i]<<" ";
//         v[s[i]-'a']+=1;
//     }
//     // now i will try to reduce my answer by changing mode

//     // lets start by max freq 1
//     // cout<<v[0]<<" ";

//         int len=0;
//     for(int i=1;i<=n;i++){
//         if(n%i!=0) continue;
//         // now lets try to calculate answer for this particular frequency
//         int curr=0;
//         vector<int>extra;
//         int count=0;
//         // count store number of characters in the list
//         // count should not exceed 36

//         for(int j=0;j<26;j++){
//             // if frequency is q
//             int freq=v[j];
//             if(freq==0) continue;
//             if(freq==i){
//                 // very good
//                 // count++;
//             }
//             else if(freq>i){
//                 // ye to jyada ho gya
//                 // count++;
//                 freq-=i;
//                 // now i have freq count left
//                 count+=(freq/i);
//                 freq=(freq%i);
//                 if(freq>0)
//                 extra.push_back(freq);
//             }
//             else{
//                 if(freq>0)
//                 extra.push_back(freq);
//             }

//         }

//         sort(extra.begin(),extra.end());
//         // cout<<endl;
//         // cout<<i<<" ";
//         // for(auto iii:extra) cout<<iii<<" ";
//         //     cout<<endl;
//         int m=extra.size();
//         int ii=0,jj=m-1;
//         while(ii<jj){
//             int cc=0;
//             cc+=extra[jj];
//             while(cc<i){
//                 cc+=extra[ii];
//                 if(cc<i)
//                 ii++;
//             }
//             if(cc>=i){
//                 // kaam ho gya
//                 count+=(i-extra[jj]);
//                 cc%=i;
//             }
//             if(cc>0)
//             extra[ii]=(cc%i);
//         else ii++;
//             jj--;
//         }                                         
//         // cout<<count<<" ";
//         if(count<27) {if(count<ans) {ans=min(ans,count); len=i;}}

//     }
//     cout<<len<<endl; 
//     // now make the string of length len
//     vector<pair<char,int>>vv;
//     fi(0,26){
//         if((v[i]<len || v[i]>len) && v[i]!=0){
//             if(v[i]<len)
//             vv.push_back(mp('a'+i,v[i]));
//             else vv.push_back(mp('a'+i,v[i]-len));
//         }
//     }
//     sort(vv.begin(),vv.end(),cmpsecond);
    
//     fi(0,vv.size()) cout<<vv[i].first<< " "<<vv[i].second<<endl;
//     
