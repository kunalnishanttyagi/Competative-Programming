class Solution {
public:
    int numDistinct(string s, string t) {
        // now here i can use dp
        long long n=s.size(),m=t.size();
        vector<vector<unsigned int>>dp(n,vector<unsigned int>(m,0));
        // now dp[i][j]=number of ways to make j length of string t such that we can choose j characters from string i
        for(int i=0;i<n;i++){
            for(int j=0;j<m;j++){
                if(s[i]==t[j]){
                    // so now i can add 1 to my answer
                    if(i-1>=0 && j>=1)
                    dp[i][j]=dp[i-1][j-1];
                    if(i-1>=0)
                    dp[i][j]=dp[i][j]+ dp[i-1][j];
                    if(j-1<0) dp[i][j]=dp[i][j]+1;
                    // dp[i][j];
                }
                else{
                    if(i>=1)
                    dp[i][j]=dp[i-1][j];
                    // dp[i][j]=dp[i-1][j]
                    else dp[i][j]=0;
                }
            }
        }
        return dp[n-1][m-1];
    }
};
